"""导出70°有源LN脊波导在1550 nm的TE0/TE1模式分布。"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm
from matplotlib.patches import Polygon
import numpy as np

from mode_pdk_sweep import (
    DEFAULT_CONFIG,
    RESULTS,
    build_case,
    load_config,
    load_lumapi,
    mode_localization_metrics,
    rib_widths_um,
    scalar,
)


def read_mode(mode: Any, solver_mode: int, central_half_width_um: float) -> dict[str, Any]:
    path = f"FDE::data::mode{solver_mode}"
    x_um = np.asarray(mode.getdata(path, "x"), dtype=float).reshape(-1) * 1e6
    y_um = np.asarray(mode.getdata(path, "y"), dtype=float).reshape(-1) * 1e6
    fields = [
        np.asarray(mode.getdata(path, component)).squeeze()
        for component in ("Ex", "Ey", "Ez")
    ]
    intensity = sum(np.abs(field) ** 2 for field in fields)
    if intensity.shape != (x_um.size, y_um.size):
        intensity = np.reshape(intensity, (x_um.size, y_um.size))
    intensity /= float(np.max(intensity))
    return {
        "x_um": x_um,
        "y_um": y_um,
        "intensity": intensity,
        "neff": float(scalar(mode.getdata(path, "neff")).real),
        "te_fraction": float(
            scalar(mode.getdata(path, "TE polarization fraction")).real
        ),
        **mode_localization_metrics(mode, path, central_half_width_um),
    }


def ln_polygons(cfg: dict[str, Any], width_um: float, etch_um: float) -> list[np.ndarray]:
    stack = cfg["stack_um"]
    waveguide = cfg["waveguide_um"]
    angle_rad = math.radians(
        float(waveguide["sidewall_angle_deg_from_horizontal"])
    )
    ln_bottom = float(stack["sin"]) + float(stack["ln_sin_interlayer_oxide"])
    slab_height = float(stack["ln"]) - etch_um
    slab_top = ln_bottom + slab_height
    isolation_top = float(waveguide["ln_isolation_width"])
    isolation_bottom = isolation_top + 2.0 * slab_height / math.tan(angle_rad)
    _, rib_bottom = rib_widths_um(width_um, etch_um, math.degrees(angle_rad))
    return [
        np.asarray(
            [
                [-0.5 * isolation_bottom, ln_bottom],
                [0.5 * isolation_bottom, ln_bottom],
                [0.5 * isolation_top, slab_top],
                [-0.5 * isolation_top, slab_top],
            ]
        ),
        np.asarray(
            [
                [-0.5 * rib_bottom, slab_top],
                [0.5 * rib_bottom, slab_top],
                [0.5 * width_um, ln_bottom + float(stack["ln"])],
                [-0.5 * width_um, ln_bottom + float(stack["ln"])],
            ]
        ),
    ]


def plot_modes(
    cfg: dict[str, Any],
    width_um: float,
    etch_um: float,
    wavelength_nm: float,
    modes: list[dict[str, Any]],
    output: Path,
) -> None:
    plt.rcParams.update(
        {
            "font.sans-serif": ["Microsoft YaHei", "SimHei", "DejaVu Sans"],
            "axes.unicode_minus": False,
            "font.size": 10,
        }
    )
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 3.8), sharex=True, sharey=True)
    image = None
    for index, (ax, result) in enumerate(zip(axes, modes)):
        image = ax.pcolormesh(
            result["x_um"],
            result["y_um"],
            result["intensity"].T,
            shading="auto",
            cmap="magma",
            norm=PowerNorm(gamma=0.45, vmin=0.0, vmax=1.0),
            rasterized=True,
        )
        for vertices in ln_polygons(cfg, width_um, etch_um):
            ax.add_patch(
                Polygon(
                    vertices,
                    closed=True,
                    facecolor="none",
                    edgecolor="#2de2e6",
                    linewidth=1.25,
                )
            )
        ax.set_title(
            f"TE{index}：有效折射率={result['neff']:.5f}\n"
            f"TE占比={100.0 * result['te_fraction']:.2f}%，"
            f"中心能量={100.0 * result['central_energy_fraction']:.2f}%"
        )
        ax.set_xlim(-2.5, 2.5)
        ax.set_ylim(0.55, 1.35)
        ax.set_xlabel("横向位置（µm）")
    axes[0].set_ylabel("竖直位置（µm）")
    if image is not None:
        colorbar = fig.colorbar(image, ax=axes.ravel().tolist(), fraction=0.028, pad=0.025)
        colorbar.set_label("归一化光强 |E|²")
    angle = float(cfg["waveguide_um"]["sidewall_angle_deg_from_horizontal"])
    fig.suptitle(
        f"{wavelength_nm:g} nm有源LN波导模式（顶宽{width_um:g} µm，侧壁与水平面{angle:g}°）",
        fontsize=13,
    )
    fig.subplots_adjust(left=0.08, right=0.91, bottom=0.16, top=0.76, wspace=0.12)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=240, facecolor="white")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--width-um", type=float, default=1.33)
    parser.add_argument("--etch-depth-um", type=float, default=0.2)
    parser.add_argument("--wavelength-nm", type=float, default=1550.0)
    parser.add_argument("--mesh-scale", type=float, default=1.5)
    parser.add_argument(
        "--output",
        type=Path,
        default=RESULTS / "有源波导_1550nm_TE0_TE1模式分布_70度.png",
    )
    args = parser.parse_args()
    cfg = load_config(args.config)
    if args.mesh_scale <= 0:
        raise ValueError("网格缩放系数必须大于零")
    cfg["fde"]["mesh_cells_x"] = round(
        cfg["fde"]["mesh_cells_x"] * args.mesh_scale
    )
    cfg["fde"]["mesh_cells_y"] = round(
        cfg["fde"]["mesh_cells_y"] * args.mesh_scale
    )
    wavelength_um = args.wavelength_nm * 1e-3
    lumapi = load_lumapi(cfg)
    mode = lumapi.MODE(hide=bool(cfg["lumerical"]["hide"]))
    try:
        build_case(
            mode,
            cfg,
            args.width_um,
            args.etch_depth_um,
            False,
            wavelength_um,
        )
        found = int(mode.findmodes())
        if found < 2:
            raise RuntimeError("未找到TE0和TE1两个有源波导模式")
        candidates = [
            read_mode(
                mode,
                solver_mode,
                float(cfg["fde"]["central_half_width_um"]),
            )
            for solver_mode in range(1, found + 1)
        ]
        candidates = [
            result
            for result in candidates
            if result["te_fraction"] >= float(cfg["fde"]["te_fraction_min"])
            and result["central_energy_fraction"]
            >= float(cfg["fde"]["central_energy_fraction_min"])
        ]
        candidates.sort(key=lambda result: result["neff"], reverse=True)
        if len(candidates) < 2:
            raise RuntimeError("未找到两个满足偏振和芯区局域条件的TE模式")
        modes = candidates[:2]
        plot_modes(
            cfg,
            args.width_um,
            args.etch_depth_um,
            args.wavelength_nm,
            modes,
            args.output,
        )
        print(f"output={args.output}")
        for index, result in enumerate(modes):
            print(
                f"mode=TE{index} neff={result['neff']:.6f} "
                f"te_fraction={result['te_fraction']:.6f} "
                f"central={result['central_energy_fraction']:.6f}"
            )
    finally:
        mode.close()


if __name__ == "__main__":
    main()
