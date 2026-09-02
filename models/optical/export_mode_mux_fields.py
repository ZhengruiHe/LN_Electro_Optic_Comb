"""导出1550 nm模式复用器关键截面的归一化光强分布。

图中上排跟踪主波导TE0直通支，下排跟踪主波导TE1向
辅助波导TE0转换的绝热本征模分支。LN1和LN2侧壁角均从
PDK配置读取，定义为与水平面的夹角。
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm
from matplotlib.patches import Polygon
import numpy as np

from mode_mux_cross_section_sweep import build_coupled_section, branch_metrics
from mode_pdk_sweep import DEFAULT_CONFIG, RESULTS, load_config, load_lumapi, scalar


def trapezoid_vertices(
    center_um: float,
    top_width_um: float,
    z_bottom_um: float,
    height_um: float,
    angle_deg_from_horizontal: float,
) -> np.ndarray:
    bottom_width_um = top_width_um + 2.0 * height_um / math.tan(
        math.radians(angle_deg_from_horizontal)
    )
    return np.asarray(
        [
            [center_um - 0.5 * bottom_width_um, z_bottom_um],
            [center_um + 0.5 * bottom_width_um, z_bottom_um],
            [center_um + 0.5 * top_width_um, z_bottom_um + height_um],
            [center_um - 0.5 * top_width_um, z_bottom_um + height_um],
        ]
    )


def section_geometry(
    cfg: dict[str, Any],
    main_width_um: float,
    auxiliary_width_um: float,
    gap_um: float,
) -> list[tuple[str, np.ndarray]]:
    stack = cfg["stack_um"]
    waveguide = cfg["waveguide_um"]
    angle_deg = float(waveguide["sidewall_angle_deg_from_horizontal"])
    etch_depth_um = float(waveguide["etch_depths"][0])
    ln_bottom_um = float(stack["sin"]) + float(stack["ln_sin_interlayer_oxide"])
    slab_height_um = float(stack["ln"]) - etch_depth_um
    slab_top_um = ln_bottom_um + slab_height_um
    main_center_um = -0.5 * gap_um - 0.5 * main_width_um
    auxiliary_center_um = 0.5 * gap_um + 0.5 * auxiliary_width_um
    return [
        (
            "LN2保留平台",
            trapezoid_vertices(
                0.0,
                float(waveguide["ln_isolation_width"]),
                ln_bottom_um,
                slab_height_um,
                angle_deg,
            ),
        ),
        (
            "LN1主脊",
            trapezoid_vertices(
                main_center_um,
                main_width_um,
                slab_top_um,
                etch_depth_um,
                angle_deg,
            ),
        ),
        (
            "LN1辅助脊",
            trapezoid_vertices(
                auxiliary_center_um,
                auxiliary_width_um,
                slab_top_um,
                etch_depth_um,
                angle_deg,
            ),
        ),
    ]


def read_mode(
    mode: Any,
    solver_mode: int,
    main_center_um: float,
    auxiliary_center_um: float,
    main_width_um: float,
    auxiliary_width_um: float,
) -> dict[str, Any]:
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
    maximum = float(np.max(intensity))
    if maximum <= 0.0:
        raise RuntimeError(f"模式{solver_mode}光强为零")
    metrics = branch_metrics(
        mode,
        path,
        main_center_um,
        auxiliary_center_um,
        main_width_um,
        auxiliary_width_um,
    )
    return {
        "x_um": x_um,
        "y_um": y_um,
        "intensity": intensity / maximum,
        "neff": float(scalar(mode.getdata(path, "neff")).real),
        **metrics,
    }


def solve_section(
    mode: Any,
    cfg: dict[str, Any],
    name: str,
    main_width_um: float,
    auxiliary_width_um: float,
    gap_um: float,
    wavelength_um: float,
) -> dict[str, Any]:
    main_center_um, auxiliary_center_um = build_coupled_section(
        mode,
        cfg,
        main_width_um,
        auxiliary_width_um,
        gap_um,
        wavelength_um,
    )
    found = int(mode.findmodes())
    if found < 3:
        raise RuntimeError(f"{name}未找到所需的前3个模式")
    return {
        "name": name,
        "main_width_um": main_width_um,
        "auxiliary_width_um": auxiliary_width_um,
        "gap_um": gap_um,
        "geometry": section_geometry(
            cfg, main_width_um, auxiliary_width_um, gap_um
        ),
        "direct": read_mode(
            mode,
            1,
            main_center_um,
            auxiliary_center_um,
            main_width_um,
            auxiliary_width_um,
        ),
        "conversion": read_mode(
            mode,
            3,
            main_center_um,
            auxiliary_center_um,
            main_width_um,
            auxiliary_width_um,
        ),
    }


def plot_results(
    cfg: dict[str, Any], sections: list[dict[str, Any]], output: Path
) -> None:
    plt.rcParams.update(
        {
            "font.sans-serif": ["Microsoft YaHei", "SimHei", "DejaVu Sans"],
            "axes.unicode_minus": False,
            "font.size": 9,
        }
    )
    fig, axes = plt.subplots(2, 3, figsize=(12.0, 5.7), sharex=True, sharey=True)
    row_keys = ("direct", "conversion")
    row_names = ("TE0直通分支", "TE1→TE0转换分支")
    image = None
    for column, section in enumerate(sections):
        for row, (key, row_name) in enumerate(zip(row_keys, row_names)):
            ax = axes[row, column]
            result = section[key]
            image = ax.pcolormesh(
                result["x_um"],
                result["y_um"],
                result["intensity"].T,
                shading="auto",
                cmap="magma",
                norm=PowerNorm(gamma=0.45, vmin=0.0, vmax=1.0),
                rasterized=True,
            )
            for geometry_name, vertices in section["geometry"]:
                ax.add_patch(
                    Polygon(
                        vertices,
                        closed=True,
                        facecolor="none",
                        edgecolor="#2de2e6",
                        linewidth=1.15,
                    )
                )
            main_fraction = 100.0 * result["main_energy_fraction"]
            auxiliary_fraction = 100.0 * result["auxiliary_energy_fraction"]
            ax.set_title(
                f"{section['name']}\n"
                f"有效折射率={result['neff']:.4f}，主路{main_fraction:.1f}%，"
                f"辅助{auxiliary_fraction:.1f}%",
                fontsize=9,
            )
            ax.set_xlim(-3.9, 3.9)
            ax.set_ylim(0.45, 1.45)
            ax.set_aspect("auto")
            if column == 0:
                ax.set_ylabel(f"{row_name}\n竖直方向（µm）")
            if row == 1:
                ax.set_xlabel("横向（µm）")
    if image is not None:
        colorbar = fig.colorbar(image, ax=axes.ravel().tolist(), fraction=0.018, pad=0.015)
        colorbar.set_label("归一化光强 |E|²")
    angle_deg = float(cfg["waveguide_um"]["sidewall_angle_deg_from_horizontal"])
    fig.suptitle(
        f"1550 nm模式复用器局部本征模分布（LN1/LN2侧壁与水平面{angle_deg:g}°）",
        fontsize=13,
    )
    fig.subplots_adjust(left=0.085, right=0.91, bottom=0.10, top=0.84, wspace=0.18, hspace=0.38)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=240, facecolor="white")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--wavelength-nm", type=float, default=1550.0)
    parser.add_argument(
        "--output",
        type=Path,
        default=RESULTS / "模式复用器_1550nm_光学模式分布.png",
    )
    args = parser.parse_args()
    cfg = load_config(args.config)
    wavelength_um = args.wavelength_nm * 1e-3
    lumapi = load_lumapi(cfg)
    mode = lumapi.MODE(hide=bool(cfg["lumerical"]["hide"]))
    try:
        sections = [
            solve_section(mode, cfg, "输入端", 1.40, 0.30, 2.50, wavelength_um),
            solve_section(mode, cfg, "避免交叉中心", 1.31, 0.33, 1.00, wavelength_um),
            solve_section(mode, cfg, "输出端", 1.10, 0.40, 2.50, wavelength_um),
        ]
        plot_results(cfg, sections, args.output)
        print(f"output={args.output}")
        for section in sections:
            for key in ("direct", "conversion"):
                result = section[key]
                print(
                    f"section={section['name']} branch={key} "
                    f"neff={result['neff']:.6f} "
                    f"main={result['main_energy_fraction']:.6f} "
                    f"auxiliary={result['auxiliary_energy_fraction']:.6f}"
                )
    finally:
        mode.close()


if __name__ == "__main__":
    main()
