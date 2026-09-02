"""导出1550 nm有源波导TE0/TE1复数电磁场，供光电重叠积分使用。"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import numpy as np

import sys


ROOT = Path(__file__).resolve().parents[2]
OPTICAL_MODELS = ROOT / "models" / "optical"
sys.path.insert(0, str(OPTICAL_MODELS))

from mode_pdk_sweep import (  # noqa: E402
    DEFAULT_CONFIG,
    build_case,
    load_config,
    load_lumapi,
    mode_localization_metrics,
    scalar,
)


RESULTS = ROOT / "results" / "eo"


def read_mode(mode: Any, solver_mode: int, central_half_width_um: float) -> dict[str, Any]:
    path = f"FDE::data::mode{solver_mode}"
    x_m = np.asarray(mode.getdata(path, "x"), dtype=float).reshape(-1)
    y_m = np.asarray(mode.getdata(path, "y"), dtype=float).reshape(-1)
    fields: dict[str, np.ndarray] = {}
    for component in ("Ex", "Ey", "Ez", "Hx", "Hy", "Hz"):
        value = np.asarray(mode.getdata(path, component)).squeeze()
        if value.shape != (x_m.size, y_m.size):
            value = np.reshape(value, (x_m.size, y_m.size))
        fields[component] = value
    return {
        "solver_mode": solver_mode,
        "x_m": x_m,
        "y_m": y_m,
        "neff": float(scalar(mode.getdata(path, "neff")).real),
        "te_fraction": float(
            scalar(mode.getdata(path, "TE polarization fraction")).real
        ),
        **mode_localization_metrics(mode, path, central_half_width_um),
        **fields,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--width-um", type=float, default=1.33)
    parser.add_argument("--etch-depth-um", type=float, default=0.2)
    parser.add_argument("--wavelength-nm", type=float, default=1550.0)
    parser.add_argument("--mesh-scale", type=float, default=1.5)
    parser.add_argument("--search-neff", type=float, default=1.8)
    parser.add_argument(
        "--output",
        type=Path,
        default=RESULTS / "active_waveguide_TE0_TE1_complex_fields_1550nm.npz",
    )
    args = parser.parse_args()
    cfg = load_config(args.config)
    cfg["fde"]["mesh_cells_x"] = round(cfg["fde"]["mesh_cells_x"] * args.mesh_scale)
    cfg["fde"]["mesh_cells_y"] = round(cfg["fde"]["mesh_cells_y"] * args.mesh_scale)
    # 含金属时最高折射率解是金属表面模，必须在有源LN模式附近定向搜索。
    cfg["fde"]["search_neff"] = args.search_neff
    cfg["fde"]["trial_modes"] = max(int(cfg["fde"]["trial_modes"]), 12)
    lumapi = load_lumapi(cfg)
    mode = lumapi.MODE(hide=bool(cfg["lumerical"]["hide"]))
    try:
        build_case(
            mode,
            cfg,
            args.width_um,
            args.etch_depth_um,
            True,
            args.wavelength_nm * 1e-3,
        )
        found = int(mode.findmodes())
        candidates = [
            read_mode(
                mode,
                solver_mode,
                float(cfg["fde"]["central_half_width_um"]),
            )
            for solver_mode in range(1, found + 1)
        ]
        for item in candidates:
            print(
                f"candidate solver_mode={item['solver_mode']} "
                f"neff={item['neff']:.9f} te_fraction={item['te_fraction']:.6f} "
                f"central={item['central_energy_fraction']:.6f} "
                f"parity={item['parity_correlation']:.6f}"
            )
        candidates = [
            item
            for item in candidates
            if item["te_fraction"] >= float(cfg["fde"]["te_fraction_min"])
            and item["central_energy_fraction"]
            >= float(cfg["fde"]["central_energy_fraction_min"])
        ]
        candidates.sort(key=lambda item: item["neff"], reverse=True)
        if len(candidates) < 2:
            raise RuntimeError("未找到两个满足偏振与芯区局域条件的TE模式")
        args.output.parent.mkdir(parents=True, exist_ok=True)
        payload: dict[str, Any] = {
            "wavelength_m": np.asarray(args.wavelength_nm * 1e-9),
            "width_um": np.asarray(args.width_um),
            "etch_depth_um": np.asarray(args.etch_depth_um),
            "sidewall_angle_deg": np.asarray(
                cfg["waveguide_um"]["sidewall_angle_deg_from_horizontal"]
            ),
            "metal_included": np.asarray(True),
        }
        for mode_index, result in enumerate(candidates[:2]):
            prefix = f"te{mode_index}_"
            for key, value in result.items():
                payload[prefix + key] = np.asarray(value)
        np.savez_compressed(args.output, **payload)
        print(f"output={args.output.resolve()}")
        for mode_index, result in enumerate(candidates[:2]):
            print(
                f"mode=TE{mode_index} solver_mode={result['solver_mode']} "
                f"neff={result['neff']:.9f} te_fraction={result['te_fraction']:.6f} "
                f"central={result['central_energy_fraction']:.6f}"
            )
    finally:
        mode.close()


if __name__ == "__main__":
    main()
