"""扫描本PDK绝热TE0/TE1模式复用器的局部耦合截面。

两条LN1脊波导位于同一个LN2保留平台上。主支路逐渐变窄，辅助支路
逐渐变宽。本脚本沿归一化传播位置逐点求解局部本征模，用于寻找
TE1(主支路)-TE0(辅助支路)的避免交叉；完整器件传播效率随后由EME求解。
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Any

import numpy as np

from mode_pdk_sweep import (
    DEFAULT_CONFIG,
    RESULTS,
    add_polygon,
    add_rect,
    ln_index_string,
    load_config,
    load_lumapi,
    rib_widths_um,
    scalar,
)


def interpolate(start: float, stop: float, position: float) -> float:
    return start + (stop - start) * position


def add_rib(
    mode: Any,
    name: str,
    center_um: float,
    top_width_um: float,
    slab_top_um: float,
    ln_top_um: float,
    etch_depth_um: float,
    sidewall_angle_deg: float,
    wavelength_um: float,
) -> None:
    _, bottom_width_um = rib_widths_um(
        top_width_um, etch_depth_um, sidewall_angle_deg
    )
    add_polygon(
        mode,
        name,
        np.asarray(
            [
                [center_um - 0.5 * bottom_width_um, slab_top_um],
                [center_um + 0.5 * bottom_width_um, slab_top_um],
                [center_um + 0.5 * top_width_um, ln_top_um],
                [center_um - 0.5 * top_width_um, ln_top_um],
            ]
        ),
        "<Object defined dielectric>",
        1,
        ln_index_string(wavelength_um),
    )


def build_coupled_section(
    mode: Any,
    cfg: dict[str, Any],
    main_width_um: float,
    auxiliary_width_um: float,
    top_gap_um: float,
    wavelength_um: float,
) -> tuple[float, float]:
    mode.switchtolayout()
    mode.deleteall()

    stack = cfg["stack_um"]
    waveguide = cfg["waveguide_um"]
    materials = cfg["material_models"]
    fde = cfg["fde"]
    etch_depth_um = float(waveguide["etch_depths"][0])
    sidewall_angle_deg = float(waveguide["sidewall_angle_deg_from_horizontal"])
    half_span_um = 0.5 * float(fde["x_span_um"])

    sin_top_um = float(stack["sin"])
    ln_bottom_um = sin_top_um + float(stack["ln_sin_interlayer_oxide"])
    ln_top_um = ln_bottom_um + float(stack["ln"])
    slab_height_um = float(stack["ln"]) - etch_depth_um
    slab_top_um = ln_bottom_um + slab_height_um
    top_oxide_um = ln_top_um + float(stack["top_cladding"])

    # 把版图顶边间隙放在x=0；左侧为主支路，右侧为辅助支路。
    main_center_um = -0.5 * top_gap_um - 0.5 * main_width_um
    auxiliary_center_um = 0.5 * top_gap_um + 0.5 * auxiliary_width_um

    add_rect(
        mode,
        "Bottom_Oxide",
        -half_span_um,
        half_span_um,
        float(fde["y_min_um"]),
        0.0,
        materials["sio2"],
        4,
    )
    add_rect(
        mode,
        "SiN_Removed_Oxide_Fill",
        -half_span_um,
        half_span_um,
        0.0,
        sin_top_um,
        materials["sio2"],
        3,
    )
    add_rect(
        mode,
        "LN_SiN_Interlayer_Oxide",
        -half_span_um,
        half_span_um,
        sin_top_um,
        ln_bottom_um,
        materials["sio2"],
        4,
    )

    ln2_top_width_um = float(waveguide["ln_isolation_width"])
    _, ln2_bottom_width_um = rib_widths_um(
        ln2_top_width_um, slab_height_um, sidewall_angle_deg
    )
    add_polygon(
        mode,
        "LN_Residual_Slab",
        np.asarray(
            [
                [-0.5 * ln2_bottom_width_um, ln_bottom_um],
                [0.5 * ln2_bottom_width_um, ln_bottom_um],
                [0.5 * ln2_top_width_um, slab_top_um],
                [-0.5 * ln2_top_width_um, slab_top_um],
            ]
        ),
        "<Object defined dielectric>",
        1,
        ln_index_string(wavelength_um),
    )
    add_rib(
        mode,
        "LN1_Main_Rib",
        main_center_um,
        main_width_um,
        slab_top_um,
        ln_top_um,
        etch_depth_um,
        sidewall_angle_deg,
        wavelength_um,
    )
    add_rib(
        mode,
        "LN1_Auxiliary_Rib",
        auxiliary_center_um,
        auxiliary_width_um,
        slab_top_um,
        ln_top_um,
        etch_depth_um,
        sidewall_angle_deg,
        wavelength_um,
    )
    add_rect(
        mode,
        "Top_Oxide",
        -half_span_um,
        half_span_um,
        ln_bottom_um,
        top_oxide_um,
        materials["sio2"],
        4,
    )

    mode.addfde()
    mode.set("solver type", "2D Z normal")
    mode.set("x min", -half_span_um * 1e-6)
    mode.set("x max", half_span_um * 1e-6)
    mode.set("y min", float(fde["y_min_um"]) * 1e-6)
    mode.set("y max", float(fde["y_max_um"]) * 1e-6)
    mode.set("x min bc", "PML")
    mode.set("x max bc", "PML")
    mode.set("y min bc", "PML")
    mode.set("y max bc", "PML")
    mode.set("mesh cells x", int(fde["mesh_cells_x"]))
    mode.set("mesh cells y", int(fde["mesh_cells_y"]))
    mode.set("wavelength", wavelength_um * 1e-6)
    mode.set("number of trial modes", int(fde["trial_modes"]))
    mode.set("calculate group index", False)
    mode.set("detailed dispersion calculation", False)
    mode.set("simulation temperature", float(cfg["temperature_c"]))
    return main_center_um, auxiliary_center_um


def branch_metrics(
    mode: Any,
    path: str,
    main_center_um: float,
    auxiliary_center_um: float,
    main_width_um: float,
    auxiliary_width_um: float,
) -> dict[str, float]:
    x_um = np.asarray(mode.getdata(path, "x"), dtype=float).reshape(-1) * 1e6
    y_um = np.asarray(mode.getdata(path, "y"), dtype=float).reshape(-1) * 1e6
    fields = [
        np.asarray(mode.getdata(path, component)).squeeze()
        for component in ("Ex", "Ey", "Ez")
    ]
    intensity = sum(np.abs(field) ** 2 for field in fields)
    if intensity.shape != (x_um.size, y_um.size):
        intensity = np.reshape(intensity, (x_um.size, y_um.size))
    weights = np.abs(np.gradient(x_um))[:, None] * np.abs(np.gradient(y_um))[None, :]
    lateral_energy = np.sum(intensity * weights, axis=1)
    total = float(np.sum(lateral_energy))

    # 窗口比脊顶各向外扩展0.45 µm，用来包含倏逝场但避免两支路窗口重叠。
    padding_um = 0.45
    main_window = np.abs(x_um - main_center_um) <= 0.5 * main_width_um + padding_um
    auxiliary_window = (
        np.abs(x_um - auxiliary_center_um)
        <= 0.5 * auxiliary_width_um + padding_um
    )
    pair_left_um = main_center_um - 0.5 * main_width_um - 0.7
    pair_right_um = auxiliary_center_um + 0.5 * auxiliary_width_um + 0.7
    pair_window = (x_um >= pair_left_um) & (x_um <= pair_right_um)
    return {
        "main_energy_fraction": float(np.sum(lateral_energy[main_window]) / total),
        "auxiliary_energy_fraction": float(
            np.sum(lateral_energy[auxiliary_window]) / total
        ),
        "pair_energy_fraction": float(np.sum(lateral_energy[pair_window]) / total),
        "peak_x_um": float(x_um[int(np.argmax(lateral_energy))]),
        "centroid_x_um": float(np.sum(x_um * lateral_energy) / total),
    }


def solve_section(
    mode: Any,
    cfg: dict[str, Any],
    position: float,
    main_width_um: float,
    auxiliary_width_um: float,
    top_gap_um: float,
    wavelength_um: float,
) -> list[dict[str, float | int]]:
    main_center_um, auxiliary_center_um = build_coupled_section(
        mode,
        cfg,
        main_width_um,
        auxiliary_width_um,
        top_gap_um,
        wavelength_um,
    )
    found = int(mode.findmodes())
    rows: list[dict[str, float | int]] = []
    for mode_index in range(1, found + 1):
        path = f"FDE::data::mode{mode_index}"
        te_fraction = float(
            scalar(mode.getdata(path, "TE polarization fraction")).real
        )
        if te_fraction < 0.5:
            continue
        neff = scalar(mode.getdata(path, "neff"))
        metrics = branch_metrics(
            mode,
            path,
            main_center_um,
            auxiliary_center_um,
            main_width_um,
            auxiliary_width_um,
        )
        rows.append(
            {
                "position": position,
                "main_width_um": main_width_um,
                "auxiliary_width_um": auxiliary_width_um,
                "top_gap_um": top_gap_um,
                "solver_mode": mode_index,
                "neff_real": float(neff.real),
                "neff_imag": float(neff.imag),
                "te_fraction": te_fraction,
                **metrics,
            }
        )
    rows.sort(key=lambda row: float(row["neff_real"]), reverse=True)
    for rank, row in enumerate(rows, start=1):
        row["neff_rank"] = rank
    return rows[:6]


def write_rows(path: Path, rows: list[dict[str, float | int]]) -> None:
    fields = [
        "position",
        "main_width_um",
        "auxiliary_width_um",
        "top_gap_um",
        "neff_rank",
        "solver_mode",
        "neff_real",
        "neff_imag",
        "te_fraction",
        "main_energy_fraction",
        "auxiliary_energy_fraction",
        "pair_energy_fraction",
        "peak_x_um",
        "centroid_x_um",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument(
        "--output",
        type=Path,
        default=RESULTS / "模式复用器局部截面扫描.csv",
    )
    parser.add_argument("--main-start-um", type=float, default=1.60)
    parser.add_argument("--main-stop-um", type=float, default=1.10)
    parser.add_argument("--aux-start-um", type=float, default=0.30)
    parser.add_argument("--aux-stop-um", type=float, default=0.80)
    parser.add_argument("--gap-um", type=float, default=0.50)
    parser.add_argument("--points", type=int, default=31)
    parser.add_argument("--scan-start", type=float, default=0.0)
    parser.add_argument("--scan-stop", type=float, default=1.0)
    parser.add_argument("--wavelength-nm", type=float, default=1550.0)
    parser.add_argument("--mesh-scale", type=float, default=1.0)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.points < 3:
        raise ValueError("扫描点数必须至少为3")
    if not 0.0 <= args.scan_start < args.scan_stop <= 1.0:
        raise ValueError("归一化扫描范围必须满足0≤起点<终点≤1")
    if min(args.main_start_um, args.main_stop_um, args.aux_start_um, args.aux_stop_um) < 0.3:
        raise ValueError("波导顶宽不能小于PDK规定的0.3 µm")
    if args.gap_um < 0.3:
        raise ValueError("波导顶边间隙不能小于PDK规定的0.3 µm")
    if args.mesh_scale <= 0:
        raise ValueError("网格缩放系数必须大于零")

    cfg = load_config(args.config)
    cfg["fde"]["mesh_cells_x"] = round(
        float(cfg["fde"]["mesh_cells_x"]) * args.mesh_scale
    )
    cfg["fde"]["mesh_cells_y"] = round(
        float(cfg["fde"]["mesh_cells_y"]) * args.mesh_scale
    )
    cfg["fde"]["trial_modes"] = max(10, int(cfg["fde"]["trial_modes"]))
    wavelength_um = args.wavelength_nm * 1e-3

    lumapi = load_lumapi(cfg)
    mode = lumapi.MODE(hide=bool(cfg["lumerical"]["hide"]))
    all_rows: list[dict[str, float | int]] = []
    try:
        print(f"mode_version={mode.version()}", flush=True)
        for position in np.linspace(args.scan_start, args.scan_stop, args.points):
            main_width_um = interpolate(
                args.main_start_um, args.main_stop_um, float(position)
            )
            auxiliary_width_um = interpolate(
                args.aux_start_um, args.aux_stop_um, float(position)
            )
            rows = solve_section(
                mode,
                cfg,
                float(position),
                main_width_um,
                auxiliary_width_um,
                args.gap_um,
                wavelength_um,
            )
            all_rows.extend(rows)
            top_three = [float(row["neff_real"]) for row in rows[:3]]
            gap_23 = top_three[1] - top_three[2] if len(top_three) >= 3 else float("nan")
            print(
                f"position={position:.4f} main={main_width_um:.4f} "
                f"aux={auxiliary_width_um:.4f} neff_top3={top_three} "
                f"delta_n23={gap_23:.6f}",
                flush=True,
            )
        write_rows(args.output, all_rows)
        rank_2 = {
            float(row["position"]): float(row["neff_real"])
            for row in all_rows
            if int(row["neff_rank"]) == 2
        }
        rank_3 = {
            float(row["position"]): float(row["neff_real"])
            for row in all_rows
            if int(row["neff_rank"]) == 3
        }
        common = sorted(set(rank_2) & set(rank_3))
        if common:
            crossing_position = min(common, key=lambda p: rank_2[p] - rank_3[p])
            print(
                f"minimum_delta_n23={rank_2[crossing_position] - rank_3[crossing_position]:.8f} "
                f"at_position={crossing_position:.4f}",
                flush=True,
            )
        print(f"output={args.output}", flush=True)
    finally:
        mode.close()


if __name__ == "__main__":
    main()
