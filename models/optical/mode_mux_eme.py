"""建立并求解本PDK的绝热TE0/TE1模式复用器EME模型。

EME的传播轴为x，横向为y，竖直方向为z。LN1和LN2的侧壁角
均定义为与水平面的夹角（当前配置为70°），并在竖直方向用多个薄层
逼近两级梯形。该脚本先用于粗网格验证拓扑，
最终尺寸必须再经过模式数、网格、波长和工艺容差收敛。
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any, Callable

import numpy as np

from mode_pdk_sweep import (
    DEFAULT_CONFIG,
    RESULTS,
    load_config,
    load_lumapi,
    zelmon_ln_indices,
)


def ln_index_string_eme(wavelength_um: float) -> str:
    """返回EME坐标x(传播)、y(晶体z)、z(竖直)下的LN张量。"""
    ordinary, extraordinary = zelmon_ln_indices(np.asarray([wavelength_um]))
    return f"{ordinary[0]:.12g};{extraordinary[0]:.12g};{ordinary[0]:.12g}"


def smooth_asymmetric_progress(u: np.ndarray, exponent: float, ratio: float) -> np.ndarray:
    """端点一阶平滑，并把目标避免交叉推到器件中部。"""
    numerator = np.power(u, exponent)
    denominator = numerator + ratio * np.power(1.0 - u, exponent)
    progress = np.divide(
        numerator,
        denominator,
        out=np.zeros_like(numerator),
        where=denominator > 0,
    )
    progress[-1] = 1.0
    return progress


def smoothstep(value: np.ndarray) -> np.ndarray:
    clipped = np.clip(value, 0.0, 1.0)
    return clipped * clipped * (3.0 - 2.0 * clipped)


def slow_crossing_progress(u: np.ndarray) -> np.ndarray:
    """在预计避免交叉的进度0.05--0.10附近降低几何变化速度。"""
    progress = np.empty_like(u)
    first = u <= 0.35
    middle = (u > 0.35) & (u <= 0.65)
    last = u > 0.65
    progress[first] = 0.05 * smoothstep(u[first] / 0.35)
    progress[middle] = 0.05 + 0.05 * smoothstep((u[middle] - 0.35) / 0.30)
    progress[last] = 0.10 + 0.90 * smoothstep((u[last] - 0.65) / 0.35)
    return progress


def plateau_gap_profile(
    u: np.ndarray, end_gap_um: float, minimum_gap_um: float
) -> np.ndarray:
    """两端平滑解耦，在器件中部保持恒定最小间隙。"""
    gap = np.full_like(u, minimum_gap_um)
    approach = u < 0.30
    separation = u > 0.70
    gap[approach] = end_gap_um - (end_gap_um - minimum_gap_um) * smoothstep(
        u[approach] / 0.30
    )
    gap[separation] = minimum_gap_um + (
        end_gap_um - minimum_gap_um
    ) * smoothstep((u[separation] - 0.70) / 0.30)
    return gap


def add_top_view_polygon(
    mode: Any,
    name: str,
    x_um: np.ndarray,
    center_um: np.ndarray,
    width_um: np.ndarray,
    z_min_um: float,
    z_max_um: float,
    index: str,
) -> None:
    left = np.column_stack((x_um, center_um - 0.5 * width_um))
    right = np.column_stack((x_um[::-1], (center_um + 0.5 * width_um)[::-1]))
    vertices = np.vstack((left, right)) * 1e-6
    mode.addpoly()
    mode.set("name", name)
    mode.set("vertices", vertices)
    mode.set("z min", z_min_um * 1e-6)
    mode.set("z max", z_max_um * 1e-6)
    mode.set("material", "<Object defined dielectric>")
    mode.set("index", index)
    mode.set("override mesh order from material database", True)
    mode.set("mesh order", 1)


def add_layered_trapezoid(
    mode: Any,
    name: str,
    x_um: np.ndarray,
    center_um: np.ndarray,
    top_width_um: np.ndarray,
    z_bottom_um: float,
    height_um: float,
    sidewall_angle_deg: float,
    vertical_slices: int,
    index: str,
) -> None:
    """用竖直薄层逼近沿传播方向宽度可变的梯形。

    sidewall_angle_deg为侧壁与水平面的夹角；LN1和LN2均调用此几何规则。
    """
    for slice_index in range(vertical_slices):
        z0_um = z_bottom_um + height_um * slice_index / vertical_slices
        z1_um = z_bottom_um + height_um * (slice_index + 1) / vertical_slices
        z_mid_um = 0.5 * (z0_um + z1_um)
        distance_from_top_um = z_bottom_um + height_um - z_mid_um
        width_um = top_width_um + 2.0 * distance_from_top_um / math.tan(
            math.radians(sidewall_angle_deg)
        )
        add_top_view_polygon(
            mode,
            f"{name}_slice_{slice_index + 1}",
            x_um,
            center_um,
            width_um,
            z0_um,
            z1_um,
            index,
        )


def mode_paths(
    length_um: float,
    samples: int,
    main_start_um: float,
    main_stop_um: float,
    auxiliary_start_um: float,
    auxiliary_stop_um: float,
    end_gap_um: float,
    minimum_gap_um: float,
    progress_exponent: float,
    progress_ratio: float,
    trajectory: str,
) -> dict[str, np.ndarray]:
    x_um = np.linspace(0.0, length_um, samples)
    u = x_um / length_um
    if trajectory == "slow_crossing":
        progress = slow_crossing_progress(u)
    else:
        progress = smooth_asymmetric_progress(u, progress_exponent, progress_ratio)
    main_width_um = main_start_um + (main_stop_um - main_start_um) * progress
    auxiliary_width_um = auxiliary_start_um + (
        auxiliary_stop_um - auxiliary_start_um
    ) * progress
    if trajectory == "slow_crossing":
        gap_um = plateau_gap_profile(u, end_gap_um, minimum_gap_um)
    else:
        # 余弦函数使间隙在两端为end_gap、中部为minimum_gap，且端点斜率为零。
        gap_um = minimum_gap_um + (end_gap_um - minimum_gap_um) * (
            0.5 + 0.5 * np.cos(2.0 * np.pi * u)
        )
    main_center_um = -0.5 * gap_um - 0.5 * main_width_um
    auxiliary_center_um = 0.5 * gap_um + 0.5 * auxiliary_width_um
    return {
        "x_um": x_um,
        "u": u,
        "progress": progress,
        "main_width_um": main_width_um,
        "auxiliary_width_um": auxiliary_width_um,
        "gap_um": gap_um,
        "main_center_um": main_center_um,
        "auxiliary_center_um": auxiliary_center_um,
    }


def build_model(mode: Any, cfg: dict[str, Any], args: argparse.Namespace) -> dict[str, np.ndarray]:
    mode.switchtolayout()
    mode.deleteall()
    stack = cfg["stack_um"]
    waveguide = cfg["waveguide_um"]
    materials = cfg["material_models"]
    sidewall_angle_deg = float(waveguide["sidewall_angle_deg_from_horizontal"])
    etch_depth_um = float(waveguide["etch_depths"][0])
    residual_slab_um = float(stack["ln"]) - etch_depth_um
    ln_bottom_um = float(stack["sin"]) + float(stack["ln_sin_interlayer_oxide"])
    slab_top_um = ln_bottom_um + residual_slab_um
    ln_top_um = ln_bottom_um + float(stack["ln"])
    paths = mode_paths(
        args.length_um,
        args.geometry_samples,
        args.main_start_um,
        args.main_stop_um,
        args.aux_start_um,
        args.aux_stop_um,
        args.end_gap_um,
        args.minimum_gap_um,
        args.progress_exponent,
        args.progress_ratio,
        args.trajectory,
    )
    x_um = paths["x_um"]
    index = ln_index_string_eme(args.wavelength_nm * 1e-3)

    platform_center_um = np.zeros_like(x_um)
    platform_top_width_um = np.full_like(
        x_um, float(waveguide["ln_isolation_width"])
    )
    add_layered_trapezoid(
        mode,
        "LN2_Residual_Platform",
        x_um,
        platform_center_um,
        platform_top_width_um,
        ln_bottom_um,
        residual_slab_um,
        sidewall_angle_deg,
        args.vertical_slices,
        index,
    )
    add_layered_trapezoid(
        mode,
        "LN1_Main_Rib",
        x_um,
        paths["main_center_um"],
        paths["main_width_um"],
        slab_top_um,
        etch_depth_um,
        sidewall_angle_deg,
        args.vertical_slices,
        index,
    )
    add_layered_trapezoid(
        mode,
        "LN1_Auxiliary_Rib",
        x_um,
        paths["auxiliary_center_um"],
        paths["auxiliary_width_um"],
        slab_top_um,
        etch_depth_um,
        sidewall_angle_deg,
        args.vertical_slices,
        index,
    )

    mode.addeme()
    mode.set("solver type", "3D")
    mode.set("x min", 0.0)
    mode.set("y", 0.0)
    mode.set("y span", args.y_span_um * 1e-6)
    mode.set("z min", float(cfg["fde"]["y_min_um"]) * 1e-6)
    mode.set("z max", float(cfg["fde"]["y_max_um"]) * 1e-6)
    mode.set("y min bc", "PML")
    mode.set("y max bc", "PML")
    mode.set("z min bc", "PML")
    mode.set("z max bc", "PML")
    mode.set("background material", materials["sio2"])
    mode.set("wavelength", args.wavelength_nm * 1e-9)
    mode.set("simulation temperature", float(cfg["temperature_c"]))
    mode.set("number of cell groups", 3)
    group_spans_um = np.asarray(
        [args.approach_length_um, args.coupling_length_um, args.separation_length_um]
    )
    if not math.isclose(float(np.sum(group_spans_um)), args.length_um):
        raise ValueError("三个EME分组长度之和必须等于总长度")
    mode.set("group spans", group_spans_um * 1e-6)
    mode.set("cells", np.asarray(args.cells, dtype=float))
    mode.set("subcell method", np.asarray([1.0, 1.0, 1.0]))
    mode.set("number of modes for all cell groups", args.modes)
    mode.set("mesh cells y", args.mesh_cells_y)
    mode.set("mesh cells z", args.mesh_cells_z)
    mode.set("convergence tolerance", 1e-9)
    mode.set("energy conservation", "make passive")

    port_mode_numbers = {
        "port_1": args.port1_modes,
        "port_2": args.port2_modes,
    }
    for port_name, selected_modes in port_mode_numbers.items():
        mode.select(f"EME::Ports::{port_name}")
        mode.set("use full simulation span", True)
        mode.set("mode selection", "user select")
        update_status = int(
            mode.updateportmodes(np.asarray(selected_modes, dtype=float))
        )
        if update_status != 1:
            raise RuntimeError(f"{port_name}端口模式更新失败")
    return paths


def write_geometry(path: Path, paths: dict[str, np.ndarray]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(paths)
    with path.open("w", newline="", encoding="utf-8-sig") as stream:
        writer = csv.writer(stream)
        writer.writerow(fields)
        for values in zip(*(paths[field] for field in fields)):
            writer.writerow(values)


def extract_s_matrix(result: Any) -> np.ndarray:
    if isinstance(result, dict):
        for key in ("S", "s"):
            if key in result:
                return np.asarray(result[key], dtype=complex).squeeze()
    array = np.asarray(result)
    if array.ndim >= 2:
        return np.asarray(array, dtype=complex).squeeze()
    raise RuntimeError("无法从EME结果中识别S矩阵")


def write_s_matrix(path: Path, matrix: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as stream:
        writer = csv.writer(stream)
        writer.writerow(["输出索引", "输入索引", "实部", "虚部", "功率"])
        for output_index in range(matrix.shape[0]):
            for input_index in range(matrix.shape[1]):
                value = matrix[output_index, input_index]
                writer.writerow(
                    [
                        output_index + 1,
                        input_index + 1,
                        float(np.real(value)),
                        float(np.imag(value)),
                        float(np.abs(value) ** 2),
                    ]
                )


def parse_cells(value: str) -> list[int]:
    cells = [int(item.strip()) for item in value.split(",")]
    if len(cells) != 3 or min(cells) < 1:
        raise argparse.ArgumentTypeError("cells必须是三个正整数，例如3,21,3")
    return cells


def parse_mode_numbers(value: str) -> list[int]:
    modes = [int(item.strip()) for item in value.split(",")]
    if len(modes) != 3 or min(modes) < 1 or len(set(modes)) != 3:
        raise argparse.ArgumentTypeError("端口模式必须是三个互不相同的正整数")
    return modes


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--wavelength-nm", type=float, default=1550.0)
    parser.add_argument("--length-um", type=float, default=750.0)
    parser.add_argument("--approach-length-um", type=float, default=100.0)
    parser.add_argument("--coupling-length-um", type=float, default=550.0)
    parser.add_argument("--separation-length-um", type=float, default=100.0)
    parser.add_argument("--main-start-um", type=float, default=1.40)
    parser.add_argument("--main-stop-um", type=float, default=1.10)
    parser.add_argument("--aux-start-um", type=float, default=0.30)
    parser.add_argument("--aux-stop-um", type=float, default=0.40)
    parser.add_argument("--end-gap-um", type=float, default=2.50)
    parser.add_argument("--minimum-gap-um", type=float, default=1.00)
    parser.add_argument("--progress-exponent", type=float, default=2.0)
    parser.add_argument("--progress-ratio", type=float, default=2.3333333333)
    parser.add_argument(
        "--trajectory",
        choices=("slow_crossing", "asymmetric"),
        default="asymmetric",
    )
    parser.add_argument("--geometry-samples", type=int, default=121)
    parser.add_argument("--vertical-slices", type=int, default=4)
    parser.add_argument("--y-span-um", type=float, default=12.0)
    parser.add_argument("--mesh-cells-y", type=int, default=240)
    parser.add_argument("--mesh-cells-z", type=int, default=120)
    parser.add_argument("--modes", type=int, default=8)
    parser.add_argument("--cells", type=parse_cells, default=[3, 21, 3])
    parser.add_argument("--port1-modes", type=parse_mode_numbers, default=[1, 3, 4])
    parser.add_argument("--port2-modes", type=parse_mode_numbers, default=[1, 3, 4])
    parser.add_argument("--build-only", action="store_true")
    parser.add_argument(
        "--project",
        type=Path,
        default=RESULTS / "模式复用器_EME_初始.lms",
    )
    parser.add_argument(
        "--geometry-output",
        type=Path,
        default=RESULTS / "模式复用器_EME_几何轨迹.csv",
    )
    parser.add_argument(
        "--s-output",
        type=Path,
        default=RESULTS / "模式复用器_EME_S矩阵.csv",
    )
    return parser.parse_args()


def validate(args: argparse.Namespace) -> None:
    if min(args.main_start_um, args.main_stop_um, args.aux_start_um, args.aux_stop_um) < 0.3:
        raise ValueError("LN1波导顶宽不能小于0.3 µm")
    if args.minimum_gap_um < 0.3 or args.end_gap_um < args.minimum_gap_um:
        raise ValueError("间隙必须满足PDK最小值，且端部间隙不得小于耦合区间隙")
    if args.geometry_samples < 5 or args.vertical_slices < 1:
        raise ValueError("几何采样数或竖直分层数无效")
    if args.modes < 4:
        raise ValueError("EME每个截面至少需要4个模式")


def main() -> None:
    args = parse_args()
    validate(args)
    cfg = load_config(args.config)
    lumapi = load_lumapi(cfg)
    mode = lumapi.MODE(hide=bool(cfg["lumerical"]["hide"]))
    try:
        print(f"mode_version={mode.version()}", flush=True)
        paths = build_model(mode, cfg, args)
        write_geometry(args.geometry_output, paths)
        args.project.parent.mkdir(parents=True, exist_ok=True)
        mode.save(str(args.project.resolve()))
        print(f"project={args.project}", flush=True)
        print(f"geometry={args.geometry_output}", flush=True)
        if args.build_only:
            print("status=build_only", flush=True)
            return

        print("eme_step=calculate_modes", flush=True)
        mode.run()
        print("eme_step=propagate", flush=True)
        mode.emepropagate()
        result_names = str(mode.getresult("EME"))
        print(f"eme_results={result_names}", flush=True)
        result: Any | None = None
        for result_name in ("power normalized user s matrix", "user s matrix"):
            try:
                result = mode.getresult("EME", result_name)
                print(f"eme_s_result={result_name}", flush=True)
                break
            except Exception:
                continue
        if result is None:
            raise RuntimeError("EME没有返回功率归一化或用户S矩阵")
        matrix = extract_s_matrix(result)
        write_s_matrix(args.s_output, matrix)
        print(f"s_matrix_shape={matrix.shape}", flush=True)
        print(f"s_output={args.s_output}", flush=True)
        if matrix.shape[0] >= 6 and matrix.shape[1] >= 6:
            print(
                "selected_powers="
                + json.dumps(
                    {
                        "port1_main_TE0_to_port2_main_TE0": float(
                            abs(matrix[3, 0]) ** 2
                        ),
                        "port1_main_TE1_to_port2_aux_TE0": float(
                            abs(matrix[4, 1]) ** 2
                        ),
                        "port1_aux_TE0_to_port2_main_TE1": float(
                            abs(matrix[5, 2]) ** 2
                        ),
                    },
                    ensure_ascii=False,
                ),
                flush=True,
            )
        mode.save(str(args.project.resolve()))
    finally:
        mode.close()


if __name__ == "__main__":
    main()
