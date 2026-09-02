"""核查行波电极与 TE0/TE1 光模的速度匹配及四程延迟容差。"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

import numpy as np


C_M_PER_S = 299_792_458.0


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def select_mode(
    rows: list[dict[str, str]], width_um: float, etch_um: float, label: str
) -> dict[str, str]:
    matches = [
        row
        for row in rows
        if math.isclose(float(row["width_um"]), width_um)
        and math.isclose(float(row["etch_depth_um"]), etch_um)
        and row["mode_label"] == label
    ]
    if len(matches) != 1:
        raise ValueError(
            f"不能唯一定位光模 {label}: width={width_um}, etch={etch_um}"
        )
    return matches[0]


def mismatch_metrics(
    frequency_hz: np.ndarray,
    n_rf: np.ndarray,
    n_group: float,
    electrode_length_m: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """返回折射率差、相位走离和相对理想同步的功率惩罚。"""
    delta_n = n_rf - n_group
    x = np.pi * frequency_hz * delta_n * electrode_length_m / C_M_PER_S
    amplitude = np.abs(np.sinc(x / np.pi))
    penalty_db = -20.0 * np.log10(np.maximum(amplitude, 1e-15))
    phase_walkoff_deg = np.degrees(2.0 * x)
    return delta_n, phase_walkoff_deg, penalty_db


def sinc_root_for_penalty(penalty_db: float) -> float:
    target = 10.0 ** (-penalty_db / 20.0)
    low, high = 0.0, math.pi
    for _ in range(100):
        middle = 0.5 * (low + high)
        value = math.sin(middle) / middle if middle else 1.0
        if value > target:
            low = middle
        else:
            high = middle
    return 0.5 * (low + high)


def maximum_length_for_penalty(
    frequency_hz: float, delta_n: float, penalty_db: float
) -> float:
    if abs(delta_n) < 1e-15:
        return math.inf
    x_limit = sinc_root_for_penalty(penalty_db)
    return x_limit * C_M_PER_S / (
        math.pi * frequency_hz * abs(delta_n)
    )


def nearest_index(values: np.ndarray, target: float) -> int:
    return int(np.argmin(np.abs(values - target)))


def fmt(value: float, digits: int = 4) -> str:
    if math.isinf(value):
        return "∞"
    return f"{value:.{digits}f}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--rf-csv",
        type=Path,
        default=Path("results/hfss/pdk_10GHz_a70_D90_h4_t10_r45_c5_w43_g5_phasecheck_500_1000.csv"),
    )
    parser.add_argument(
        "--optical-csv",
        type=Path,
        default=Path("results/optical/pdk_1550_active_waveguide_a70_selected.csv"),
    )
    parser.add_argument("--waveguide-width-um", type=float, default=1.33)
    parser.add_argument("--etch-depth-um", type=float, default=0.2)
    parser.add_argument("--target-frequency-ghz", type=float, default=10.0)
    parser.add_argument("--band-start-ghz", type=float, default=8.0)
    parser.add_argument("--band-stop-ghz", type=float, default=12.0)
    parser.add_argument("--electrode-length-cm", type=float, default=1.0)
    parser.add_argument("--signal-width-um", type=float, default=43.0)
    parser.add_argument("--inner-gap-um", type=float, default=5.0)
    parser.add_argument("--t-neck-length-um", type=float, default=4.0)
    parser.add_argument("--t-cap-width-um", type=float, default=2.0)
    parser.add_argument("--t-neck-width-um", type=float, default=10.0)
    parser.add_argument("--t-cap-length-um", type=float, default=45.0)
    parser.add_argument("--t-unit-gap-um", type=float, default=5.0)
    parser.add_argument("--trunk-gap-um", type=float, default=17.0)
    parser.add_argument("--loop-phase-tolerance-deg", type=float, default=5.0)
    parser.add_argument("--loop1-periods", type=float, default=2.0)
    parser.add_argument("--loop2-periods", type=float, default=2.5)
    parser.add_argument("--loop3-periods", type=float, default=2.0)
    parser.add_argument(
        "--rf-adaptive-converged",
        action="store_true",
        help="仅在HFSS自适应网格已按设定标准收敛时使用",
    )
    parser.add_argument(
        "--rf-port-extra-mode-cleared",
        action="store_true",
        help="仅在端口额外传播/慢衰减模式已排除时使用",
    )
    parser.add_argument(
        "--rf-port-cross-coupling-max-db",
        type=float,
        default=None,
        help="指定工作带内CPW模与奇模之间最差交叉耦合，单位dB",
    )
    parser.add_argument(
        "--output-dir", type=Path, default=Path("results/system/velocity_matching_10GHz")
    )
    args = parser.parse_args()

    if args.t_cap_length_um <= 0 or args.t_unit_gap_um < 0:
        raise ValueError("T帽长度必须为正，T单元间隙不能为负")
    duty_cycle = args.t_cap_length_um / (
        args.t_cap_length_um + args.t_unit_gap_um
    )

    rf_rows = read_csv(args.rf_csv)
    optical_rows = read_csv(args.optical_csv)
    frequency_ghz = np.asarray(
        [float(row["frequency_ghz"]) for row in rf_rows], dtype=float
    )
    frequency_hz = frequency_ghz * 1e9
    n_rf = np.asarray([float(row["n_rf"]) for row in rf_rows], dtype=float)
    loss_db_per_cm = np.asarray(
        [float(row["loss_db_per_cm"]) for row in rf_rows], dtype=float
    )
    impedance_real = np.asarray(
        [float(row["long_abcd_zc_real_ohm"]) for row in rf_rows], dtype=float
    )
    s11_db = np.asarray(
        [float(row["long_s11_db"]) for row in rf_rows], dtype=float
    )

    te0 = select_mode(
        optical_rows, args.waveguide_width_um, args.etch_depth_um, "TE0"
    )
    te1 = select_mode(
        optical_rows, args.waveguide_width_um, args.etch_depth_um, "TE1"
    )
    sidewall_angle_deg = float(
        te0.get("sidewall_angle_deg_from_horizontal", 90.0)
    )
    rib_top_width_um = float(te0.get("rib_top_width_um", args.waveguide_width_um))
    rib_bottom_width_um = float(
        te0.get("rib_bottom_width_um", args.waveguide_width_um)
    )
    mode_group_indices = {"TE0": float(te0["ng"]), "TE1": float(te1["ng"])}
    electrode_length_m = args.electrode_length_cm * 1e-2

    response_rows: list[dict[str, Any]] = []
    mode_metrics: dict[str, dict[str, np.ndarray]] = {}
    for label, n_group in mode_group_indices.items():
        delta_n, phase_deg, penalty_db = mismatch_metrics(
            frequency_hz, n_rf, n_group, electrode_length_m
        )
        mode_metrics[label] = {
            "delta_n": delta_n,
            "phase_deg": phase_deg,
            "penalty_db": penalty_db,
        }

    for index in range(len(frequency_ghz)):
        response_rows.append(
            {
                "frequency_ghz": frequency_ghz[index],
                "n_rf": n_rf[index],
                "ng_te0": mode_group_indices["TE0"],
                "ng_te1": mode_group_indices["TE1"],
                "delta_n_te0": mode_metrics["TE0"]["delta_n"][index],
                "delta_n_te1": mode_metrics["TE1"]["delta_n"][index],
                "phase_walkoff_te0_deg_per_cm": mode_metrics["TE0"][
                    "phase_deg"
                ][index],
                "phase_walkoff_te1_deg_per_cm": mode_metrics["TE1"][
                    "phase_deg"
                ][index],
                "mismatch_penalty_te0_db_per_cm": mode_metrics["TE0"][
                    "penalty_db"
                ][index],
                "mismatch_penalty_te1_db_per_cm": mode_metrics["TE1"][
                    "penalty_db"
                ][index],
                "rf_loss_db_per_cm": loss_db_per_cm[index],
                "impedance_real_ohm": impedance_real[index],
                "s11_db": s11_db[index],
            }
        )

    target_index = nearest_index(frequency_ghz, args.target_frequency_ghz)
    solved_target_ghz = float(frequency_ghz[target_index])
    solved_target_hz = solved_target_ghz * 1e9
    band_mask = (
        (frequency_ghz >= args.band_start_ghz)
        & (frequency_ghz <= args.band_stop_ghz)
    )
    if not np.any(band_mask):
        raise ValueError("射频数据不包含指定频带")

    mode_summary: dict[str, dict[str, float]] = {}
    for label, n_group in mode_group_indices.items():
        delta_target = float(mode_metrics[label]["delta_n"][target_index])
        mode_summary[label] = {
            "ng": n_group,
            "delta_n_at_target": delta_target,
            "phase_walkoff_deg_per_cm_at_target": float(
                mode_metrics[label]["phase_deg"][target_index]
            ),
            "mismatch_penalty_db_per_cm_at_target": float(
                mode_metrics[label]["penalty_db"][target_index]
            ),
            "maximum_abs_delta_n_in_band": float(
                np.max(np.abs(mode_metrics[label]["delta_n"][band_mask]))
            ),
            "maximum_penalty_db_per_cm_in_band": float(
                np.max(mode_metrics[label]["penalty_db"][band_mask])
            ),
            "one_db_mismatch_length_cm_at_target": 100.0
            * maximum_length_for_penalty(solved_target_hz, delta_target, 1.0),
        }

    loop_multipliers = [
        args.loop1_periods,
        args.loop2_periods,
        args.loop3_periods,
    ]
    loop_labels = ["回路1", "回路2", "回路3"]
    loop_ng = [
        mode_group_indices["TE0"],
        0.5 * (mode_group_indices["TE0"] + mode_group_indices["TE1"]),
        mode_group_indices["TE1"],
    ]
    target_hz = args.target_frequency_ghz * 1e9
    loop_rows: list[dict[str, Any]] = []
    for label, multiplier, n_group in zip(loop_labels, loop_multipliers, loop_ng):
        delay_s = multiplier / target_hz
        length_m = C_M_PER_S * delay_s / n_group
        permitted_delay_s = (
            args.loop_phase_tolerance_deg / 360.0 / target_hz
        )
        permitted_length_m = C_M_PER_S * permitted_delay_s / n_group
        loop_rows.append(
            {
                "loop": label,
                "delay_periods": multiplier,
                "effective_ng": n_group,
                "target_delay_ps": delay_s * 1e12,
                "target_length_mm": length_m * 1e3,
                "phase_tolerance_deg": args.loop_phase_tolerance_deg,
                "permitted_delay_error_ps": permitted_delay_s * 1e12,
                "permitted_length_error_um": permitted_length_m * 1e6,
            }
        )

    cross_section_rows: list[dict[str, Any]] = []
    for optical_row in optical_rows:
        n_group = float(optical_row["ng"])
        delta_n = float(n_rf[target_index] - n_group)
        _, phase_deg, penalty_db = mismatch_metrics(
            np.asarray([solved_target_hz]),
            np.asarray([n_rf[target_index]]),
            n_group,
            electrode_length_m,
        )
        te_fraction = float(optical_row["te_fraction"])
        cross_section_rows.append(
            {
                "width_um": float(optical_row["width_um"]),
                "etch_depth_um": float(optical_row["etch_depth_um"]),
                "mode_label": optical_row["mode_label"],
                "te_fraction": te_fraction,
                "ng": n_group,
                "delta_n_at_target": delta_n,
                "phase_walkoff_deg_per_cm": float(phase_deg[0]),
                "mismatch_penalty_db_per_cm": float(penalty_db[0]),
                "mode_tracking_flag": "通过" if te_fraction >= 0.9 else "需复核",
            }
        )

    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    response_path = output_dir / "速度匹配频率响应.csv"
    loop_path = output_dir / "四程回路长度容差.csv"
    cross_section_path = output_dir / "光学截面敏感性.csv"
    summary_path = output_dir / "速度匹配摘要.json"
    report_path = output_dir / "速度匹配核查报告.md"
    write_csv(response_path, response_rows)
    write_csv(loop_path, loop_rows)
    write_csv(cross_section_path, cross_section_rows)

    rf_validated = (
        args.rf_adaptive_converged and args.rf_port_extra_mode_cleared
    )
    validation_issues: list[str] = []
    if not args.rf_adaptive_converged:
        validation_issues.append("HFSS自适应网格未按设定标准收敛")
    if not args.rf_port_extra_mode_cleared:
        validation_issues.append("端口额外模式尚未排除")

    limitations = [
        f"PDK未公开LN1/LN2实际刻蚀深度与侧壁角；当前{args.etch_depth_um * 1000:g} nm刻蚀是扫描候选，{sidewall_angle_deg:g}°水平夹角是用户指定的名义值",
        "射频介电常数张量、损耗角正切、金属粗糙度和测试背面边界仍需代工或实测确认",
        "HFSS损耗来自相位验证网格，尚未进行最终皮肤深度、airbox和端口尺寸收敛",
        f"回路容差只针对{args.target_frequency_ghz:g} GHz射频包络相位，不是1550 nm载波绝对相位容差",
    ]
    limitations.extend(validation_issues)
    summary = {
        "status": (
            "工程筛选结果，射频自适应已收敛且端口模式已复核；PDK关键参数仍未签核"
            if rf_validated
            else "仅限趋势筛选：射频自适应未收敛且端口存在额外模式提示"
        ),
        "rf_adaptive_converged": args.rf_adaptive_converged,
        "rf_port_extra_mode_cleared": args.rf_port_extra_mode_cleared,
        "rf_port_cross_coupling_max_db_in_selected_band": args.rf_port_cross_coupling_max_db,
        "rf_input": str(args.rf_csv),
        "optical_input": str(args.optical_csv),
        "selected_optical_geometry_um": {
            "width": args.waveguide_width_um,
            "etch_depth": args.etch_depth_um,
            "sidewall_angle_deg_from_horizontal": sidewall_angle_deg,
            "rib_top_width": rib_top_width_um,
            "rib_bottom_width": rib_bottom_width_um,
        },
        "electrode_geometry_um": {
            "signal_width": args.signal_width_um,
            "inner_gap": args.inner_gap_um,
            "t_neck_length": args.t_neck_length_um,
            "t_cap_width": args.t_cap_width_um,
            "t_neck_width": args.t_neck_width_um,
            "t_cap_length": args.t_cap_length_um,
            "t_unit_gap": args.t_unit_gap_um,
            "trunk_gap": args.trunk_gap_um,
            "t_cap_duty_cycle": duty_cycle,
        },
        "target_frequency_requested_ghz": args.target_frequency_ghz,
        "nearest_solved_frequency_ghz": solved_target_ghz,
        "electrode_length_cm": args.electrode_length_cm,
        "n_rf_at_target": float(n_rf[target_index]),
        "rf_loss_db_per_cm_at_target": float(loss_db_per_cm[target_index]),
        "impedance_real_ohm_at_target": float(impedance_real[target_index]),
        "s11_db_at_target": float(s11_db[target_index]),
        "band_ghz": [args.band_start_ghz, args.band_stop_ghz],
        "n_rf_range_in_band": [
            float(np.min(n_rf[band_mask])),
            float(np.max(n_rf[band_mask])),
        ],
        "impedance_real_range_ohm_in_band": [
            float(np.min(impedance_real[band_mask])),
            float(np.max(impedance_real[band_mask])),
        ],
        "mode_summary": mode_summary,
        "loop_phase_tolerance_deg": args.loop_phase_tolerance_deg,
        "loop_length_definition": "从上一程参考面到下一程参考面的总群时延路径，包含调制臂、模式复用器、过渡、弯曲和被动回环",
        "limitations": limitations,
    }
    with summary_path.open("w", encoding="utf-8") as stream:
        json.dump(summary, stream, ensure_ascii=False, indent=2)

    impedance_ok = 45.0 <= impedance_real[target_index] <= 55.0
    velocity_ok = all(
        abs(item["delta_n_at_target"]) <= 0.1
        for item in mode_summary.values()
    )
    numerical_verdict = (
        "同时满足当前初筛的速度匹配（TE0/TE1均为 |Δn|≤0.1）和45–55 Ω阻抗要求"
        if impedance_ok and velocity_ok
        else "尚未同时满足当前初筛的速度匹配（|Δn|≤0.1）和45–55 Ω阻抗要求"
    )
    verdict = (
        numerical_verdict
        if rf_validated
        else f"数值趋势上{numerical_verdict}，但{'；'.join(validation_issues)}，不能判定通过"
    )
    if rf_validated:
        cross_mode_text = (
            f"{args.band_start_ghz:g}–{args.band_stop_ghz:g} GHz内CPW模与奇模的最差交叉耦合为 {args.rf_port_cross_coupling_max_db:.2f} dB。"
            if args.rf_port_cross_coupling_max_db is not None
            else "端口模式已经复核。"
        )
        validation_note = (
            f"> 数值验收：两长度HFSS模型均达到自适应收敛标准；{cross_mode_text}"
        )
    else:
        validation_note = f"> 复核警告：{'；'.join(validation_issues)}。扫频完成不等于数值已经验收。"
    report_lines = [
        "# 速度匹配核查报告",
        "",
        f"当前结论：{args.signal_width_um:g} µm信号主干、{args.inner_gap_um:g} µm内调制间隙、{duty_cycle * 100:.1f}%占空比T形加载的候选，在目标频率附近{verdict}；同时PDK关键截面和材料参数尚未签核，因此当前结果只能保留为候选筛选证据。",
        "",
        validation_note,
        "",
        "## 一、选定截面与射频结果",
        "",
        f"- 光波导：梯形顶宽 {rib_top_width_um:g} µm、底宽 {rib_bottom_width_um:.4f} µm、LN刻蚀 {args.etch_depth_um * 1000:g} nm、侧壁与水平面夹角 {sidewall_angle_deg:g}°；",
        f"- 电极有效长度：{args.electrode_length_cm:g} cm；",
        f"- 设计目标频率：{args.target_frequency_ghz:.2f} GHz；",
        f"- 最近原始扫频点：{solved_target_ghz:.2f} GHz（它只是扫频网格中离{args.target_frequency_ghz:.2f} GHz最近的点）；",
        f"- 射频有效折射率：{n_rf[target_index]:.4f}；",
        f"- 特性阻抗实部：{impedance_real[target_index]:.2f} Ω；",
        f"- 射频筛选损耗：{loss_db_per_cm[target_index]:.2f} dB/cm；",
        f"- 长线模型端口回波：{s11_db[target_index]:.2f} dB。",
        f"- T形尺寸：帽横向宽度 {args.t_cap_width_um:g} µm、颈沿传播方向宽度 {args.t_neck_width_um:g} µm、帽沿传播方向长度 {args.t_cap_length_um:g} µm、周期空隙 {args.t_unit_gap_um:g} µm。",
        "",
        "## 二、TE0与TE1速度匹配",
        "",
        f"| 光模 | 群折射率 | {solved_target_ghz:.2f} GHz折射率差 | 1 cm相位走离 | 1 cm效率惩罚 | 1 dB失配长度 |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for label in ("TE0", "TE1"):
        item = mode_summary[label]
        report_lines.append(
            f"| {label} | {item['ng']:.4f} | {item['delta_n_at_target']:+.4f} | "
            f"{item['phase_walkoff_deg_per_cm_at_target']:+.2f}° | "
            f"{item['mismatch_penalty_db_per_cm_at_target']:.4f} dB | "
            f"{fmt(item['one_db_mismatch_length_cm_at_target'], 2)} cm |"
        )
    report_lines.extend(
        [
            "",
            f"{args.band_start_ghz:g}–{args.band_stop_ghz:g} GHz内，射频有效折射率范围为 {np.min(n_rf[band_mask]):.4f}–{np.max(n_rf[band_mask]):.4f}。TE0和TE1的具体相位走离及1 cm效率惩罚以上表数值为准。",
            "",
            "## 三、四程回路的射频包络容差",
            "",
            "表中的目标长度是从上一程参考面到下一程参考面的总光程，包含调制臂、模式复用器、锥形过渡、弯曲和外部回环，不是只量蛇形直线。",
            "",
            f"以下容差按{args.target_frequency_ghz:g} GHz下每段回路相位误差不超过 ±{args.loop_phase_tolerance_deg:g}°计算。它约束的是调制包络到达时间，不要求1550 nm光载波的绝对相位锁定。",
            "",
            "| 回路 | 目标延迟 | 有效群折射率 | 目标长度 | 允许时延误差 | 允许长度误差 |",
            "|---|---:|---:|---:|---:|---:|",
        ]
    )
    for row in loop_rows:
        report_lines.append(
            f"| {row['loop']} | {row['delay_periods']:g}T | {row['effective_ng']:.4f} | "
            f"{row['target_length_mm']:.3f} mm | ±{row['permitted_delay_error_ps']:.3f} ps | "
            f"±{row['permitted_length_error_um']:.1f} µm |"
        )
    report_lines.extend(
        [
            "",
            "## 四、截面与工艺敏感性",
            "",
            f"本次核查采用顶宽 {rib_top_width_um:g} µm、底宽 {rib_bottom_width_um:.4f} µm、刻蚀 {args.etch_depth_um * 1000:g} nm、侧壁与水平面夹角 {sidewall_angle_deg:g}°的梯形MODE结果，并核对TE0/TE1两个模式。当前输入不覆盖LN刻蚀深度和侧壁角的工艺容差，因此这里只能确认名义截面的速度匹配，不能据此冻结三段回路长度。`光学截面敏感性.csv` 会记录当前输入中的模式，并把TE偏振分数低于0.9的点标为“需复核”。",
            "",
            "## 五、签核前仍需确认",
            "",
            "1. 向代工方取得LN1/LN2实际刻蚀深度、侧壁角和厚度容差；",
            "2. 用确认后的截面重新计算TE0/TE1群折射率，并对选定点做波长步长与网格收敛；",
            "3. 用确认后的LN微波介电张量、金属电导率/粗糙度和背面边界复算HFSS；",
            "4. 对最终候选进行5/25/50 GHz多频自适应、皮肤深度网格和1 cm全长模型；",
            "5. 按模式复用器、弯曲和交叉的逐段群时延重新回标三个回路长度。",
            "",
            "## 六、输出文件",
            "",
            "- `速度匹配频率响应.csv`：20–30 GHz分析所需的逐频数据；",
            "- `四程回路长度容差.csv`：三个回路的目标时延、长度及相位容差；",
            "- `光学截面敏感性.csv`：当前光学输入中各模式的速度失配；",
            "- `速度匹配摘要.json`：机器可读摘要。",
        ]
    )
    report_path.write_text("\n".join(report_lines) + "\n", encoding="utf-8")

    print(f"response={response_path}")
    print(f"loops={loop_path}")
    print(f"cross_sections={cross_section_path}")
    print(f"summary={summary_path}")
    print(f"report={report_path}")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
