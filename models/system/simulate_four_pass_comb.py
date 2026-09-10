"""合成 PDK 光模、HFSS 射频线与四程延迟的初步电光梳模型。"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import jv


C_M_PER_S = 299_792_458.0
VALIDATED_RF_DATASETS = {
    "pdk_10GHz_a70_D90_h4_t10_r45_c5_w43_g5_phasecheck_500_1000.csv",
    "pdk_10GHz_a70_regular_w43_g5_phasecheck_500_1000.csv",
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def select_optical_mode(
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
        raise ValueError(f"不能唯一定位光模 {label}: width={width_um}, etch={etch_um}")
    return matches[0]


def write_csv(path: Path, rows: list[dict[str, float | int]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def contiguous_band(
    frequency_ghz: np.ndarray, response: np.ndarray, center_ghz: float
) -> tuple[float, float]:
    center_index = int(np.argmin(np.abs(frequency_ghz - center_ghz)))
    threshold = response[center_index] / math.sqrt(2.0)
    left = center_index
    right = center_index
    while left > 0 and response[left - 1] >= threshold:
        left -= 1
    while right + 1 < len(response) and response[right + 1] >= threshold:
        right += 1
    return frequency_ghz[left], frequency_ghz[right]


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
    parser.add_argument("--sidewall-angle-deg-from-horizontal", type=float, default=70.0)
    parser.add_argument("--target-frequency-ghz", type=float, default=10.0)
    parser.add_argument("--loop1-periods", type=float, default=2.0)
    parser.add_argument("--loop2-periods", type=float, default=2.5)
    parser.add_argument("--loop3-periods", type=float, default=2.0)
    parser.add_argument("--electrode-length-cm", type=float, default=1.0)
    parser.add_argument("--signal-width-um", type=float, default=43.0)
    parser.add_argument("--inner-gap-um", type=float, default=5.0)
    parser.add_argument("--t-neck-length-um", type=float, default=4.0)
    parser.add_argument("--t-cap-width-um", type=float, default=2.0)
    parser.add_argument("--t-neck-width-um", type=float, default=10.0)
    parser.add_argument("--t-cap-length-um", type=float, default=45.0)
    parser.add_argument("--t-unit-gap-um", type=float, default=5.0)
    parser.add_argument("--trunk-gap-um", type=float, default=17.0)
    parser.add_argument("--continuous-vpi-l-v-cm", type=float, default=3.0)
    parser.add_argument(
        "--eo-overlap-json",
        type=Path,
        default=Path("results/eo/电光重叠_半波电压摘要.json"),
        help="存在时优先采用Maxwell二维场与MODE光场的重叠积分结果",
    )
    parser.add_argument("--rf-power-dbm", type=float, default=28.0)
    parser.add_argument(
        "--rf-validation-status",
        choices=("auto", "validated", "screening"),
        default="auto",
        help="auto仅认可脚本内登记的双长度相位验证数据；其余输入按筛选数据处理",
    )
    parser.add_argument("--output-dir", type=Path, default=Path("results/system"))
    args = parser.parse_args()

    rf_rows = read_csv(args.rf_csv)
    rf_frequency_ghz = np.asarray(
        [float(row["frequency_ghz"]) for row in rf_rows]
    )
    rf_n = np.asarray([float(row["n_rf"]) for row in rf_rows])
    rf_loss_db_cm = np.asarray(
        [float(row["loss_db_per_cm"]) for row in rf_rows]
    )
    rf_s11_db = np.asarray([float(row["long_s11_db"]) for row in rf_rows])
    if args.rf_validation_status == "auto":
        rf_validated = args.rf_csv.name in VALIDATED_RF_DATASETS
    else:
        rf_validated = args.rf_validation_status == "validated"
    if rf_validated and "screening_valid" in rf_rows[0]:
        rf_validated = all(int(float(row["screening_valid"])) == 1 for row in rf_rows)

    optical_rows = read_csv(args.optical_csv)
    te0 = select_optical_mode(
        optical_rows,
        args.waveguide_width_um,
        args.etch_depth_um,
        "TE0",
    )
    te1 = select_optical_mode(
        optical_rows,
        args.waveguide_width_um,
        args.etch_depth_um,
        "TE1",
    )
    ng_te0 = float(te0["ng"])
    ng_te1 = float(te1["ng"])

    # 第一、三回路应为整数周期，中间回路应为半整数周期，用于补偿
    # GSG上下调制区的pi极性翻转。周期数通过命令行传入；历史默认值
    # 仍保留为2T、2.5T、2T，实际版图候选使用独立输出目录保存。
    # 这里用TE0、两模平均、TE1作为三个回路的初步有效群折射率，
    # 版图精修时需对各分段群时延积分。
    loop_multipliers = np.asarray(
        [args.loop1_periods, args.loop2_periods, args.loop3_periods]
    )
    loop_ng = np.asarray([ng_te0, 0.5 * (ng_te0 + ng_te1), ng_te1])
    target_hz = args.target_frequency_ghz * 1e9
    loop_delays_s = loop_multipliers / target_hz
    loop_lengths_m = C_M_PER_S * loop_delays_s / loop_ng
    cumulative_delays_s = np.asarray(
        [0.0, loop_delays_s[0], sum(loop_delays_s[:2]), sum(loop_delays_s)]
    )
    polarity = np.asarray([1.0, 1.0, -1.0, -1.0])
    pass_ng = np.asarray([ng_te0, ng_te1, ng_te0, ng_te1])

    frequency_ghz = np.linspace(
        max(1.0, float(rf_frequency_ghz.min())),
        min(30.0, float(rf_frequency_ghz.max())),
        2901,
    )
    frequency_hz = frequency_ghz * 1e9
    n_rf = np.interp(frequency_ghz, rf_frequency_ghz, rf_n)
    loss_db_cm = np.interp(frequency_ghz, rf_frequency_ghz, rf_loss_db_cm)
    s11_db = np.interp(frequency_ghz, rf_frequency_ghz, rf_s11_db)

    if args.t_cap_length_um <= 0 or args.t_unit_gap_um < 0:
        raise ValueError("T帽长度必须为正，T单元间隙不能为负")
    if args.inner_gap_um <= 0 or args.trunk_gap_um <= 0:
        raise ValueError("电极间隙必须为正")

    # T 帽覆盖区使用内调制间隙，其余主干区域用一阶 E~V/g 近似。
    # 这个因子只用于 Vpi 初估，不替代电光场重叠积分。
    duty_cycle = args.t_cap_length_um / (
        args.t_cap_length_um + args.t_unit_gap_um
    )
    inner_gap_um = args.inner_gap_um
    trunk_gap_um = args.trunk_gap_um
    eo_field_factor = duty_cycle + (1.0 - duty_cycle) * (
        inner_gap_um / trunk_gap_um
    )

    electrode_length_m = args.electrode_length_cm * 1e-2
    loss_np_total = np.maximum(loss_db_cm, 0.0) * args.electrode_length_cm / 8.685889638
    rf_average = np.ones_like(loss_np_total)
    nonzero = loss_np_total > 1e-12
    rf_average[nonzero] = (
        1.0 - np.exp(-loss_np_total[nonzero])
    ) / loss_np_total[nonzero]
    reflection = 10.0 ** (s11_db / 20.0)
    mismatch_voltage = np.sqrt(np.maximum(0.0, 1.0 - reflection**2))

    coherent_sum = np.zeros_like(frequency_hz, dtype=complex)
    eo_overlap_used = args.eo_overlap_json.is_file()
    eo_overlap = None
    if eo_overlap_used:
        eo_overlap = json.loads(args.eo_overlap_json.read_text(encoding="utf-8"))
        mode_beta = {
            item["mode"]: float(item["weighted_beta_per_v_per_m"])
            for item in eo_overlap["modes"]
        }
        pass_beta_per_v_per_m = np.asarray(
            [mode_beta["TE0"], mode_beta["TE1"], mode_beta["TE0"], mode_beta["TE1"]]
        )
    else:
        pass_beta_per_v_per_m = np.full(
            4, math.pi / (args.continuous_vpi_l_v_cm * 1e-2)
        )
    electrode_type = (
        str(eo_overlap.get("electrode_type") or "segmented_t")
        if eo_overlap_used
        else "segmented_t"
    )
    electrode_label = "普通CPW" if electrode_type == "regular" else "T形电极"
    passive_lengths_m = loop_lengths_m - electrode_length_m
    if np.any(passive_lengths_m <= 0.0):
        raise ValueError("有源电极长度已经超过至少一段四程总回路长度")
    coherent_phase_per_v = np.zeros_like(frequency_hz, dtype=complex)
    for pass_index in range(4):
        mismatch_x = (
            np.pi
            * frequency_hz
            * (n_rf - pass_ng[pass_index])
            * electrode_length_m
            / C_M_PER_S
        )
        single_pass = np.sinc(mismatch_x / np.pi) * np.exp(-1j * mismatch_x)
        pass_phasor = (
            polarity[pass_index]
            * single_pass
            * np.exp(-1j * 2.0 * np.pi * frequency_hz * cumulative_delays_s[pass_index])
        )
        coherent_sum += pass_phasor
        coherent_phase_per_v += (
            pass_beta_per_v_per_m[pass_index]
            * electrode_length_m
            * pass_phasor
        )

    if eo_overlap_used:
        phase_efficiency_rad_per_v = (
            rf_average * mismatch_voltage * np.abs(coherent_phase_per_v)
        )
        vpi_v = math.pi / np.maximum(phase_efficiency_rad_per_v, 1e-12)
        voltage_gain = np.abs(coherent_phase_per_v) / max(
            pass_beta_per_v_per_m[0] * electrode_length_m, 1e-12
        )
    else:
        voltage_gain = (
            eo_field_factor * rf_average * mismatch_voltage * np.abs(coherent_sum)
        )
        vpi_v = (
            args.continuous_vpi_l_v_cm
            / args.electrode_length_cm
            / np.maximum(voltage_gain, 1e-12)
        )
        phase_efficiency_rad_per_v = math.pi / vpi_v
    band_low, band_high = contiguous_band(
        frequency_ghz, phase_efficiency_rad_per_v, args.target_frequency_ghz
    )

    response_rows: list[dict[str, float | int]] = []
    for index in range(len(frequency_ghz)):
        response_rows.append(
            {
                "frequency_ghz": frequency_ghz[index],
                "n_rf": n_rf[index],
                "rf_loss_db_per_cm": loss_db_cm[index],
                "s11_db": s11_db[index],
                "coherent_pass_enhancement": abs(coherent_sum[index]),
                "total_voltage_gain_vs_continuous_single_pass": voltage_gain[index],
                "phase_efficiency_rad_per_v": phase_efficiency_rad_per_v[index],
                "estimated_vpi_v": vpi_v[index],
            }
        )

    center_index = int(np.argmin(np.abs(frequency_ghz - args.target_frequency_ghz)))
    power_w = 10.0 ** ((args.rf_power_dbm - 30.0) / 10.0)
    voltage_peak_v = math.sqrt(2.0 * 50.0 * power_w)
    beta_rad = math.pi * voltage_peak_v / vpi_v[center_index]
    orders = np.arange(-100, 101)
    amplitude = jv(orders, beta_rad)
    relative_db = 20.0 * np.log10(np.maximum(np.abs(amplitude), 1e-30))
    relative_db -= relative_db.max()
    comb_rows = [
        {
            "order": int(order),
            "offset_ghz": float(order * args.target_frequency_ghz),
            "relative_power_db": float(relative_db[index]),
        }
        for index, order in enumerate(orders)
    ]
    line_counts = {
        f"above_{abs(threshold)}dB": int(np.sum(relative_db >= threshold))
        for threshold in (-20.0, -30.0, -40.0, -60.0)
    }

    output_dir = args.output_dir
    geometry_tag = (
        f"regular_w{args.signal_width_um:g}_g{args.inner_gap_um:g}"
        if electrode_type == "regular"
        else (
            f"h{args.t_neck_length_um:g}_r{args.t_cap_length_um:g}_"
            f"c{args.t_unit_gap_um:g}_w{args.signal_width_um:g}_g{args.inner_gap_um:g}"
        )
    ).replace(".", "p")
    frequency_tag = f"{args.target_frequency_ghz:g}GHz".replace(".", "p")
    response_path = output_dir / f"four_pass_{frequency_tag}_{geometry_tag}_response.csv"
    comb_path = output_dir / f"four_pass_{frequency_tag}_{geometry_tag}_comb_{args.rf_power_dbm:g}dBm.csv"
    summary_path = output_dir / f"four_pass_{frequency_tag}_{geometry_tag}_summary.json"
    report_path = output_dir / f"four_pass_{frequency_tag}_{geometry_tag}_报告.md"
    write_csv(response_path, response_rows)
    write_csv(comb_path, comb_rows)
    limitations = [
        "RF loss comes from the phase-verification mesh without final skin-depth convergence",
        "loop group indices are effective approximations and require routed segment-by-segment delay extraction",
        "comb spectrum is an ideal pure phase-modulation Bessel spectrum without multiplexer reflections or optical propagation loss",
    ]
    if eo_overlap_used:
        limitations.append(
            "Vpi uses a Maxwell2D/MODE overlap integral, but r33 and r13 are Ansys example starting values rather than foundry-confirmed coefficients"
        )
        if electrode_type == "segmented_t":
            limitations.append(
                "The electrostatic period average uses 90% cap and 10% trunk cross-sections and omits the three-dimensional neck-edge correction"
            )
    else:
        limitations.append(
            "Vpi uses the 3 Vcm reference and a gap-weighted field approximation because no electro-optic overlap result was supplied"
        )
    if not rf_validated:
        limitations.append(
            "RF adaptive passes did not converge for the present input data; nRF, impedance, and loss are screening trends only"
        )
    if not rf_validated:
        limitations.append(
            "HFSS port mode isolation has not been cleared and must be revalidated"
        )
    electrode_geometry_um = (
        {
            "signal_width": args.signal_width_um,
            "signal_ground_gap": args.inner_gap_um,
        }
        if electrode_type == "regular"
        else {
            "signal_width": args.signal_width_um,
            "inner_gap": args.inner_gap_um,
            "t_neck_length": args.t_neck_length_um,
            "t_cap_width": args.t_cap_width_um,
            "t_neck_width": args.t_neck_width_um,
            "t_cap_length": args.t_cap_length_um,
            "t_unit_gap": args.t_unit_gap_um,
            "trunk_gap": args.trunk_gap_um,
            "t_cap_duty_cycle": duty_cycle,
        }
    )
    summary = {
        "status": (
            "preliminary_engineering_estimate_not_tapeout_ready"
            if rf_validated
            else "screening_only_rf_not_numerically_validated"
        ),
        "rf_adaptive_converged": rf_validated,
        "rf_port_extra_mode_cleared": rf_validated,
        "target_frequency_ghz": args.target_frequency_ghz,
        "electrode_type": electrode_type,
        "active_electrode_length_mm": args.electrode_length_cm * 10.0,
        "optical_geometry_um": {
            "width": args.waveguide_width_um,
            "etch_depth": args.etch_depth_um,
        },
        "electrode_geometry_um": electrode_geometry_um,
        "ng_te0": ng_te0,
        "ng_te1": ng_te1,
        "loop_delay_multipliers": loop_multipliers.tolist(),
        "loop_length_definition": "total group-delay path from one pass reference plane to the next, including the preceding modulation arm, multiplexers, tapers, bends and passive delay routing",
        "loop_lengths_mm": (loop_lengths_m * 1e3).tolist(),
        "passive_lengths_after_subtracting_active_arm_mm": (
            passive_lengths_m * 1e3
        ).tolist(),
        "eo_overlap_used": eo_overlap_used,
        "eo_overlap_source": str(args.eo_overlap_json) if eo_overlap_used else None,
        "eo_overlap_mode_vpi_l_v_cm": (
            {
                item["mode"]: item["weighted_vpi_l_v_cm"]
                for item in eo_overlap["modes"]
            }
            if eo_overlap_used
            else None
        ),
        "eo_field_factor_gap_approximation": (
            None if eo_overlap_used else eo_field_factor
        ),
        "continuous_reference_vpi_l_v_cm": (
            None if eo_overlap_used else args.continuous_vpi_l_v_cm
        ),
        "n_rf_at_target": float(n_rf[center_index]),
        "rf_loss_db_per_cm_at_target": float(loss_db_cm[center_index]),
        "coherent_pass_enhancement_at_target": float(abs(coherent_sum[center_index])),
        "voltage_gain_at_target": float(voltage_gain[center_index]),
        "phase_efficiency_rad_per_v_at_target": float(
            phase_efficiency_rad_per_v[center_index]
        ),
        "estimated_vpi_at_target_v": float(vpi_v[center_index]),
        "estimated_3db_band_ghz": [float(band_low), float(band_high)],
        "rf_power_dbm": args.rf_power_dbm,
        "rf_peak_voltage_v": voltage_peak_v,
        "modulation_index_rad": float(beta_rad),
        "modulation_index_pi": float(beta_rad / math.pi),
        "comb_line_counts_relative_to_strongest": line_counts,
        "limitations": limitations,
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    with summary_path.open("w", encoding="utf-8") as stream:
        json.dump(summary, stream, ensure_ascii=False, indent=2)

    report_lines = [
        f"# 四程电光梳工程估计（{electrode_label}）",
        "",
        f"当前结论：HFSS传输参数、MODE群折射率与光电场重叠积分已接入四程相干模型；在{args.target_frequency_ghz:g} GHz、{args.rf_power_dbm:g} dBm射频驱动下，理想纯相位调制模型给出{line_counts['above_20.0dB']}根高于最强梳齿−20 dB的谱线。含射频沿程损耗与回波后的有效半波电压工程预测为{vpi_v[center_index]:.3f} V。",
        "",
        "## 输入结构",
        "",
        f"- LN波导：顶宽 {args.waveguide_width_um:g} µm、刻蚀 {args.etch_depth_um * 1000:g} nm、侧壁与水平面夹角 {args.sidewall_angle_deg_from_horizontal:g}°；",
        f"- 光学群折射率：TE0为 {ng_te0:.4f}，TE1为 {ng_te1:.4f}；",
        f"- 有源电极长度：{args.electrode_length_cm * 10.0:g} mm；",
        (
            f"- 普通CPW：{args.signal_width_um:g} µm信号线、{args.inner_gap_um:g} µm连续间隙；"
            if electrode_type == "regular"
            else f"- T形电极：{args.signal_width_um:g} µm信号主干、{args.inner_gap_um:g} µm内间隙、{args.t_cap_length_um:g}/{args.t_unit_gap_um:g} µm加载/空隙、占空比 {duty_cycle * 100:.1f}%；"
        ),
        f"- 四程回路：{loop_multipliers[0]:g}T、{loop_multipliers[1]:g}T、{loop_multipliers[2]:g}T，对应名义长度 {loop_lengths_m[0] * 1e3:.3f}、{loop_lengths_m[1] * 1e3:.3f}、{loop_lengths_m[2] * 1e3:.3f} mm。",
        f"- 扣除有源电极后的无源长度：{passive_lengths_m[0] * 1e3:.3f}、{passive_lengths_m[1] * 1e3:.3f}、{passive_lengths_m[2] * 1e3:.3f} mm。",
        "",
        f"## {args.target_frequency_ghz:g} GHz结果",
        "",
        f"- 射频有效折射率：{n_rf[center_index]:.4f}；",
        f"- 射频损耗：{loss_db_cm[center_index]:.2f} dB/cm；",
        f"- 四程相干增强：{abs(coherent_sum[center_index]):.3f}倍；",
        f"- 四程总相位效率：{phase_efficiency_rad_per_v[center_index]:.3f} rad/V；",
        f"- 有效半波电压工程估计：{vpi_v[center_index]:.3f} V；",
        f"- 模型3 dB射频工作区间：{band_low:.2f}–{band_high:.2f} GHz；",
        f"- {args.rf_power_dbm:g} dBm对应50 Ω负载峰值电压：{voltage_peak_v:.3f} V，调制指数：{beta_rad:.3f} rad。",
        "",
        "## 理想梳齿数",
        "",
        f"- 高于最强梳齿−20 dB：{line_counts['above_20.0dB']}根；",
        f"- 高于最强梳齿−30 dB：{line_counts['above_30.0dB']}根；",
        f"- 高于最强梳齿−40 dB：{line_counts['above_40.0dB']}根。",
        "",
        "这里统计的是理想贝塞尔相位调制谱，不包含模式复用器反射、波导与弯曲损耗、四程之间的幅度不均衡以及光电探测带宽。",
        "",
        "## 仍需完成的签核仿真",
        "",
        "1. 用代工确认或实测的r33、r13以及顶氧化层厚度做半波电压工艺角点；",
        "2. 对金属皮肤深度、粗糙度和材料损耗做网格与工艺角点收敛；",
        "3. 将模式复用器、弯曲、交叉和实际布线路径逐段加入四程群时延；",
        "4. 用实际光学传播损耗、耦合损耗和各程幅度不平衡重新计算可观测梳齿数。",
    ]
    report_path.write_text("\n".join(report_lines) + "\n", encoding="utf-8")

    print(f"response={response_path}")
    print(f"comb={comb_path}")
    print(f"summary={summary_path}")
    print(f"report={report_path}")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
