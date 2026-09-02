"""把Maxwell二维准静电场与MODE复数光场做Pockels重叠积分。"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import sys
from typing import Any

import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import numpy as np
from scipy.interpolate import RegularGridInterpolator


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "models" / "optical"))
from mode_pdk_sweep import zelmon_ln_indices  # noqa: E402


C_M_PER_S = 299_792_458.0
EPSILON_0_F_PER_M = 8.854_187_812_8e-12


def read_rf_interpolator(path: Path) -> RegularGridInterpolator:
    data = np.genfromtxt(path, delimiter=",", names=True, encoding="utf-8-sig")
    x_um = np.unique(data["x_um"])
    y_um = np.unique(data["y_um"])
    ex = data["Ex_V_per_m"].reshape(x_um.size, y_um.size)
    if np.isnan(ex).any():
        raise ValueError(f"射频场中存在NaN：{path}")
    return RegularGridInterpolator(
        (x_um, y_um), ex, bounds_error=False, fill_value=0.0
    )


def integrate_2d(values: np.ndarray, x_m: np.ndarray, y_m: np.ndarray) -> float:
    return float(np.trapezoid(np.trapezoid(values, y_m, axis=1), x_m, axis=0))


def ln_mask(x_um: np.ndarray, y_um: np.ndarray) -> np.ndarray:
    angle = math.radians(70.0)
    slab_half_width = 3.6 + (0.9 - y_um) / math.tan(angle)
    rib_half_width = 0.665 + (1.1 - y_um) / math.tan(angle)
    slab = (
        (y_um >= 0.7)
        & (y_um <= 0.9)
        & (np.abs(x_um) <= slab_half_width)
    )
    rib = (
        (y_um >= 0.9)
        & (y_um <= 1.1)
        & (np.abs(x_um) <= rib_half_width)
    )
    return slab | rib


def ln_polygons() -> list[np.ndarray]:
    angle = math.radians(70.0)
    slab_bottom_width = 7.2 + 2.0 * 0.2 / math.tan(angle)
    rib_bottom_width = 1.33 + 2.0 * 0.2 / math.tan(angle)
    return [
        np.asarray(
            [
                [-slab_bottom_width / 2.0, 0.7],
                [slab_bottom_width / 2.0, 0.7],
                [3.6, 0.9],
                [-3.6, 0.9],
            ]
        ),
        np.asarray(
            [
                [-rib_bottom_width / 2.0, 0.9],
                [rib_bottom_width / 2.0, 0.9],
                [0.665, 1.1],
                [-0.665, 1.1],
            ]
        ),
    ]


def mode_overlap(
    optical: Any,
    mode_name: str,
    rf_cap: RegularGridInterpolator,
    rf_trunk: RegularGridInterpolator,
    rail_center_um: float,
    duty_cycle: float,
    r33_m_per_v: float,
    r13_m_per_v: float,
) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    prefix = mode_name + "_"
    x_m = optical[prefix + "x_m"]
    y_m = optical[prefix + "y_m"]
    x_um, y_um = np.meshgrid(x_m * 1e6, y_m * 1e6, indexing="ij")
    points = np.column_stack(((x_um + rail_center_um).ravel(), y_um.ravel()))
    ex_cap = rf_cap(points).reshape(x_um.shape)
    ex_trunk = rf_trunk(points).reshape(x_um.shape)
    ex_weighted = duty_cycle * ex_cap + (1.0 - duty_cycle) * ex_trunk

    ex_opt = optical[prefix + "Ex"]
    ey_opt = optical[prefix + "Ey"]
    ez_opt = optical[prefix + "Ez"]
    hx_opt = optical[prefix + "Hx"]
    hy_opt = optical[prefix + "Hy"]
    poynting_z = 0.5 * np.real(
        ex_opt * np.conj(hy_opt) - ey_opt * np.conj(hx_opt)
    )
    optical_power_w = integrate_2d(poynting_z, x_m, y_m)
    if optical_power_w <= 0:
        raise ValueError(f"{mode_name}的纵向光功率不是正数")

    no, ne = zelmon_ln_indices(np.asarray([float(optical["wavelength_m"]) * 1e6]))
    no = float(no[0])
    ne = float(ne[0])
    mask = ln_mask(x_um, y_um)
    r33_density = np.zeros_like(x_um, dtype=float)
    r13_density = np.zeros_like(x_um, dtype=float)
    r33_density[mask] = (
        EPSILON_0_F_PER_M
        * ne**4
        * r33_m_per_v
        * np.abs(ex_opt[mask]) ** 2
    )
    r13_density[mask] = (
        EPSILON_0_F_PER_M
        * no**4
        * r13_m_per_v
        * (
            np.abs(ey_opt[mask]) ** 2
            + np.abs(ez_opt[mask]) ** 2
        )
    )
    omega_rad_s = 2.0 * math.pi * C_M_PER_S / float(optical["wavelength_m"])

    def beta_terms(field: np.ndarray) -> tuple[float, float, float]:
        r33_integral = integrate_2d(r33_density * field, x_m, y_m)
        r13_integral = integrate_2d(r13_density * field, x_m, y_m)
        scale = omega_rad_s / (4.0 * optical_power_w)
        beta_r33 = scale * r33_integral
        beta_r13 = scale * r13_integral
        return beta_r33, beta_r13, beta_r33 + beta_r13

    cap_r33, cap_r13, cap_beta = beta_terms(ex_cap)
    trunk_r33, trunk_r13, trunk_beta = beta_terms(ex_trunk)
    weighted_r33, weighted_r13, weighted_beta = beta_terms(ex_weighted)

    def vpi_l_v_cm(beta_per_v_m: float) -> float:
        return 100.0 * math.pi / abs(beta_per_v_m)

    result = {
        "mode": mode_name.upper(),
        "solver_mode": int(optical[prefix + "solver_mode"]),
        "neff": float(optical[prefix + "neff"]),
        "te_fraction": float(optical[prefix + "te_fraction"]),
        "central_energy_fraction": float(
            optical[prefix + "central_energy_fraction"]
        ),
        "optical_power_normalization_w": optical_power_w,
        "rail_center_um": rail_center_um,
        "cap_beta_per_v_per_m": cap_beta,
        "trunk_beta_per_v_per_m": trunk_beta,
        "weighted_beta_per_v_per_m": weighted_beta,
        "cap_vpi_l_v_cm": vpi_l_v_cm(cap_beta),
        "trunk_vpi_l_v_cm": vpi_l_v_cm(trunk_beta),
        "weighted_vpi_l_v_cm": vpi_l_v_cm(weighted_beta),
        "weighted_r33_beta_fraction": float(weighted_r33 / weighted_beta),
        "weighted_r13_beta_fraction": float(weighted_r13 / weighted_beta),
        "ln_mean_ex_cap_v_per_m_per_v": float(np.mean(ex_cap[mask])),
        "ln_mean_ex_trunk_v_per_m_per_v": float(np.mean(ex_trunk[mask])),
        "ln_mean_ex_weighted_v_per_m_per_v": float(np.mean(ex_weighted[mask])),
    }
    arrays = {
        "x_um": x_um,
        "y_um": y_um,
        "rf_ex_weighted": ex_weighted,
        "optical_intensity": (
            np.abs(ex_opt) ** 2 + np.abs(ey_opt) ** 2 + np.abs(ez_opt) ** 2
        ),
    }
    return result, arrays


def plot_overlap(
    results: list[dict[str, Any]],
    arrays: list[dict[str, np.ndarray]],
    output: Path,
) -> None:
    plt.rcParams.update(
        {
            "font.sans-serif": ["Microsoft YaHei", "SimHei", "DejaVu Sans"],
            "axes.unicode_minus": False,
            "font.size": 10,
        }
    )
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 3.9), sharex=True, sharey=True)
    image = None
    for ax, result, data in zip(axes, results, arrays):
        x_um = data["x_um"][:, 0]
        y_um = data["y_um"][0, :]
        rf = data["rf_ex_weighted"] / 1e5
        intensity = data["optical_intensity"]
        intensity = intensity / np.max(intensity)
        image = ax.pcolormesh(
            x_um,
            y_um,
            rf.T,
            shading="auto",
            cmap="RdBu_r",
            vmin=0.55,
            vmax=0.9,
        )
        ax.contour(
            x_um,
            y_um,
            intensity.T,
            levels=[0.1, 0.5, 0.9],
            colors=["#222222", "#ffffff", "#ffd400"],
            linewidths=[0.8, 1.0, 1.2],
        )
        for vertices in ln_polygons():
            ax.add_patch(
                Polygon(
                    vertices,
                    closed=True,
                    facecolor="none",
                    edgecolor="#00e5ff",
                    linewidth=1.0,
                )
            )
        ax.set_xlim(-2.5, 2.5)
        ax.set_ylim(0.62, 1.22)
        ax.set_xlabel("相对波导中心的横向位置（µm）")
        ax.set_title(
            f"{result['mode']}：VπL={result['weighted_vpi_l_v_cm']:.2f} V·cm\n"
            f"r33贡献={100.0 * result['weighted_r33_beta_fraction']:.1f}%"
        )
    axes[0].set_ylabel("竖直位置（µm）")
    if image is not None:
        bar = fig.colorbar(image, ax=axes.ravel().tolist(), fraction=0.03, pad=0.025)
        bar.set_label("90%/10%加权横向电场（10⁵ V/m，每1 V线电压）")
    fig.suptitle("PDK有源LN波导：光场等强线与T形电极横向场重叠", fontsize=13)
    fig.subplots_adjust(left=0.08, right=0.90, bottom=0.16, top=0.76, wspace=0.10)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=240, facecolor="white")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--optical-fields",
        type=Path,
        default=ROOT / "results" / "eo" / "active_waveguide_TE0_TE1_complex_fields_1550nm.npz",
    )
    parser.add_argument(
        "--cap-field",
        type=Path,
        default=ROOT / "results" / "eo" / "T形电极_cap_单位电压电场.csv",
    )
    parser.add_argument(
        "--trunk-field",
        type=Path,
        default=ROOT / "results" / "eo" / "T形电极_trunk_单位电压电场.csv",
    )
    parser.add_argument("--duty-cycle", type=float, default=0.9)
    parser.add_argument("--rail-center-um", type=float, default=30.0)
    parser.add_argument("--r33-pm-per-v", type=float, default=30.9)
    parser.add_argument("--r13-pm-per-v", type=float, default=9.6)
    parser.add_argument("--electrode-length-cm", type=float, default=1.0)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "results" / "eo")
    args = parser.parse_args()
    if not 0.0 <= args.duty_cycle <= 1.0:
        raise ValueError("占空比必须位于0到1之间")
    optical = np.load(args.optical_fields)
    rf_cap = read_rf_interpolator(args.cap_field)
    rf_trunk = read_rf_interpolator(args.trunk_field)
    results: list[dict[str, Any]] = []
    arrays: list[dict[str, np.ndarray]] = []
    for mode_name in ("te0", "te1"):
        result, field_arrays = mode_overlap(
            optical,
            mode_name,
            rf_cap,
            rf_trunk,
            args.rail_center_um,
            args.duty_cycle,
            args.r33_pm_per_v * 1e-12,
            args.r13_pm_per_v * 1e-12,
        )
        results.append(result)
        arrays.append(field_arrays)

    vpi_single_pass = [
        item["weighted_vpi_l_v_cm"] / args.electrode_length_cm
        for item in results
    ]
    ideal_four_pass_vpi = 1.0 / (
        2.0 / vpi_single_pass[0] + 2.0 / vpi_single_pass[1]
    )
    summary = {
        "status": "pdk_geometry_overlap_complete_material_coefficients_not_foundry_confirmed",
        "method": "Maxwell2D_unit_line_voltage_field_plus_MODE_complex_field_first_order_Pockels_perturbation",
        "wavelength_nm": float(optical["wavelength_m"]) * 1e9,
        "waveguide_top_width_um": float(optical["width_um"]),
        "ln1_etch_depth_um": float(optical["etch_depth_um"]),
        "ln_sidewall_angle_deg_from_horizontal": float(optical["sidewall_angle_deg"]),
        "active_region_sin_fully_removed_and_sio2_refilled": True,
        "metal_included_in_optical_mode": bool(optical["metal_included"]),
        "t_cap_duty_cycle": args.duty_cycle,
        "electrode_length_cm": args.electrode_length_cm,
        "electro_optic_coefficients_pm_per_v": {
            "r33": args.r33_pm_per_v,
            "r13": args.r13_pm_per_v,
            "provenance": "Ansys_2023_LN_phase_modulator_CHARGE_example_starting_values_not_foundry_data",
        },
        "modes": results,
        "single_pass_vpi_for_selected_length_v": {
            "TE0": vpi_single_pass[0],
            "TE1": vpi_single_pass[1],
        },
        "ideal_four_pass_vpi_before_rf_loss_and_reflection_v": ideal_four_pass_vpi,
        "four_pass_sequence": ["TE0", "TE1", "TE0", "TE1"],
        "limitations": [
            "r33与r13采用Ansys示例文献起始值，不是本PDK代工实测值",
            "二维场按帽区90%与主干区10%加权，未单独积分10 µm颈部三维边缘场",
            "半波电压尚未乘入HFSS沿程损耗、端口失配和四程光学插损",
            "最终流片前还需电极尺寸、顶氧化层厚度和LN电光系数的工艺角点扫描",
        ],
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    summary_path = args.output_dir / "电光重叠_半波电压摘要.json"
    report_path = args.output_dir / "电光重叠与半波电压报告.md"
    figure_path = args.output_dir / "TE0_TE1光场与T形电极横向场重叠.png"
    with summary_path.open("w", encoding="utf-8") as stream:
        json.dump(summary, stream, ensure_ascii=False, indent=2)
    plot_overlap(results, arrays, figure_path)
    report_lines = [
        "# 电光重叠与半波电压结果",
        "",
        "当前PDK几何下的光电重叠积分已经完成。1 cm单程的TE0半波电压为"
        f"{vpi_single_pass[0]:.3f} V，TE1为{vpi_single_pass[1]:.3f} V；"
        f"两种模式各通过两次且相位完全同相时，理想四程半波电压为{ideal_four_pass_vpi:.3f} V。",
        "",
        "## 计算输入",
        "",
        f"- LN波导顶宽：{float(optical['width_um']):.2f} µm；",
        f"- LN1刻蚀：{float(optical['etch_depth_um']) * 1000:.0f} nm；",
        f"- LN1/LN2侧壁与水平面夹角：{float(optical['sidewall_angle_deg']):.0f}°；",
        "- 有源区下方SiN：完全移除并用SiO2回填；",
        f"- T形帽/主干权重：{args.duty_cycle * 100:.0f}%/{(1.0 - args.duty_cycle) * 100:.0f}%；",
        f"- 电光系数：r33={args.r33_pm_per_v:g} pm/V、r13={args.r13_pm_per_v:g} pm/V。",
        "",
        "## 结果",
        "",
        "| 模式 | 帽区VπL | 主干区VπL | 周期加权VπL | r33贡献 |",
        "|---|---:|---:|---:|---:|",
    ]
    for item in results:
        report_lines.append(
            f"| {item['mode']} | {item['cap_vpi_l_v_cm']:.3f} V·cm | "
            f"{item['trunk_vpi_l_v_cm']:.3f} V·cm | "
            f"{item['weighted_vpi_l_v_cm']:.3f} V·cm | "
            f"{100.0 * item['weighted_r33_beta_fraction']:.2f}% |"
        )
    report_lines.extend(
        [
            "",
            "这里的单位电压是信号电极相对两侧地电极的线电压。左右调制区的横向场符号相反、幅度相等；四程回路利用极性补偿使四次相位调制同相叠加。",
            "",
            "## 结论层级",
            "",
            "几何、电场和光场来自当前本机求解；电光系数采用Ansys官方示例脚本中的起始值，不属于PDK保证参数。因此该结果可以替代原先简单的间隙比例估算，但仍是工程预测，不是流片签核或实测半波电压。",
            "",
            "官方方法参考：https://optics.ansys.com/hc/en-us/articles/19435937674387-Thin-Film-Lithium-Niobate-Electro-Optic-Phase-Modulator",
        ]
    )
    report_path.write_text("\n".join(report_lines) + "\n", encoding="utf-8")
    print(
        f"TE0_vpi_l_v_cm={results[0]['weighted_vpi_l_v_cm']:.9f} "
        f"TE1_vpi_l_v_cm={results[1]['weighted_vpi_l_v_cm']:.9f} "
        f"ideal_four_pass_vpi_v={ideal_four_pass_vpi:.9f}"
    )
    print(f"summary={summary_path}")
    print(f"report={report_path}")
    print(f"figure={figure_path}")


if __name__ == "__main__":
    main()
