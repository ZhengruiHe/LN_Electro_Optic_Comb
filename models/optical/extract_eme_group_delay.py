"""从已求解的EME工程提取模式复用器复传输相位和群时延。

该脚本不重新建立几何，也不覆盖原工程。它调用EME波长扫描，保存目标
通道的复数传输、展开相位以及由相位斜率得到的群时延，供四程回路回标。
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

import numpy as np

from mode_pdk_sweep import DEFAULT_CONFIG, load_config, load_lumapi


C_M_PER_S = 299_792_458.0
ROOT = Path(__file__).resolve().parents[2]
DEFAULT_PROJECT = ROOT / "results" / "optical" / "模式复用器_70度_EME_精细检查.lms"
DEFAULT_OUTPUT = ROOT / "results" / "optical" / "模式复用器_70度_EME_群时延.csv"
DEFAULT_SUMMARY = ROOT / "results" / "optical" / "模式复用器_70度_EME_群时延摘要.json"


CHANNELS = {
    "TE0_正向_s41": "s41",
    "TE0_反向_s14": "s14",
    "TE1转辅助_正向_s52": "s52",
    "辅助转TE1_反向_s25": "s25",
    "独立通道_正向_s63": "s63",
    "独立通道_反向_s36": "s36",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--project", type=Path, default=DEFAULT_PROJECT)
    parser.add_argument("--center-nm", type=float, default=1550.0)
    parser.add_argument("--span-nm", type=float, default=2.0)
    parser.add_argument("--points", type=int, default=9)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    return parser.parse_args()


def scalar_series(value: Any) -> np.ndarray:
    return np.asarray(value, dtype=complex).reshape(-1)


def main() -> None:
    args = parse_args()
    if not args.project.is_file():
        raise FileNotFoundError(f"找不到EME工程：{args.project}")
    if args.points < 5 or args.span_nm <= 0:
        raise ValueError("波长扫描至少需要5点，扫描宽度必须为正数")

    start_nm = args.center_nm - 0.5 * args.span_nm
    stop_nm = args.center_nm + 0.5 * args.span_nm
    cfg = load_config(args.config)
    lumapi = load_lumapi(cfg)
    mode = lumapi.MODE(
        filename=str(args.project.resolve()), hide=bool(cfg["lumerical"]["hide"])
    )
    try:
        print(f"mode_version={mode.version()}", flush=True)
        mode.setemeanalysis("wavelength sweep", 1)
        mode.setemeanalysis("start wavelength", start_nm * 1e-9)
        mode.setemeanalysis("stop wavelength", stop_nm * 1e-9)
        mode.setemeanalysis("number of wavelength points", args.points)
        mode.setemeanalysis("calculate group delays", 1)
        mode.emesweep("wavelength sweep")
        result = mode.getemesweep("S_wavelength_sweep")
    finally:
        mode.close()

    wavelength_m = np.asarray(result["wavelength"], dtype=float).reshape(-1)
    frequency_hz = C_M_PER_S / wavelength_m
    omega = 2.0 * np.pi * frequency_hz
    rows: list[dict[str, float]] = []
    channel_data: dict[str, dict[str, np.ndarray]] = {}
    for label, result_name in CHANNELS.items():
        transmission = scalar_series(result[result_name])
        raw_phase = np.unwrap(np.angle(transmission))
        # MODE 2023 R2导出的正向EME传输相位随角频率增加而增加，和本项目
        # 统一采用的exp(-i*beta*L)传播相位约定相反。先对传输系数取共轭，
        # 再按tau=-d(arg(S))/d(omega)提取，可得到与750um传播时间一致的
        # 正群时延；同时保留原始相位，便于复核求解器符号约定。
        phase = np.unwrap(np.angle(np.conj(transmission)))
        group_delay_s = -np.gradient(phase, omega, edge_order=2)
        channel_data[label] = {
            "transmission": transmission,
            "raw_phase": raw_phase,
            "phase": phase,
            "group_delay_s": group_delay_s,
        }

    for index in range(wavelength_m.size):
        row: dict[str, float] = {
            "wavelength_nm": float(wavelength_m[index] * 1e9),
            "frequency_thz": float(frequency_hz[index] * 1e-12),
        }
        for label, data in channel_data.items():
            value = data["transmission"][index]
            row[f"{label}_实部"] = float(np.real(value))
            row[f"{label}_虚部"] = float(np.imag(value))
            row[f"{label}_功率"] = float(abs(value) ** 2)
            row[f"{label}_求解器原始展开相位_rad"] = float(
                data["raw_phase"][index]
            )
            row[f"{label}_统一约定展开相位_rad"] = float(data["phase"][index])
            row[f"{label}_群时延_ps"] = float(data["group_delay_s"][index] * 1e12)
        rows.append(row)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8-sig") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    center_index = int(np.argmin(np.abs(wavelength_m * 1e9 - args.center_nm)))
    center_channels = {
        label: {
            "result_name": CHANNELS[label],
            "power": float(abs(data["transmission"][center_index]) ** 2),
            "group_delay_ps": float(data["group_delay_s"][center_index] * 1e12),
            "equivalent_group_index_for_750um": float(
                data["group_delay_s"][center_index] * C_M_PER_S / 750e-6
            ),
        }
        for label, data in channel_data.items()
    }
    summary = {
        "status": "EME已求解工程的名义波长扫描；不是工艺容差结果",
        "project": str(args.project.resolve()),
        "center_wavelength_nm": args.center_nm,
        "span_nm": args.span_nm,
        "points": args.points,
        "method": (
            "MODE正向S参数相位符号与项目传播约定相反；对传输系数取共轭后，"
            "按负的相位对角频率斜率计算正群时延"
        ),
        "channels_at_center": center_channels,
        "loop_mux_delay_ps": {
            "回路1_TE0正向加辅助转TE1反向": center_channels["TE0_正向_s41"]["group_delay_ps"]
            + center_channels["辅助转TE1_反向_s25"]["group_delay_ps"],
            "回路2_TE1转辅助正向加TE0反向": center_channels["TE1转辅助_正向_s52"]["group_delay_ps"]
            + center_channels["TE0_反向_s14"]["group_delay_ps"],
            "回路3_TE0正向加辅助转TE1反向": center_channels["TE0_正向_s41"]["group_delay_ps"]
            + center_channels["辅助转TE1_反向_s25"]["group_delay_ps"],
        },
        "limitations": [
            "群时延对应750um模式复用器本体，不包含新增端口锥形段",
            "波长扫描复用已求解EME基底，未对每个波长重新进行独立网格收敛",
        ],
    }
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.summary.write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(f"output={args.output}", flush=True)
    print(f"summary={args.summary}", flush=True)
    print(json.dumps(summary, ensure_ascii=False, indent=2), flush=True)


if __name__ == "__main__":
    main()
