"""只改变MUX最小耦合边缘间隙并运行一份独立EME；不覆盖基线结果。"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2), encoding="utf-8")


def read_s(path: Path) -> list[dict]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--minimum-gap-um", type=float, default=1.00)
    parser.add_argument("--aux-start-um", type=float, default=.35)
    parser.add_argument("--time-limit-s", type=int, default=1800)
    args = parser.parse_args()
    if args.minimum_gap_um < .3:
        raise ValueError("最小耦合边缘间隙不能低于0.30um PDK下限")
    if not .3 <= args.aux_start_um <= .4:
        raise ValueError("辅助脊起始宽度应在0.30至0.40um之间")
    args.output_dir.mkdir(parents=True, exist_ok=False)
    config = {
        "wavelength_nm": 1550., "length_um": 750., "approach_length_um": 100.,
        "coupling_length_um": 550., "separation_length_um": 100.,
        "main_start_um": 1.4, "main_stop_um": 1.1, "aux_start_um": .3,
        "aux_stop_um": .4, "end_gap_um": 2.5, "minimum_gap_um": args.minimum_gap_um,
        "aux_start_um": args.aux_start_um,
        "progress_exponent": 2., "progress_ratio": 2.3333333333,
        "geometry_samples": 121, "vertical_slices": 4, "y_span_um": 12.,
        "mesh_cells_y": 240, "mesh_cells_z": 120, "modes": 8,
        "cells": [3, 31, 3], "port1_modes": [1, 3, 4], "port2_modes": [1, 3, 4],
        "change_scope": "仅aux_start_um；不改变耦合间隙、主波导、长度、端口模式或侧壁",
        "pdk_interpretation": "minimum_gap_um是LN1芯层边缘间隙，不是中心线间距",
    }
    baseline_geom = ROOT / "results/optical/模式复用器_70度_EME_精细检查_几何.csv"
    baseline_s = ROOT / "results/optical/模式复用器_70度_EME_精细检查_S矩阵.csv"
    baseline_project = ROOT / "results/optical/模式复用器_70度_EME_精细检查.lms"
    script = ROOT / "models/optical/mode_mux_eme.py"
    state = {"status": "starting", "controller_pid": __import__('os').getpid(),
             "created_at": time.strftime('%Y-%m-%d %H:%M:%S'), "config": config,
             "baseline_geometry_sha256": digest(baseline_geom),
             "baseline_s_sha256": digest(baseline_s), "baseline_project_sha256": digest(baseline_project),
             "solver_script_sha256": digest(script), "full_chain_pass": False}
    save(args.output_dir / "状态.json", state)
    tag = f"aux{args.aux_start_um:g}".replace('.', 'p')
    geometry = args.output_dir / f"MUX_{tag}_gap{args.minimum_gap_um:g}um_几何.csv"
    smatrix = args.output_dir / f"MUX_{tag}_gap{args.minimum_gap_um:g}um_S矩阵.csv"
    project = args.output_dir / f"MUX_{tag}_gap{args.minimum_gap_um:g}um_EME.lms"
    command = [sys.executable, "-X", "utf8", str(script),
               "--minimum-gap-um", str(args.minimum_gap_um), "--end-gap-um", "2.5",
               "--length-um", "750", "--approach-length-um", "100",
               "--coupling-length-um", "550", "--separation-length-um", "100",
               "--main-start-um", "1.4", "--main-stop-um", "1.1",
               "--aux-start-um", str(args.aux_start_um), "--aux-stop-um", ".4",
               "--geometry-samples", "121", "--vertical-slices", "4",
               "--y-span-um", "12", "--mesh-cells-y", "240", "--mesh-cells-z", "120",
               "--modes", "8", "--cells", "3,31,3", "--port1-modes", "1,3,4",
               "--port2-modes", "1,3,4", "--wavelength-nm", "1550",
               "--project", str(project), "--geometry-output", str(geometry), "--s-output", str(smatrix)]
    log = args.output_dir / "EME控制器输出.log"
    try:
        state.update(status="running", command_scope="aux_width_only")
        save(args.output_dir / "状态.json", state)
        with log.open("x", encoding="utf-8") as stream:
            child = subprocess.Popen(command, cwd=ROOT, stdout=stream, stderr=subprocess.STDOUT)
            state.update(child_pid=child.pid, child_started_at=time.strftime('%Y-%m-%d %H:%M:%S'))
            save(args.output_dir / "状态.json", state)
            code = child.wait(timeout=args.time_limit_s)
        state.update(child_pid=None, returncode=code)
        if code:
            raise RuntimeError(f"EME脚本返回{code}，详见{log}")
        if not geometry.is_file() or not smatrix.is_file() or not project.is_file():
            raise RuntimeError("EME完成但缺少几何、S矩阵或LMS工程")
        rows = read_s(geometry)
        gaps = [float(row["gap_um"]) for row in rows]
        if not gaps or abs(min(gaps) - args.minimum_gap_um) > 1e-9:
            raise RuntimeError("输出几何最小间隙与请求值不一致")
        srows = read_s(smatrix)
        powers = sorted((float(row["功率"]), int(row["输出索引"]), int(row["输入索引"]))
                        for row in srows if int(row["输出索引"]) != int(row["输入索引"]))[-12:][::-1]
        result = {
            "status": "eme_completed",
            "wavelength_nm": 1550., "minimum_edge_gap_um": min(gaps),
            "maximum_edge_gap_um": max(gaps), "geometry_rows": len(rows),
            "auxiliary_start_width_um": args.aux_start_um,
            "auxiliary_stop_width_um": .4,
            "selected_port_modes": {"port_1": [1, 3, 4], "port_2": [1, 3, 4]},
            "strongest_off_diagonal_channels": [
                {"power": p, "output_index": o, "input_index": i} for p, o, i in powers
            ],
            "interpretation": "先保存完整6x6端口S矩阵；输出/输入索引需结合端口场形确认，不能仅按矩阵编号称TE1转换",
            "baseline_geometry_sha256": state["baseline_geometry_sha256"],
            "source_baseline_unchanged": digest(baseline_geom) == state["baseline_geometry_sha256"] and digest(baseline_s) == state["baseline_s_sha256"],
            "full_chain_pass": False,
            "tolerance_sweep": False,
        }
        save(args.output_dir / "结果.json", result)
        state.update(status="completed", fdtd_status="not_applicable_EME", result=result)
        save(args.output_dir / "状态.json", state)
    except subprocess.TimeoutExpired:
        state.update(status="failed_timeout", child_pid=None, error="超过EME控制器时间上限；未强制终止其子进程")
        save(args.output_dir / "状态.json", state)
        raise
    except Exception as exc:
        state.update(status="failed", child_pid=None, error=str(exc))
        save(args.output_dir / "状态.json", state)
        raise


if __name__ == "__main__":
    main()
