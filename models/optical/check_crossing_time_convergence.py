"""独立交叉器时间收敛检查：仅延长时间，不改变几何、材料、网格或停止阈值。

已有数据只读，新工程单独保存。达到时间上限仍可生成比较报告，
但不得把结果稳定与达到能量停止门槛混为一谈，更不等于完整交叉器通过。
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys
import time

from crossing_fdtd import MODEL_VERSION, assert_wavelength
from scan_crossing_design import metrics

ROOT = Path(__file__).resolve().parents[2]
WORKER = Path(__file__).with_name("crossing_fdtd.py")


def validate_baseline(raw, params, time_fs):
    if raw.get("model_version") != MODEL_VERSION or not raw.get("spectral_validation_pass"):
        raise ValueError("基准不是已核实真实频点的当前模型")
    if raw.get("fdtd_status") not in (1, 2):
        raise ValueError("基准未正常生成时间窗结果")
    if not math.isfinite(time_fs) or time_fs <= params["time_fs"]:
        raise ValueError("新时间必须有限且大于基准时间")
    if raw["source"] != params["source"]:
        raise ValueError("基准输入方向不一致")
    for name, port in raw["ports"].items():
        assert_wavelength({"lambda": [port["actual_wavelength_nm"] * 1e-9]},
                          params["wavelength_nm"], name)


def compare_results(old, new):
    source = old["source"]
    if source != new["source"] or old["ports"].keys() != new["ports"].keys():
        raise ValueError("新旧端口不匹配")
    target = {"north": "south", "south": "north", "east": "west", "west": "east"}[source]
    a, b = old["ports"][target], new["ports"][target]
    sa, sb = complex(**a["S"]), complex(**b["S"])
    phase = math.degrees(math.atan2((sb / sa).imag, (sb / sa).real))
    delta = abs(a["power_db"] - b["power_db"])
    energy_ok = new["fdtd_status"] == 2
    return {
        "insertion_loss_change_db": delta,
        "through_phase_change_deg": phase,
        "port_power_db_changes": {name: new["ports"][name]["power_db"] - value["power_db"]
                                  for name, value in old["ports"].items()},
        "time_window_loss_stable": delta <= .02,
        "energy_shutoff_reached": energy_ok,
        "time_check_pass": delta <= .02 and energy_ok,
        "criteria_note": "插损变化≤0.02 dB且达到原能量停止门槛；相位变化另行报告，不代表网格/边界通过",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--time-fs", type=float, default=5000.)
    args = parser.parse_args()
    baseline = args.baseline_dir.resolve()
    old_file = baseline / "结果.json"
    param_file = baseline / "参数与几何.json"
    old = json.loads(old_file.read_text(encoding="utf-8"))
    params = json.loads(param_file.read_text(encoding="utf-8"))["args"]
    validate_baseline(old, params, args.time_fs)
    # 只替换计算时间和输出位置；不从其他扫描目录猜测或拼接参数。
    new_params = {k: v for k, v in params.items() if k not in ("output_dir", "build_only")}
    new_params["time_fs"] = args.time_fs
    args.output_dir.mkdir(parents=True, exist_ok=False)
    project_dir = args.output_dir / f"{params['source']}_{args.time_fs:g}fs"
    state_file = args.output_dir / "确认状态.json"
    state = {"state": "running", "pid": os.getpid(), "completed": [],
             "plan": [{"name": project_dir.name, "params": new_params}],
             "current": project_dir.name, "confirmed_nominal_pass": False,
             "baseline": str(baseline), "scope": "独立交叉器单方向时间收敛；不做工艺容差",
             "started_local": time.strftime("%Y-%m-%d %H:%M:%S"),
             "worker_sha256": hashlib.sha256(WORKER.read_bytes()).hexdigest(),
             "baseline_result_sha256": hashlib.sha256(old_file.read_bytes()).hexdigest(),
             "baseline_parameters_sha256": hashlib.sha256(param_file.read_bytes()).hexdigest()}

    def save():
        state["updated_local"] = time.strftime("%Y-%m-%d %H:%M:%S")
        state_file.write_text(json.dumps(state, ensure_ascii=False, indent=2), encoding="utf-8")

    save()
    command = [sys.executable, "-X", "utf8", "-u", str(WORKER)]
    for key, value in new_params.items():
        if value is not None:
            command.extend(["--" + key.replace("_", "-"), str(value)])
    command.extend(["--output-dir", str(project_dir.resolve())])
    try:
        print("开始时间收敛检查：" + project_dir.name, flush=True)
        child = subprocess.Popen(command, cwd=ROOT)
        state["child_pid"] = child.pid
        save()
        while True:
            try:
                code = child.wait(timeout=30)
                break
            except subprocess.TimeoutExpired:
                save()
        state.pop("child_pid", None)
        if code:
            raise subprocess.CalledProcessError(code, command)
        if hashlib.sha256(WORKER.read_bytes()).hexdigest() != state["worker_sha256"]:
            raise ValueError("执行期间建模脚本改变，不自动认可结果")
        new = json.loads((project_dir / "结果.json").read_text(encoding="utf-8"))
        validate_baseline(new, new_params, args.time_fs + 1)
        state["completed"].append(metrics(project_dir))
        state["checks"] = compare_results(old, new)
        # 日志只保留明确的数值进度行；不复制许可证、主机或凭据上下文。
        progress = []
        for log in project_dir.glob("*_p0.log"):
            for line in log.read_text(encoding="utf-8", errors="replace").splitlines():
                if re.match(r"^\s*\d+(?:\.\d+)?% complete\.", line) and "Auto Shutoff:" in line:
                    progress.append(line)
        state["last_numeric_progress"] = progress[-5:]
        state["state"] = "completed" if state["checks"]["time_check_pass"] else "needs_numerical_review"
        state["remaining"] = ["另一晶向", "直波导归一化", "网格与边界收敛",
                              "时间检查未通过时定位残余场；不自动放宽门槛或改成黑盒通过"]
        state.pop("current", None)
        save()
        print(json.dumps(state["checks"], ensure_ascii=False, indent=2), flush=True)
    except Exception as error:
        state["state"] = "failed"
        state["error"] = str(error)
        save()
        raise


if __name__ == "__main__":
    main()
