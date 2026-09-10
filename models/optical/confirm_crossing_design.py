"""在真实1550nm重新确认当前交叉器候选，依次检查基准、边界、网格和GDS格点。

禁止复用频点错误的v1结果；所有结果写新目录，一次只运行一个FDTD。
不做制造容差或器件几何优化；每组只改变明确列出的数值设置。
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import msvcrt
import os
from pathlib import Path
import subprocess
import sys
import time

from crossing_fdtd import MODEL_VERSION
from scan_crossing_design import metrics
from crossing_pdk import RULE_VERSION

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = Path(__file__).with_name("crossing_fdtd.py")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--width-um", type=float, default=2.5)
    p.add_argument("--half-length-um", type=float, default=16.)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--reuse-pdk-report", type=Path, help="复用几何相同且GDS哈希一致的本地PDK核查，避免依赖打开的KLayout")
    p.add_argument("--previous-status", type=Path, help="只记录中断批次位置，不覆盖、不自动信任其running状态")
    a = p.parse_args()
    a.output_dir.mkdir(parents=True, exist_ok=True)
    # Windows文件锁在进程退出时释放，防止重复求解而无需杀进程或删历史工程。
    with (a.output_dir/"运行锁").open("a+b") as lock:
        lock.seek(0)
        try:
            msvcrt.locking(lock.fileno(), msvcrt.LK_NBLCK, 1)
        except OSError as error:
            raise RuntimeError("同一输出目录已有确认批次在运行") from error
        run_batch(a)


def run_batch(a):
    state_file = a.output_dir/"确认状态.json"
    if state_file.exists():
        raise FileExistsError("已有确认记录；请指定新目录，不覆盖旧批次")
    variants = [("exact", {}), ("straight", {}),
                ("boundary", {"xy_padding_um":8., "z_min_um":-2.3, "z_max_um":4.2}),
                ("mesh60", {"mesh_nm":60.}), ("mask1nm", {"mask_grid_nm":1.})]
    plan = []
    for name, changes in variants:
        for source, axis in (("north", "y"), ("west", "x")):
            params = {"width_x_um":a.width_um, "width_y_um":a.width_um,
                      "half_length_um":a.half_length_um, "port_width_um":.7,
                      "mesh_nm":80., "mesh_z_nm":25., "slices":8,
                      "wavelength_nm":1550., "time_fs":2500., "shutoff":1e-7,
                      "xy_padding_um":4., "z_min_um":-1.3, "z_max_um":3.2,
                      "mask_grid_nm":0., "source":source, **changes}
            if name == "straight":
                params["reference_axis"] = axis
            plan.append({"name":name+"_"+source, "variant":name, "params":params})
    state = {"state":"running", "pid":os.getpid(), "model_version":MODEL_VERSION,
             "plan":plan, "completed":[], "confirmed_nominal_pass":False,
             "physical_geometry_um":{"center_width_x":a.width_um,"center_width_y":a.width_um,
                                     "half_length":a.half_length_um,"port_width":.7},
             "scope":"1550nm名义参数确认；不做工艺容差；不等于整片PDK签核"}
    state["started_local"] = time.strftime("%Y-%m-%d %H:%M:%S")
    state["previous_status"] = str(a.previous_status.resolve()) if a.previous_status else None
    state["script_sha256_at_start"] = hashlib.sha256(SCRIPT.read_bytes()).hexdigest()
    def save():
        state["updated_local"] = time.strftime("%Y-%m-%d %H:%M:%S")
        state_file.write_text(json.dumps(state,ensure_ascii=False,indent=2),encoding="utf-8")
    save()
    try:
        if a.reuse_pdk_report:
            pr = json.loads(a.reuse_pdk_report.read_text(encoding="utf-8"))
            matching = [c for c in pr["cases"] if c["width_x_um"]==a.width_um
                        and c["width_y_um"]==a.width_um and c["half_length_um"]==a.half_length_um
                        and c["port_width_um"]==.7 and c.get("slab_center_width_um",7.2)==7.2]
            if (pr.get("rule_version")!=RULE_VERSION or not pr.get("all_known_geometry_rules_pass")
                    or len(matching)!=1 or not matching[0]["known_geometry_rules_pass"]
                    or any(v!=0 for v in matching[0]["checks"].values())
                    or hashlib.sha256(Path(pr["gds"]).read_bytes()).hexdigest()!=pr["gds_sha256"]):
                raise ValueError("已有PDK报告与当前候选不匹配，或GDS/规则版本已变化")
            state["pdk_preflight"] = str(a.reuse_pdk_report.resolve())
        else:
            pdk_dir = a.output_dir/"PDK前置"
            subprocess.run([sys.executable,"-X","utf8",str(SCRIPT.with_name("crossing_pdk.py")),
                            "--widths",str(a.width_um),"--half-lengths",str(a.half_length_um),
                            "--output-dir",str(pdk_dir)],cwd=ROOT,check=True)
            state["pdk_preflight"] = str((pdk_dir/"PDK前置检查.json").resolve())
        for case in plan:
            state["current"] = case["name"]
            save()
            directory = a.output_dir/case["name"]
            cmd = [sys.executable,"-X","utf8","-u",str(SCRIPT)]
            for key,value in case["params"].items():
                cmd.extend(["--"+key.replace('_','-'),str(value)])
            cmd.extend(["--output-dir",str(directory)])
            print("RUN "+case["name"],flush=True)
            if hashlib.sha256(SCRIPT.read_bytes()).hexdigest()!=state["script_sha256_at_start"]:
                raise RuntimeError("求解器脚本在批次运行中发生改变，停止以避免混用模型")
            child = subprocess.Popen(cmd,cwd=ROOT)
            state["child_pid"] = child.pid
            save()
            while True:
                try:
                    returncode = child.wait(timeout=30)
                    break
                except subprocess.TimeoutExpired:
                    save()  # 每30秒刷新，配合PID检查识别陈旧running记录。
            state.pop("child_pid",None)
            if returncode:
                raise subprocess.CalledProcessError(returncode,cmd)
            raw = json.loads((directory/"结果.json").read_text(encoding="utf-8"))
            if raw.get("model_version") != MODEL_VERSION or not raw.get("spectral_validation_pass"):
                raise ValueError("真实频点或模型版本未通过")
            row = metrics(directory)
            if row["fdtd_status"] != 2:
                raise ValueError("求解未达到能量自动停止门槛")
            row.update({"name":case["name"],"variant":case["variant"],
                        "spectral_validation_pass":True,
                        "actual_wavelengths_nm":{k:v["actual_wavelength_nm"] for k,v in raw["ports"].items()}})
            state["completed"].append(row)
            save()
        by_name = {r["name"]:r for r in state["completed"]}
        checks = {}
        for source in ("north","west"):
            base = by_name["exact_"+source]
            ref = by_name["straight_"+source]
            excess = -10*math.log10(base["transmission"]/ref["transmission"])
            deltas = {v:abs(by_name[v+"_"+source]["insertion_loss_db"]-base["insertion_loss_db"])
                      for v in ("boundary","mesh60","mask1nm")}
            variants_ok = all(r["insertion_loss_db"]<=.3 and r["reflection_db"]<=-25
                              and r["max_crosstalk_db"]<=-25 and r["outgoing_modal_power_sum"]<=1.01
                              for r in state["completed"] if r["source"]==source and r["variant"]!="straight")
            checks[source] = {"normalized_excess_loss_db":excess,"changes_db":deltas,
                              "pass":(-.02<=excess<=.3 and abs(ref["insertion_loss_db"])<=.03
                                       and all(v<=.02 for v in deltas.values()) and variants_ok)}
        state["checks"] = checks
        state["confirmed_nominal_pass"] = all(c["pass"] for c in checks.values())
        state["criteria_note"] = "双晶向额外插损≤0.3dB；数值/边界/格点变化≤0.02dB；项目门槛不是官方标准"
        state["remaining"] = ["PML层数独立检查", "1549-1551nm真实材料色散及群时延",
                              "掩膜外未刻材料与硅衬底影响", "完整级联与整片PDK掩膜修正"]
        state["state"] = "completed"
        state.pop("current",None)
        save()
    except Exception as error:
        state["state"] = "failed"
        state["error"] = str(error)
        save()
        raise


if __name__ == "__main__":
    main()
