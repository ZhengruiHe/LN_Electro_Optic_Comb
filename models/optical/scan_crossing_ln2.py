"""独立交叉器LN2中心宽度扫描，端口平台固定7.2um，两级200nm/70度不变。

先匹配真实PDK几何报告，再串行求解两个晶向与同长直波导。
时间上限结束的结果仅列为初筛，不因单例status=1而丢弃后续宽度比较。
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time

import numpy as np

from crossing_fdtd import GEOMETRY_VERSION, MODEL_VERSION, assert_wavelength
from crossing_pdk import RULE_VERSION
from scan_crossing_design import metrics

ROOT = Path(__file__).resolve().parents[2]
WORKER = Path(__file__).with_name("crossing_fdtd.py")


def make_plan(widths, time_fs):
    if not widths or len(widths) != len(set(widths)) or 7.2 not in widths:
        raise ValueError("宽度列表需不重复，并包含7.2um基准")
    if not all(math.isfinite(v) and .3 <= v < 12.8 for v in widths):
        raise ValueError("初扫中心宽度需为有限数且小于固定等宽区总长12.8um")
    if not math.isfinite(time_fs) or time_fs <= 0:
        raise ValueError("时间窗需为有限正数")
    cases = []
    for width in widths:
        for source in ("north", "west"):
            cases.append({"name": f"S{width:g}_{source}", "variant": "crossing",
                          "params": {"slab_center_width_um": width, "source": source}})
    for source, axis in (("north", "y"), ("west", "x")):
        cases.append({"name": "straight_"+source, "variant": "straight",
                      "params": {"slab_center_width_um": 7.2, "source": source, "reference_axis": axis}})
    common = {"width_x_um":2.5, "width_y_um":2.5, "half_length_um":16.,
              "port_width_um":.7, "slab_flat_half_length_um":6.4,
              "mesh_nm":80., "mesh_z_nm":25., "slices":8, "mask_grid_nm":0.,
              "xy_padding_um":4., "z_min_um":-1.3, "z_max_um":3.2,
              "wavelength_nm":1550., "time_fs":time_fs, "shutoff":1e-7}
    for case in cases:
        case["params"] = {**common, **case["params"]}
    return cases


def verify_preflight(path, widths):
    report = json.loads(path.read_text(encoding="utf-8"))
    if report.get("rule_version") != RULE_VERSION or not report.get("all_known_geometry_rules_pass"):
        raise ValueError("PDK几何报告未通过或规则版本不同")
    if hashlib.sha256(Path(report["gds"]).read_bytes()).hexdigest() != report["gds_sha256"]:
        raise ValueError("PDK检查后的GDS已变化")
    for width in widths:
        matching = [c for c in report["cases"] if c.get("slab_center_width_um") == width
                    and c["width_x_um"] == 2.5 and c["width_y_um"] == 2.5
                    and c["half_length_um"] == 16. and c["port_width_um"] == .7
                    and c.get("slab_flat_half_length_um") == 6.4
                    and c.get("geometry_version") == GEOMETRY_VERSION]
        if len(matching) != 1 or not matching[0]["known_geometry_rules_pass"] or any(matching[0]["checks"].values()):
            raise ValueError(f"LN2={width:g}um无匹配的PDK通过证据")
    return {"report":str(path.resolve()), "report_sha256":hashlib.sha256(path.read_bytes()).hexdigest(),
            "gds":report["gds"], "gds_sha256":report["gds_sha256"]}


def summarize(rows):
    references = {r["source"]:r for r in rows if r["variant"] == "straight"}
    designs = {}
    for row in rows:
        if row["variant"] != "crossing":
            continue
        ref = references.get(row["source"])
        entry = {"raw_insertion_loss_db":row["insertion_loss_db"],
                 "reflection_db":row["reflection_db"], "max_crosstalk_db":row["max_crosstalk_db"],
                 "energy_shutoff_reached":row["energy_shutoff_reached"],
                 "normalized_excess_loss_db":None}
        if ref:
            entry["normalized_excess_loss_db"] = -10*math.log10(row["transmission"]/ref["transmission"])
            entry["reference_energy_shutoff_reached"] = ref["energy_shutoff_reached"]
            entry["reference_insertion_loss_db"] = ref["insertion_loss_db"]
        designs.setdefault(str(row["slab_center_width_um"]), {})[row["source"]] = entry
    ranking = []
    for width, axes in designs.items():
        if all(s in axes and axes[s]["normalized_excess_loss_db"] is not None for s in ("north","west")):
            ranking.append({"slab_center_width_um":float(width),
                            "worst_normalized_excess_loss_db":max(v["normalized_excess_loss_db"] for v in axes.values()),
                            "all_energy_shutoffs_reached":all(v["energy_shutoff_reached"] and
                                                              v["reference_energy_shutoff_reached"] for v in axes.values())})
    return {"by_width":designs, "preliminary_ranking":sorted(ranking,key=lambda r:r["worst_normalized_excess_loss_db"]),
            "confirmed_nominal_pass":False,
            "note":"仅名义参数初筛；时间、网格、边界、端口模式和GDS格点仍需确认，不自动选择最终版图"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--slab-widths", default="3.6,10.8,7.2")
    parser.add_argument("--time-fs", type=float, default=2500.)
    parser.add_argument("--pdk-report", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    widths = [float(v) for v in args.slab_widths.split(',')]
    plan = make_plan(widths, args.time_fs)
    evidence = verify_preflight(args.pdk_report, widths)
    args.output_dir.mkdir(parents=True, exist_ok=False)
    state_file = args.output_dir/"确认状态.json"
    dependencies = [WORKER, Path(__file__), WORKER.with_name("crossing_pdk.py"),
                    WORKER.with_name("scan_crossing_design.py"), WORKER.with_name("mode_pdk_sweep.py"),
                    WORKER.with_name("mode_pdk_config.json")]
    hashes = {str(p.resolve()):hashlib.sha256(p.read_bytes()).hexdigest() for p in dependencies}
    state = {"state":"running", "pid":os.getpid(), "plan":plan, "completed":[],
             "confirmed_nominal_pass":False, "pdk_preflight":evidence, "dependency_sha256":hashes,
             "started_local":time.strftime("%Y-%m-%d %H:%M:%S"),
             "scope":"独立LN交叉器LN2局部宽度初筛；端口固定；不做工艺容差，不更新整片GDS"}

    def save():
        state["updated_local"] = time.strftime("%Y-%m-%d %H:%M:%S")
        temporary = state_file.with_suffix(".tmp")
        temporary.write_text(json.dumps(state,ensure_ascii=False,indent=2),encoding="utf-8")
        temporary.replace(state_file)

    save()
    try:
        for case in plan:
            for path, expected in hashes.items():
                if hashlib.sha256(Path(path).read_bytes()).hexdigest() != expected:
                    raise RuntimeError("批次运行期间依赖文件变化，停止以避免混用模型："+Path(path).name)
            verify_preflight(args.pdk_report,widths)
            state["current"] = case["name"]
            directory = args.output_dir/case["name"]
            command = [sys.executable,"-X","utf8","-u",str(WORKER)]
            for key,value in case["params"].items():
                command.extend(["--"+key.replace('_','-'),str(value)])
            command.extend(["--output-dir",str(directory.resolve())])
            print("RUN "+case["name"],flush=True)
            child = subprocess.Popen(command,cwd=ROOT)
            state["child_pid"] = child.pid
            save()
            while True:
                try:
                    code = child.wait(timeout=30)
                    break
                except subprocess.TimeoutExpired:
                    save()
            state.pop("child_pid",None)
            if code:
                raise subprocess.CalledProcessError(code,command)
            raw = json.loads((directory/"结果.json").read_text(encoding="utf-8"))
            if (raw.get("model_version") != MODEL_VERSION or not raw.get("spectral_validation_pass")
                    or raw["geometry"].get("geometry_version") != GEOMETRY_VERSION
                    or raw["fdtd_status"] not in (1,2)):
                raise ValueError("未产生当前模型的有效时间窗结果")
            for name, port in raw["ports"].items():
                assert_wavelength({"lambda":[port["actual_wavelength_nm"]*1e-9]},1550.,name)
                if not math.isfinite(port["power"]) or port["power"] < 0:
                    raise ValueError("端口功率无效")
            row = metrics(directory)
            profile_file = directory/f"端口模式_{raw['source']}.npz"
            with np.load(profile_file) as profile:
                assert_wavelength(profile,1550.,"source profile")
                e2 = np.sum(np.abs(profile["E1"].squeeze())**2,axis=-1)
                transverse = "x" if raw["source"] == "north" else "y"
                t, z = profile[transverse].ravel(), profile["z"].ravel()
                peak = np.unravel_index(np.argmax(e2),e2.shape)
                row["source_mode_peak_transverse_um"] = float(t[peak[0]]*1e6)
                row["source_mode_peak_z_um"] = float(z[peak[1]]*1e6)
            row.update({"name":case["name"], "variant":case["variant"],
                        "slab_center_width_um":case["params"]["slab_center_width_um"],
                        "energy_shutoff_reached":raw["fdtd_status"] == 2,
                        "result_level":"时间窗初筛，完整数值收敛待确认"})
            state["completed"].append(row)
            state["checks"] = summarize(state["completed"])
            save()
        state["state"] = "screening_completed"
        state.pop("current",None)
        save()
        print(json.dumps(state["checks"],ensure_ascii=False,indent=2),flush=True)
    except Exception as error:
        state["state"] = "failed"
        state["error"] = str(error)
        save()
        raise


if __name__ == "__main__":
    main()
