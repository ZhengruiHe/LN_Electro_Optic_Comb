"""逐例运行交叉器名义参数扫描，并自动对前列候选验证第二晶向。

本机一次仅运行一个FDTD，避免16GB内存机器被并发求解耗尽。
已有完成结果按完整参数匹配复用；中断时不覆盖已有工程。
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
FDTD_SCRIPT = Path(__file__).with_name("crossing_fdtd.py")
SOURCE_ROOT = ROOT / "results/optical/交叉器验证"


def parameter_match(data, width, length, source, reference, mesh, mesh_z, slices):
    from crossing_fdtd import MODEL_VERSION
    if data.get("model_version") != MODEL_VERSION:
        return False  # 历史v1实际频点不等于标注波长，不允许当作1550nm缓存。
    a = data["args"]
    if a.get("slab_center_width_um", 7.2) != 7.2:
        return False  # 原LN1扫描只能复用7.2um固定平台，不能混入LN2扫描结果。
    return all(a.get(k) == v for k, v in {
        "width_x_um": width, "width_y_um": width, "half_length_um": length,
        "source": source, "reference_axis": reference, "mesh_nm": mesh,
        "mesh_z_nm": mesh_z, "slices": slices, "port_width_um": .7,
        "wavelength_nm": 1550., "time_fs": 1500., "shutoff": 1e-5,
        "xy_padding_um":4., "z_min_um":-1.3, "z_max_um":3.2, "mask_grid_nm":0.}.items())


def metrics(directory: Path) -> dict:
    r = json.loads((directory / "结果.json").read_text(encoding="utf-8"))
    params = json.loads((directory / "参数与几何.json").read_text(encoding="utf-8"))["args"]
    opposite = {"west":"east", "east":"west", "north":"south", "south":"north"}
    s = r["source"]
    others = [p for p in r["ports"] if p not in (s, opposite[s])]
    return {"directory": str(directory.resolve()), "source": s,
            "center_width_um": params["width_x_um"], "half_length_um": params["half_length_um"],
            "mesh_nm": params["mesh_nm"], "mesh_z_nm": params["mesh_z_nm"], "slices": params["slices"],
            "fdtd_status": r["fdtd_status"],
            "transmission": r["ports"][opposite[s]]["power"],
            "insertion_loss_db": -r["ports"][opposite[s]]["power_db"],
            "reflection_db": r["ports"][s]["power_db"],
            "max_crosstalk_db": max([r["ports"][p]["power_db"] for p in others], default=-300),
            "outgoing_modal_power_sum": sum(p["power"] for p in r["ports"].values())}


def completed_match(width, length, source, reference=None, mesh=100., mesh_z=50., slices=4):
    for parameter_file in SOURCE_ROOT.rglob("参数与几何.json"):
        if not (parameter_file.parent / "结果.json").exists():
            continue
        data = json.loads(parameter_file.read_text(encoding="utf-8"))
        if parameter_match(data, width, length, source, reference, mesh, mesh_z, slices):
            result = metrics(parameter_file.parent)
            if result["fdtd_status"] == 2:
                return result
    return None


def run_case(output, width, length, source, reference=None, mesh=100., mesh_z=50., slices=4):
    cached = completed_match(width, length, source, reference, mesh, mesh_z, slices)
    if cached:
        print("REUSE "+cached["directory"], flush=True)
        return cached
    name = f"W{width:g}_L{length:g}_{source}_d{mesh:g}_z{mesh_z:g}_s{slices}"
    if reference:
        name += "_straight_"+reference
    directory = output / name
    if directory.exists():
        directory = output / (name+"_retry_"+time.strftime("%Y%m%d_%H%M%S"))
    command = [sys.executable, "-X", "utf8", "-u", str(FDTD_SCRIPT),
               "--width-x-um", str(width), "--width-y-um", str(width),
               "--half-length-um", str(length), "--source", source,
               "--mesh-nm", str(mesh), "--mesh-z-nm", str(mesh_z), "--slices", str(slices),
               "--output-dir", str(directory)]
    if reference:
        command += ["--reference-axis", reference]
    print("RUN "+name, flush=True)
    subprocess.run(command, cwd=ROOT, check=True)
    result = metrics(directory)
    if result["fdtd_status"] != 2:
        raise RuntimeError("FDTD未因能量衰减收敛退出，需延长时间窗："+str(directory))
    return result


def write_status(path, data):
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--widths", default="1.5,2,2.5,3")
    p.add_argument("--half-lengths", default="8,12,16")
    p.add_argument("--top-candidates", type=int, default=2)
    p.add_argument("--output-dir", type=Path, default=SOURCE_ROOT/"参数扫描_v1")
    p.add_argument("--refine", action="store_true", help="对双晶向最优候选继续网格/侧壁离散及直波导基准检查")
    args = p.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    # 先检查完整参数族，再启动昂贵求解；未通过手册几何约束的扫描直接停止。
    pdk_dir = args.output_dir / ("PDK前置_"+time.strftime("%Y%m%d_%H%M%S"))
    subprocess.run([sys.executable, "-X", "utf8", str(FDTD_SCRIPT.with_name("crossing_pdk.py")),
                    "--widths", args.widths, "--half-lengths", args.half_lengths,
                    "--output-dir", str(pdk_dir)], cwd=ROOT, check=True)
    status_file = args.output_dir/"扫描状态.json"
    rows, second, refined = [], [], []
    state = {"state":"running", "stage":"parameter_scan", "pid":__import__('os').getpid(),
             "widths_um":[float(x) for x in args.widths.split(',')],
             "half_lengths_um":[float(x) for x in args.half_lengths.split(',')],
             "primary":rows, "second_axis":second, "refinement":refined,
             "nominal_crossing_pass":False,
             "pdk_preflight":str((pdk_dir/"PDK前置检查.json").resolve()),
             "note":"只做名义几何优化；不做制造容差，论文损耗不作为当前结果"}
    try:
        for length in state["half_lengths_um"]:
            for width in state["widths_um"]:
                state["current"] = {"width":width,"half_length":length,"source":"west"}
                write_status(status_file, state)
                rows.append(run_case(args.output_dir, width, length, "west"))
                write_status(status_file, state)
        candidates = sorted(rows, key=lambda v:v["insertion_loss_db"])[:args.top_candidates]
        state["stage"] = "second_crystal_axis"
        for candidate in candidates:
            width, length = candidate["center_width_um"], candidate["half_length_um"]
            state["current"] = {"width":width,"half_length":length,"source":"north"}
            write_status(status_file, state)
            other = run_case(args.output_dir, width, length, "north")
            second.append({"west":candidate, "north":other,
                           "worst_loss_db":max(candidate["insertion_loss_db"],other["insertion_loss_db"])})
            write_status(status_file, state)
        best = min(second, key=lambda v:v["worst_loss_db"])
        state["best_coarse"] = best
        if args.refine:
            state["stage"] = "nominal_numerical_convergence"
            width, length = best["west"]["center_width_um"], best["west"]["half_length_um"]
            for source, axis in (("west","x"),("north","y")):
                # 先分辨横向网格误差，再改变竖直网格与侧壁薄层数。
                for mesh, mesh_z, slices in ((80.,50.,4),(80.,25.,8)):
                    state["current"] = {"width":width,"half_length":length,"source":source,
                                        "mesh_nm":mesh,"mesh_z_nm":mesh_z,"slices":slices}
                    write_status(status_file,state)
                    refined.append(run_case(args.output_dir,width,length,source,mesh=mesh,mesh_z=mesh_z,slices=slices))
                    write_status(status_file,state)
                reference = run_case(args.output_dir,width,length,source,reference=axis,mesh=80.,mesh_z=25.,slices=8)
                state["reference_"+source] = reference
                write_status(status_file,state)
            checks = {}
            for source in ("west","north"):
                items = [r for r in refined if r["source"] == source]
                fine = items[-1]
                reference = state["reference_"+source]
                excess = -10*math.log10(fine["transmission"]/reference["transmission"])
                delta_xy = abs(items[0]["insertion_loss_db"]-best[source]["insertion_loss_db"])
                delta_z = abs(items[1]["insertion_loss_db"]-items[0]["insertion_loss_db"])
                checks[source] = {"fine_insertion_loss_db":fine["insertion_loss_db"],
                                  "straight_normalized_excess_loss_db":excess,
                                  "mesh_xy_change_db":delta_xy,"mesh_z_sidewall_change_db":delta_z,
                                  "max_crosstalk_db":fine["max_crosstalk_db"],"reflection_db":fine["reflection_db"],
                                  "pass": (-.02<=excess<=.3 and delta_xy<=.03 and delta_z<=.03
                                           and fine["max_crosstalk_db"]<=-25 and fine["reflection_db"]<=-25
                                           and fine["outgoing_modal_power_sum"]<=1.01
                                           and abs(reference["insertion_loss_db"])<=.03)}
            state["checks"] = checks
            state["criteria_note"] = "每晶向额外插损≤0.3dB、串扰/反射≤-25dB、两级数值变化≤0.03dB；项目初筛标准，非代工保证"
            state["nominal_crossing_pass"] = all(c["pass"] for c in checks.values())
            design = Path(refined[-1]["directory"])/"参数与几何.json"
            state["selected_design"] = str(design.resolve())
            state["design_sha256"] = hashlib.sha256(design.read_bytes()).hexdigest()
            state["remaining"] = ["1549-1551nm材料色散及复传输群时延", "完整四程光学级联", "官方DRC与功能通过后的容差"]
        state["stage"] = "finished"
        state["state"] = "completed"
        write_status(status_file,state)
        with (args.output_dir/"名义参数扫描.csv").open("w",encoding="utf-8-sig",newline="") as stream:
            writer=csv.DictWriter(stream,fieldnames=list(rows[0]))
            writer.writeheader(); writer.writerows(rows)
    except Exception as error:
        state["state"]="failed"
        state["error"]=str(error)
        write_status(status_file,state)
        raise


if __name__ == "__main__":
    main()
