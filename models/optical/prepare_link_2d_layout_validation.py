"""按最新版GDS提取四组近邻芯层，核对路线后供二维联合仿真使用。

20/1不是残余LN材料边界。本脚本保留既有7.2um等效平台假设，
只把20/0实际多边形作为芯层输入，不修改GDS、不启动光学求解器。
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
from shapely.affinity import translate
from shapely.geometry import Polygon, box
from shapely.ops import unary_union

ROOT = Path(__file__).resolve().parents[2]
PAIRS = ("left_lower", "left_upper", "right_lower", "right_upper")


def write_new(path, data):
    with path.open("x", encoding="utf-8") as stream:
        json.dump(data, stream, ensure_ascii=False, indent=2)


def identify_actual_pair(core, expected, roi):
    """选择目标两路但保留原GDS边界；较近的第三路不能静默排除。"""
    a, b = core.intersection(roi), expected.intersection(roi)
    selected = a.intersection(b.buffer(1.0))
    other = a.difference(b.buffer(1.0))
    gap = None if other.is_empty else selected.distance(other)
    difference = selected.hausdorff_distance(b)
    passed = difference <= .010 and (gap is None or gap >= 10.)
    return {
        "actual_core_identified": passed,
        "missing_core_area_um2": b.difference(a.buffer(.002)).area,
        "additional_core_area_um2": a.difference(b.buffer(.002)).area,
        "route_to_actual_max_shape_difference_um": difference,
        "excluded_other_route_min_edge_gap_um": gap,
        "excluded_distant_routes_not_jointly_simulated": not other.is_empty,
        "actual_core_area_um2": a.area,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gds", type=Path, required=True)
    parser.add_argument("--routes", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=False)
    jobs, configurations = [], {}
    for pair in PAIRS:
        path = args.output_dir / (pair + "_路线截取.json")
        proc = subprocess.run([sys.executable, "-X", "utf8", str(ROOT / "models/optical/prepare_link_2d_pair.py"),
                               "--routes", str(args.routes.resolve()), "--pair", pair,
                               "--output", str(path.resolve())], capture_output=True)
        if proc.returncode:
            raise RuntimeError(proc.stderr.decode("utf-8", errors="replace"))
        config = json.loads(path.read_text(encoding="utf-8"))
        configurations[pair] = config
        x0, y0 = config["origin_global_um"]
        xmin, ymin, xmax, ymax = config["bounds_um"]
        jobs.append({"name": pair, "cell": "EO4P_15MM_GSG150_HORIZONTAL_DRAFT",
                     "bounds_um": [xmin+x0-2, ymin+y0-2, xmax+x0+2, ymax+y0+2]})
    payload = args.output_dir / "原生提取请求.json"
    actual_path = args.output_dir / "实际GDS局部芯层.json"
    write_new(payload, {"gds": str(args.gds.resolve()), "jobs": jobs, "output": str(actual_path.resolve())})
    result = subprocess.run(["C:/Users/PC/AppData/Roaming/KLayout/klayout_app.exe", "-b", "-r",
                             str(ROOT / "models/layout/klayout_extract_optical_roi.py"),
                             "-rd", "payload_file=" + str(payload.resolve())], capture_output=True)
    if result.returncode:
        raise RuntimeError(result.stderr.decode("utf-8", errors="replace"))
    actual = json.loads(actual_path.read_text(encoding="utf-8"))
    checks, ready = [], []
    for job in actual["jobs"]:
        pair = job["name"]
        config = configurations[pair]
        x0, y0 = config["origin_global_um"]
        core = translate(unary_union([Polygon(v["exterior"], v["holes"]) for v in job["core_polygons"]]),
                         xoff=-x0, yoff=-y0)
        expected = unary_union([Polygon(v) for v in config["core_polygons_um"]])
        # 只在实际求解域内比较；PML外的实际延伸长度可以不同。
        roi = box(*config["bounds_um"])
        # 右侧原路径截取曾添加一小段斜线，最大偏离8nm；不把该旧形状导入求解器。
        # 10nm仅为识别对应芯层的核对阈值。实际仿真始终使用原生GDS多边形。
        check = identify_actual_pair(core, expected, roi)
        checks.append({"pair": pair, **check, "sim_domain_um": config["bounds_um"]})
        if not check["actual_core_identified"]:
            continue
        core = core.intersection(expected.buffer(1.0))
        pieces = list(core.geoms) if hasattr(core, "geoms") else [core]
        if any(v.geom_type != "Polygon" or v.interiors for v in pieces):
            raise ValueError("实际芯层包含孔洞或非多边形，不能静默填充")
        config["core_polygons_um"] = [list(v.exterior.coords) for v in pieces]
        config.update(source_gds=str(args.gds.resolve()), source_gds_sha256=actual["sha256"],
                      source_routes_sha256=hashlib.sha256(args.routes.read_bytes()).hexdigest(),
                      core_source="实际GDS20/0原生多边形，仅选择对应两路、平移和仿真域外裁剪",
                      geometry_verified_against_gds=True, verification=checks[-1])
        config["notes"].append("20/1和10/2没有被当作残余LN材料；LN2不在本轮验证范围")
        output = args.output_dir / (pair + "_实际GDS二维配置.json")
        write_new(output, config)
        ready.append({"pair": pair, "config": str(output.resolve())})
    report = {"scope": "实际芯层与近邻二维输入一致性，不是光学通过", "gds_sha256": actual["sha256"],
              "source_unchanged": actual["source_unchanged"], "checks": checks,
              "all_geometry_pass": all(v["actual_core_identified"] for v in checks),
              "ready": ready, "full_chain_pass": False}
    write_new(args.output_dir / "几何准备核对.json", report)
    print(json.dumps(report, ensure_ascii=False, indent=2))
    if not report["all_geometry_pass"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
