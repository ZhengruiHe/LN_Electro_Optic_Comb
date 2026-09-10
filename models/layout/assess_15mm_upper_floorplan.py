"""15mm不缩短：检查有源段、原MUX、GSG与右端欧拉逃逸段的上方摆位。

仅做几何可行性分析，不改源GDS、不生成最终合并版、不启动场求解。
完整回程、输入输出与时延闭合仍需另外检查。
"""
import argparse
import json
import math
from pathlib import Path
import numpy as np
from shapely.affinity import translate, scale
from shapely.geometry import LineString, box, mapping
from shapely.ops import unary_union
from shapely.prepared import prep
from euler_routes import (append_euler, bend_displacement,
                          radius_for_semicircle_displacement, partial_euler_local)
from search_fixed_block_placement import geometry


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--turn-pattern", choices=("+++", "++-", "+--", "---"), default="+++")
    p.add_argument("--compress-outputs", action="store_true")
    p.add_argument("--compression-length-um", type=float, default=400.)
    p.add_argument("--external-taper-length-um", type=float, default=200.)
    a = p.parse_args()
    a.output_dir.mkdir(parents=True, exist_ok=False)
    base = Path("results/layout/合并YSJ_版图确认_20260909_v1")
    old = json.loads((base/"YSJ源版图几何.json").read_text(encoding="utf-8"))
    new = json.loads(Path("results/layout/双端GSG150_窗口修正_v18_11/版图预览几何.json").read_text(encoding="utf-8"))
    core = geometry(new, (20, 0)).intersection(box(850, -43, 18150, 43))
    if not 100 <= a.external_taper_length_um <= 200:
        raise ValueError("本轮只分析100至200um外端接续长度；不改750um MUX本体")
    if a.external_taper_length_um != 200:
        ratio = a.external_taper_length_um/200.
        core = unary_union([
            scale(core.intersection(box(850, -43, 1050, 43)), xfact=ratio, yfact=1., origin=(1050, 0)),
            core.intersection(box(1050, -43, 17950, 43)),
            scale(core.intersection(box(17950, -43, 18150, 43)), xfact=ratio, yfact=1., origin=(17950, 0)),
        ])
    left_port = 1050-a.external_taper_length_um
    right_port = 17950+a.external_taper_length_um
    metal = geometry(new, (42, 0))
    old_core, old_metal = geometry(old, (20, 0)), geometry(old, (42, 0))
    old_core_p, old_metal_p = prep(old_core), prep(old_metal)
    # 嵌套右端U弯：三个出端口用共同纵向中心，避免相同半径平移后相交。
    d80 = bend_displacement(math.pi, 80)[1]
    first_d = d80 + 3.55
    starts = [29.85, 26.6, -29.85]
    compressed = [29.85, 26.6, 20.1] if a.compress_outputs else starts
    patterns = list(a.turn_pattern)
    names = ["loop1", "loop2", "loop3"]
    if a.compress_outputs and a.turn_pattern != "+++":
        raise ValueError("压缩出口版本当前只检查全部上翻")
    starts = starts + [-26.6]
    compressed = compressed + [23.35 if a.compress_outputs else -26.6]
    patterns = patterns + ["+" if a.compress_outputs else a.turn_pattern[2]]
    names += ["output"]
    centers = {}
    for sign in ("+", "-"):
        members = [y for y, s in zip(compressed, patterns) if s == sign]
        if members:
            centers[sign] = (max(members)+first_d/2 if sign == "+" else min(members)-first_d/2)
    escapes, escape_meta, escape_core = {}, {}, []
    for i, (name, y, target_y, sign) in enumerate(zip(names, starts, compressed, patterns), 1):
        radius = radius_for_semicircle_displacement(2*abs(centers[sign]-target_y), 80)
        points = [[right_port, y]]
        bends = []
        heading = 0.
        if a.compress_outputs:
            delta = target_y-y
            if abs(delta) < 1e-9:
                points.append([right_port+a.compression_length_um, y])
            else:
                angle = 2*math.atan2(delta, a.compression_length_um)
                _, trial = partial_euler_local(angle, 80.)
                rs = 80.*a.compression_length_um/(2*trial["displacement_x_um"])
                for _ in range(3):
                    _, trial = partial_euler_local(angle, rs)
                    rs *= a.compression_length_um/(2*trial["displacement_x_um"])
                if rs < 80:
                    raise ValueError("出口S弯曲率半径小于80um")
                heading = append_euler(points, heading, angle, rs, bends, f"出口{name}渐移1")
                heading = append_euler(points, heading, -angle, rs, bends, f"出口{name}渐移2")
                if abs(points[-1][1]-target_y)>1e-6:
                    raise ValueError("出口S弯高度不闭合")
        points.append([points[-1][0]+5., points[-1][1]])
        append_euler(points, heading, math.pi if sign == "+" else -math.pi, radius, bends, f"右端嵌套欧拉U弯{i}")
        line = LineString(points)
        escapes[name] = points
        escape_core.append(line.buffer(.35, cap_style=2, join_style=2))
        escape_meta[name] = {"radius_um": radius, "return_y_um": points[-1][1],
                                   "length_um": line.length, "bbox_um": list(line.bounds)}
    combined = unary_union([core, *escape_core])
    escape_pair_clear = all(not escape_core[i].intersects(escape_core[j])
                            for i in range(len(escape_core)) for j in range(i))
    checked = 0
    body_only = []
    with_escape = []
    dx_min = -10900+150-left_port
    for dx in np.arange(dx_min, -10379., 20.):
        for dy in np.arange(1330., 1710.01, 10.):
            checked += 1
            moved_m = translate(metal, float(dx), float(dy))
            if old_metal_p.intersects(moved_m):
                continue
            moved_c = translate(core, float(dx), float(dy))
            if old_core_p.intersects(moved_c):
                continue
            row = {"dx_um": float(dx), "dy_um": float(dy),
                   "left_MUX_port_x_um": left_port+dx, "right_MUX_port_x_um": right_port+dx,
                   "M1_clearance_um": old_metal.distance(moved_m),
                   "core_clearance_um": old_core.distance(moved_c)}
            body_only.append(row)
            if combined.bounds[3]+dy+8.6 > 1900:
                continue
            moved_all = translate(combined, float(dx), float(dy))
            if old_core_p.intersects(moved_all):
                continue
            with_escape.append({**row, "core_with_escape_clearance_um": old_core.distance(moved_all),
                                "highest_window_y_um": combined.bounds[3]+dy+8.6})
    result = {
        "status": "floorplan_analysis_not_full_routing", "electrode_length_um": 15000,
        "turn_pattern": a.turn_pattern,
        "compressed_outputs": a.compress_outputs,
        "external_taper_length_um": a.external_taper_length_um,
        "external_taper_geometry_unchanged": a.external_taper_length_um == 200,
        "external_taper_optical_performance_validated": False,
        "checked_positions": checked, "body_and_GSG_clear_count": len(body_only),
        "body_GSG_and_right_escape_clear_count": len(with_escape),
        "right_escape_pairwise_no_overlap": escape_pair_clear,
        "right_escape_geometry": escape_meta,
        "body_candidates": body_only, "escape_candidates": with_escape,
        "right_escape_local_points_um": escapes,
        "limitations": ["只检查原MUX/有源核心/GSG及试算右端U弯", "尚未设计三条完整回程及输入输出",
                        "没有改变原YSJ或源四程GDS", "不将未找到本族候选解释为所有15mm方案均不可行"],
    }
    with (a.output_dir/"15mm上方模块可行性.json").open("x", encoding="utf-8") as f:
        json.dump(result, f, ensure_ascii=False, indent=2)
    brief = {k:v for k,v in result.items() if k not in
             ("body_candidates", "escape_candidates", "right_escape_local_points_um")}
    brief["best_body_candidates"] = sorted(body_only, key=lambda r:r["core_clearance_um"], reverse=True)[:5]
    brief["best_escape_candidates"] = sorted(with_escape, key=lambda r:r["core_with_escape_clearance_um"], reverse=True)[:5]
    print(json.dumps(brief, ensure_ascii=False, indent=2), flush=True)


if __name__ == "__main__":
    main()
