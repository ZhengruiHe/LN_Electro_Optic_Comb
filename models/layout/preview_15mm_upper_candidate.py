"""生成15mm上方重排局部候选与已检查回程走廊；不是最终GDS。"""
import argparse
import hashlib
import json
from pathlib import Path
from shapely.affinity import translate, scale
from shapely.geometry import box, LineString
from shapely.ops import unary_union
from search_fixed_block_placement import geometry, preview_one_placement


def polygons(shape):
    parts = list(shape.geoms) if hasattr(shape, "geoms") else [shape]
    return [{"exterior": list(p.exterior.coords),
             "holes": [list(ring.coords) for ring in p.interiors]}
            for p in parts if p.geom_type == "Polygon"]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--candidate-dir", type=Path, required=True)
    a = p.parse_args()
    base = Path("results/layout/合并YSJ_版图确认_20260909_v1")
    old = json.loads((base/"YSJ源版图几何.json").read_text(encoding="utf-8"))
    new = json.loads(Path("results/layout/双端GSG150_窗口修正_v18_11/版图预览几何.json").read_text(encoding="utf-8"))
    study = json.loads((a.candidate_dir/"15mm上方模块可行性.json").read_text(encoding="utf-8"))
    dx, dy = -11700., 1600.
    chosen = next(r for r in study["escape_candidates"] if r["dx_um"] == dx and r["dy_um"] == dy)
    core = geometry(new, (20, 0)).intersection(box(850, -43, 18150, 43))
    core = unary_union([
        scale(core.intersection(box(850, -43, 1050, 43)), xfact=.5, yfact=1., origin=(1050, 0)),
        core.intersection(box(1050, -43, 17950, 43)),
        scale(core.intersection(box(17950, -43, 18150, 43)), xfact=.5, yfact=1., origin=(17950, 0)),
    ])
    fixed = core
    corridors = {}
    checks = {}
    for name, raw in study["right_escape_local_points_um"].items():
        points = raw + [[3700., raw[-1][1]]]
        line = LineString(points)
        shape = line.buffer(.35, cap_style=2, join_style=2)
        corridors[name] = shape
        checks[name] = {"simple_centerline": line.is_simple,
                        "YSJ_LN1_overlap_area_um2": geometry(old, (20, 0)).intersection(translate(shape, dx, dy)).area,
                        "corridor_centerline_length_um": line.length,
                        "global_stop_um": [-8000., raw[-1][1]+dy]}
    overlap_pairs = {}
    names = list(corridors)
    for i, name in enumerate(names):
        for other in names[:i]:
            overlap_pairs[name+"__"+other] = corridors[name].intersection(corridors[other]).area
    core = unary_union([fixed, *corridors.values()])
    data = [{"layer": [20, 0], "polygons": polygons(core)},
            {"layer": [42, 0], "polygons": polygons(geometry(new, (42, 0)))}]
    view = preview_one_placement(old, data, dx, dy, a.candidate_dir/"15mm保留候选_左侧延时待布.png",
                                 "15mm保留候选：粉色/蓝色为我们的局部布局；左侧延时与输入输出尚未完成")
    old_budget = json.loads(Path("results/layout/欧拉回路时延预算_v3_无斜直线/欧拉回路_10GHz方向时延预算.json").read_text(encoding="utf-8"))
    delay = []
    for row in old_budget["loops"]:
        reclaimed = row["fixed_two_external_taper_delay_ps"]/2
        delay.append({"loop": row["loop"], "target_total_delay_ps": row["target_total_delay_ps"],
                      "external_tapers_shortened_delay_reduction_ps_estimate": reclaimed,
                      "new_required_passive_route_delay_ps_estimate": row["required_passive_route_delay_ps"]+reclaimed,
                      "basis": "接续拉锥平均ng沿用旧值；不是100um拉锥的新场求解"})
    result = {"status": "partial_geometric_candidate_not_final_layout", "chosen_position": chosen,
              "source_electrode_length_um": 15000, "MUX_body_length_um": 750,
              "external_tapers_candidate_length_um": 100,
              "view_checks": view, "corridor_checks": checks,
              "corridor_pair_overlap_area_um2": overlap_pairs,
              "delay_targets_retained": delay,
              "full_loop_delay_closed": False,
              "remaining": ["左侧三段延时及返回MUX的路线", "输入端与最终输出端", "新100um外接续拉锥的功能依据",
                            "新路线工艺窗口检查与实际GDS回读", "不能以局部候选声称已生成最终合并版"]}
    with (a.candidate_dir/"局部候选与时延保留检查.json").open("x", encoding="utf-8") as f:
        json.dump(result, f, ensure_ascii=False, indent=2)
    with (a.candidate_dir/"局部候选预览几何.json").open("x", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    if any(row["YSJ_LN1_overlap_area_um2"] > 1e-9 or not row["simple_centerline"] for row in checks.values()) or any(v > 1e-9 for v in overlap_pairs.values()):
        raise RuntimeError("局部通道仍存在几何问题")


if __name__ == "__main__":
    main()
