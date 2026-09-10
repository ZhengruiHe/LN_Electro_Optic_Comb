"""用原生GDS中的实心矩形见证，证明当前两块不能仅靠整体平移无同层冲突。

只读参数：ysj_gds、ours_gds、report_json。适用于本轮固定原坐标、
不缩放不旋转、所有四程工艺/黑盒边界放在21.8×3.8mm内的范围。
不改GDS，不把波导与金属的跨层投影重叠算作冲突。
"""
import hashlib
import json
from pathlib import Path
import pya


def load(path):
    layout = pya.Layout()
    layout.read(str(path))
    if abs(layout.dbu-.001) > 1e-12 or len(layout.top_cells()) != 1:
        raise ValueError("见证复核要求单顶层和1nm格点")
    return layout, layout.top_cell()


def region(layout, top, key):
    return pya.Region(top.begin_shapes_rec(layout.layer(*key))).merged()


def rectangle(bounds):
    return pya.Region(pya.DBox(*bounds).to_itype(.001))


old_path, new_path = Path(ysj_gds), Path(ours_gds)
hashes = {"YSJ": hashlib.sha256(old_path.read_bytes()).hexdigest(),
          "four_pass": hashlib.sha256(new_path.read_bytes()).hexdigest()}
old_layout, old_top = load(old_path)
new_layout, new_top = load(new_path)
old_core = region(old_layout, old_top, (20, 0))
new_core = region(new_layout, new_top, (20, 0))
old_metal = region(old_layout, old_top, (42, 0))
new_metal = region(new_layout, new_top, (42, 0))
old_metal_rect = [50., -1090., 70., 1190.]
new_metal_rect = [2000., 38.5, 17000., 138.5]
old_core_rect = [3807.213, -1143.686, 3808.413, -997.877]
new_core_rect = [2150., 203.431, 18155., 204.131]
contained = {
    "YSJ_vertical_M1": (rectangle(old_metal_rect)-old_metal).is_empty(),
    "four_pass_15mm_upper_ground": (rectangle(new_metal_rect)-new_metal).is_empty(),
    "YSJ_vertical_LN1_branch": (rectangle(old_core_rect)-old_core).is_empty(),
    "four_pass_horizontal_LN1_return": (rectangle(new_core_rect)-new_core).is_empty(),
}
boxes = [new_top.dbbox(li) for li in new_layout.layer_indices()
         if new_layout.get_info(li).layer != 101 and not new_top.dbbox(li).empty()]
bounds = [min(b.left for b in boxes), min(b.bottom for b in boxes),
          max(b.right for b in boxes), max(b.top for b in boxes)]
frame = [-10900., -1900., 10900., 1900.]
shift = [frame[0]-bounds[0], frame[1]-bounds[1],
         frame[2]-bounds[2], frame[3]-bounds[3]]
metal_x_all = (new_metal_rect[0]+shift[2] <= old_metal_rect[0]
               and new_metal_rect[2]+shift[0] >= old_metal_rect[2])
metal_cannot_go_above = new_metal_rect[1]+shift[3] < old_metal_rect[3]
metal_clear_dy_max = min(shift[3], old_metal_rect[1]-new_metal_rect[3])
core_x_all = (new_core_rect[0]+shift[2] <= old_core_rect[0]
              and new_core_rect[2]+shift[0] >= old_core_rect[2])
core_y_envelope = [new_core_rect[1]+shift[1], new_core_rect[3]+metal_clear_dy_max]
core_cross_all_metal_clear_dy = (
    shift[1] <= metal_clear_dy_max
    and core_y_envelope[0] > old_core_rect[1]
    and core_y_envelope[1] < old_core_rect[3]
)
proof_flags = {
    **contained,
    "metal_overlap_x_unavoidable_for_all_dx": metal_x_all,
    "cannot_move_ground_above_YSJ_M1_within_block": metal_cannot_go_above,
    "LN1_overlap_x_unavoidable_for_all_dx": core_x_all,
    "LN1_crossing_for_every_remaining_dy": core_cross_all_metal_clear_dy,
}
dx, dy = -10360., -1330.
translation = pya.Trans(round(dx/.001), round(dy/.001))
core_intersection = old_core & new_core.transformed(translation)
metal_intersection = old_metal & new_metal.transformed(translation)
result = {
    "status": "rigid_translation_impossible" if all(proof_flags.values()) else "proof_not_established",
    "scope": "YSJ固定、四程只整体平移、不缩放旋转；同层20/0和42/0分别判重",
    "frame_um": frame, "source_sha256": hashes,
    "moving_effective_bbox_um": bounds,
    "translation_range_um": {"dx": [shift[0], shift[2]], "dy": [shift[1], shift[3]]},
    "native_containment_and_proof_checks": proof_flags,
    "witness_rectangles_um": {"YSJ_M1": old_metal_rect, "four_pass_M1": new_metal_rect,
                              "YSJ_LN1": old_core_rect, "four_pass_LN1": new_core_rect},
    "necessary_metal_clear_dy_range_um": [shift[1], metal_clear_dy_max],
    "four_pass_return_y_envelope_under_that_constraint_um": core_y_envelope,
    "example_placement": {"dx_um": dx, "dy_um": dy,
        "LN1_overlap_area_um2": core_intersection.area()*1e-6,
        "LN1_overlap_regions": core_intersection.size(),
        "M1_overlap_area_um2": metal_intersection.area()*1e-6,
        "LN1_intersection_boxes_um": [[poly.bbox().left*.001, poly.bbox().bottom*.001,
            poly.bbox().right*.001, poly.bbox().top*.001] for poly in core_intersection.each()]},
    "not_claimed": ["未证明重新规划四程回路后也无法放下", "未将跨层金属/光波导投影重合判为失败",
                    "未验证跨层重合对光损耗或电光相位的影响", "没有修改或输出合并GDS"],
    "sources_unchanged": {"YSJ": hashlib.sha256(old_path.read_bytes()).hexdigest() == hashes["YSJ"],
                          "four_pass": hashlib.sha256(new_path.read_bytes()).hexdigest() == hashes["four_pass"]},
}
with Path(report_json).open("x", encoding="utf-8") as stream:
    json.dump(result, stream, ensure_ascii=False, indent=2)
print(json.dumps(result, ensure_ascii=False, indent=2))
if not all(proof_flags.values()):
    raise RuntimeError("见证条件未全部成立，不得据此声称无解")
