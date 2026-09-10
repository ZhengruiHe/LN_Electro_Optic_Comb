"""将自定义交叉器的实际轮廓接入四程GDS，始终写入新文件。

使用独立KLayout内存对象读取历史基线，不清空/重载用户活动页。
损耗验证尚未完成时只允许显式输出带“待验证”标签的工作稿。
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from klink import KLinkClient
from shapely.affinity import translate
from shapely.geometry import LineString, Polygon, box
from shapely.ops import unary_union

from check_crossing_candidate import build
from route_geometry import check_route_network

ROOT = Path(__file__).resolve().parents[2]
OLD_TOP = "EO4P_10G_15MM_3T_3P5T_3T_NOMINAL"
NEW_TOP = "EO4P_10G_15MM_CUSTOM_CROSSING"


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--design", type=Path, required=True, help="FDTD参数与几何JSON")
    p.add_argument("--validation", type=Path, help="交叉器独立验证报告")
    p.add_argument("--draft", action="store_true", help="明确生成待验证工作稿，不宣布光学通过")
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--fix-join-notches", action="store_true",
                   help="仅在四处MUX外端汇合处加入明确矩形包覆补片，不修改LN核心")
    a = p.parse_args()
    if not a.draft:
        raise ValueError("PDK复核发现历史基线缺少LN2显式刻蚀映射；修正整片两次刻蚀及严格DRC之前，仅允许--draft输出。光学低损耗不能解除此门槛。")
    # 工作稿也不接受与几何不一致的验证记录；最终导出门槛在上方明确阻断。
    if a.validation:
        validation = json.loads(a.validation.read_text(encoding="utf-8"))
        if not validation.get("nominal_crossing_pass"):
            raise ValueError("所附交叉器验证报告未通过；无验证工作稿请省略--validation")
        if validation.get("design_sha256") != hashlib.sha256(a.design.read_bytes()).hexdigest():
            raise ValueError("验证报告与拟使用几何不一致")
    design = json.loads(a.design.read_text(encoding="utf-8"))
    if design["args"]["reference_axis"] is not None:
        raise ValueError("直波导基准不是交叉器")
    if design["args"]["port_width_um"] != .7:
        raise ValueError("当前四程接口必须是0.7um")
    a.output_dir.mkdir(parents=True, exist_ok=True)
    tag = "待验证工作稿" if a.draft else "交叉器名义验证版_整链未签核"
    output = (a.output_dir / f"四程10GHz_15mm_自定义交叉_{tag}.gds").resolve()
    clean = (a.output_dir / f"四程10GHz_15mm_核心工艺层_{tag}.gds").resolve()
    if output.exists() or clean.exists():
        raise FileExistsError("不覆盖历史GDS，请指定新输出目录")
    baseline = ROOT / "results/layout/四程电光梳_10GHz_15mm_3T_3.5T_3T_名义_工作版.gds"
    budget = json.loads((ROOT / "results/layout/四程电光梳_10GHz_15mm_3T_3.5T_3T_名义_检查.json").read_text(encoding="utf-8"))
    lengths = [v["actual_centerline_length_um"] for v in budget["loop_centerline_checks"]]
    routes = build(lengths, True)
    crossing = {"routes": ["input", "loop2"], "center_um": [300., 250.], "angle_deg": 90}
    topology = check_route_network(routes, [crossing])
    if not topology["topology_screen_pass"]:
        raise ValueError("非预期光路交叉尚未关闭")
    # 实际中心扩展轮廓来自FDTD；在仿真域外延长相同0.7um端口到45um单元边缘。
    unit = box(-22.5, -22.5, 22.5, 22.5)
    core = unary_union([Polygon(design["core_top_polygon_um"]),
                        box(-22.5, -.35, 22.5, .35), box(-.35, -22.5, .35, 22.5)]).intersection(unit)
    slab = unary_union([box(-22.5, -3.6, 22.5, 3.6), box(-3.6, -22.5, 3.6, 22.5)])
    trench = unary_union([box(-22.5, -8.6, 22.5, 8.6), box(-8.6, -22.5, 8.6, 22.5)])
    aperture = translate(unit, 300, 250)
    shape_groups = {"20/0": [], "20/1": [], "10/2": []}
    for layer, width, custom in (("20/0", .7, core), ("20/1", 7.2, slab), ("10/2", 17.2, trench)):
        geometries = []
        for points in routes.values():
            # 先缓冲为真实宽度再切除窗口，保证边界覆盖恰与交叉器端面一致。
            shape = LineString(points).buffer(width/2, cap_style=2, join_style=1, quad_segs=8)
            geometries.append(shape.difference(aperture))
        geometries.append(translate(custom, 300, 250))
        merged = unary_union(geometries)
        polys = [merged] if merged.geom_type == "Polygon" else list(merged.geoms)
        # 长蛇形包络可能形成环孔；KLayout区域支持显式孔。
        shape_groups[layer] = [{"exterior": list(poly.exterior.coords)[:-1],
                                "holes": [list(r.coords)[:-1] for r in poly.interiors]} for poly in polys]
    payload = {"baseline": str(baseline.resolve()), "output": str(output), "clean": str(clean),
               "old_top": OLD_TOP, "new_top": NEW_TOP, "shapes": shape_groups,
               "fix_joins": a.fix_join_notches}
    # pya仅在独立内存layout中操作；磁盘原文件不修改，用户活动页不受影响。
    code = "payload = " + repr(payload) + "\n" + '''
import json
candidate_layout = pya.Layout()
candidate_layout.read(payload["baseline"])
candidate_top = candidate_layout.cell(payload["old_top"])
if candidate_top is None:
    raise RuntimeError("基线顶层不匹配")
candidate_top.name = payload["new_top"]
removed_counts = {}
for layer_key, polygons in payload["shapes"].items():
    ln, dt = map(int, layer_key.split("/"))
    li = candidate_layout.layer(ln, dt)
    removed_counts[layer_key] = candidate_top.shapes(li).size()
    candidate_top.shapes(li).clear()
    for polygon in polygons:
        dp = pya.DPolygon([pya.DPoint(*xy) for xy in polygon["exterior"]])
        for hole in polygon["holes"]:
            dp.insert_hole([pya.DPoint(*xy) for xy in hole])
        candidate_top.shapes(li).insert(dp.to_itype(candidate_layout.dbu))
# 只平移明确标记的旧输入端面黑盒，不动输出耦合器/有源电极/MUX。
ec_moved = 0
for inst in candidate_top.each_inst():
    if "edge_coupler_9um_y_1550_ln" in inst.cell.name:
        pos = inst.dtrans.disp
        if abs(pos.x-300) < 0.002 and abs(pos.y-30.15) < 0.002:
            trans = inst.dcplx_trans
            trans.disp = pya.DVector(0,250)
            inst.dcplx_trans = trans
            ec_moved += 1
if ec_moved != 1:
    raise RuntimeError("未唯一识别基线输入耦合器，拒绝继续")
# 清理旧位置的调试文字，在新图中只保留明确状态。
li_text = candidate_layout.layer(101,0)
candidate_top.shapes(li_text).clear()
candidate_top.shapes(li_text).insert(pya.DText("CUSTOM_CROSSING_1550_NOMINAL_NOT_TAPEOUT", 1800,970).to_itype(candidate_layout.dbu))
candidate_top.shapes(li_text).insert(pya.DText("X_SIM=LN_Y; Y_SIM=LN_Z", 200,310).to_itype(candidate_layout.dbu))
changes = {}
if payload["fix_joins"]:
    # 由初稿DRC交点坐标确定的显式补片。避免布尔闭运算产生纳米级孤立碎片。
    patch_boxes = {"20/1":[(798,33.5,804,35),(18178,30,18183,32)],
                   "10/2":[(760,38.5,766,40),(18200,35,18206,39)]}
    for layer_key, upper_boxes in patch_boxes.items():
        ln,dt = map(int,layer_key.split("/"))
        li = candidate_layout.layer(ln,dt)
        region = pya.Region(candidate_top.begin_shapes_rec(li)).merged()
        patches = pya.Region()
        for x0,y0,x1,y1 in upper_boxes:
            for coords in ((x0,y0,x1,y1),(x0,-y1,x1,-y0)):
                rect = pya.DBox(*coords).to_itype(candidate_layout.dbu)
                patches.insert(rect)
                candidate_top.shapes(li).insert(rect)
        added = patches-region
        changes[layer_key] = {"added_area_um2": added.area()*candidate_layout.dbu**2,
                              "removed_area_um2": 0,
                              "scope": "仅MUX外端两对分支汇合窗口"}
candidate_layout.write(payload["output"])
# 工艺层导出以独立副本展平，不在工作版中丢弃器件层级和耦合器黑盒。
process_layout = pya.Layout()
process_layout.read(payload["output"])
process_top = process_layout.cell(payload["new_top"])
process_top.flatten(True)
allowed = {(20,0),(20,1),(10,2),(42,0)}
for li in list(process_layout.layer_indices()):
    inf = process_layout.get_info(li)
    if (inf.layer,inf.datatype) not in allowed:
        process_layout.delete_layer(li)
process_layout.write(payload["clean"])
{"removed_legacy_route_shapes":removed_counts, "input_coupler_moved":ec_moved,
 "process_bbox_um":[process_top.dbbox().left,process_top.dbbox().bottom,process_top.dbbox().right,process_top.dbbox().top],
 "join_mask_changes":changes}
'''
    c = KLinkClient()
    c.connect()
    try:
        result = c.call("exec.python", {"code": code})
        if result.get("exception"):
            raise RuntimeError(str(result["exception"]))
        details = result["return_value"]
        info = c.call("layout.file_info", {"path": str(output), "detail": "counts"})
        clean_info = c.call("layout.file_info", {"path": str(clean), "detail": "counts"})
    finally:
        c.close()
    report = {"status":tag, "working_gds":str(output), "process_gds":str(clean),
              "top_cell":NEW_TOP, "design":str(a.design.resolve()),
              "design_sha256":hashlib.sha256(a.design.read_bytes()).hexdigest(),
              "validation":str(a.validation.resolve()) if a.validation else None,
              "topology":topology, "kLayout":details, "working_info":info, "process_info":clean_info,
              "limitations":["中心线等长尚不等于晶向积分后的真实时延相等，10GHz最终时延需回填",
                             "模式复用器新接续、弯曲、焊盘及完整级联仍未签核",
                             "工艺层文件不包含PDK端面耦合器真实掩膜",
                             "若启用包络接缝清理，局部MUX外端包覆补片已同步进入工作稿/工艺层；完整级联需使用该几何"]}
    (a.output_dir/"版图检查.json").write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({"gds":str(output),"process_gds":str(clean),"details":details},ensure_ascii=False,indent=2))


if __name__ == "__main__":
    main()
