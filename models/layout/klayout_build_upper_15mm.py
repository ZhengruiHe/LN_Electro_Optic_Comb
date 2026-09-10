"""KLayout原生构建15mm上方功能版图候选；只写新文件，保持源单元不变。"""
import hashlib
import json
from pathlib import Path
import pya


with open(payload_file, encoding="utf-8") as stream:
    p = json.load(stream)
source = pya.Layout(); source.read(p["source_gds"])
layout = pya.Layout(); layout.dbu = source.dbu
if abs(layout.dbu-.001)>1e-12:
    raise ValueError("需要1nm格点")
top = layout.create_cell(p["top_name"])
shift = pya.DCplxTrans(1, 0, False, *p["shift_um"])
source_top = source.top_cell()
copied = {}
copy_normalizations = {}


def digest_cell(lib, cell):
    values = [(str(lib.get_info(li)), sorted(str(s) for s in cell.shapes(li).each()))
              for li in lib.layer_indices() if not cell.shapes(li).is_empty()]
    insts = sorted((i.cell.name, str(i.dcplx_trans)) for i in cell.each_inst())
    return hashlib.sha256(json.dumps([sorted(values), insts]).encode()).hexdigest()


def copy_cell(name):
    if name not in copied:
        src = source.cell(name)
        if src is None:
            raise ValueError("源单元缺失："+name)
        dest = layout.create_cell(name)
        dest.copy_tree(src)
        if digest_cell(source, src) != digest_cell(layout, dest):
            equality={}
            for sli in source.layer_indices():
                key=source.get_info(sli)
                before=pya.Region(src.shapes(sli)).merged()
                after=pya.Region(dest.shapes(layout.layer(key))).merged()
                texts_before=sorted(str(s) for s in src.shapes(sli).each() if s.is_text())
                texts_after=sorted(str(s) for s in dest.shapes(layout.layer(key)).each() if s.is_text())
                equality[str(key)]=(before ^ after).is_empty() and texts_before==texts_after
            equality['instances']=sorted((i.cell.name,str(i.dcplx_trans)) for i in src.each_inst())==sorted((i.cell.name,str(i.dcplx_trans)) for i in dest.each_inst())
            if not all(equality.values()):
                raise ValueError("单元复制改变了图形/文字/层级："+name+str(equality))
            copy_normalizations[name]={'reason':'原生复制去除冗余共线顶点；几何异或为空，文字及实例相同','checks':equality}
        copied[name] = dest
    return copied[name]


def reg(lib, cell, key):
    return pya.Region(cell.begin_shapes_rec(lib.layer(*key))).merged()


def rect(bounds):
    return pya.Region(pya.DBox(*bounds).to_itype(layout.dbu))


def polygon(record):
    poly = pya.DPolygon([pya.DPoint(*xy) for xy in record["exterior"]])
    for hole in record["holes"]:
        poly.insert_hole([pya.DPoint(*xy) for xy in hole])
    return poly.to_itype(layout.dbu)


mux_source_name = "MUX_TE01_A70_TAIL030_EXTERNAL_WINDOWS_WIDENED_UNVERIFIED"
mux = copy_cell(mux_source_name)
active_um=float(p.get('active_taper_um',200.))
external_um=float(p.get('external_taper_um',100.))
revised=active_um!=200. or external_um!=100.
normalize_mux_windows=bool(p.get('normalize_mux_windows_to_pdk',False))
normalize_mux_windows=bool(p.get('normalize_mux_windows_to_pdk',False))
mux.name = f"MUX_TE01_A70_A{active_um:g}_E{external_um:g}_BODY750_DRAFT" if revised else "MUX_TE01_A70_EXT100_BODY750_DRAFT"
mux_body_checks = {}
aux_tail_checks = {}

def mapped_x_region(region,fn):
    result=pya.Region()
    for poly in region.each():
        def point(pt):return pya.Point(round(fn(pt.x)),pt.y)
        new_poly=pya.Polygon([point(pt) for pt in poly.each_point_hull()])
        for h in range(poly.holes()):new_poly.insert_hole([point(pt) for pt in poly.each_point_hole(h)])
        result.insert(new_poly)
    return result

for li in list(layout.layer_indices()):
    info = layout.get_info(li)
    if mux.shapes(li).is_empty():
        continue
    text_items = [s.text.dup() for s in mux.shapes(li).each() if s.is_text()]
    if any(not (s.is_polygon() or s.is_box() or s.is_path() or s.is_text()) for s in mux.shapes(li).each()):
        raise ValueError("MUX包含未处理的非几何对象，拒绝静默缩放")
    old = pya.Region(mux.shapes(li)).merged()
    unchanged_zone = rect([0. if revised else -100000., -100000., 750., 100000.])
    external = old-unchanged_zone
    external=old & rect([750.,-100000.,100000.,100000.])
    new_external=mapped_x_region(external,lambda x:750000+(x-750000)*external_um/200.)
    negative=old & rect([-100000.,-100000.,0.,100000.])
    if (info.layer,info.datatype)==(20,0) and revised:
        main=negative & rect([-100000.,-3.,0.,-1.])
        tail=negative-main
        new_negative=mapped_x_region(main,lambda x:x*active_um/200.)+tail
        aux_tail_checks[str(info)]=((new_negative-main.merged()) & rect([-100000.,0.,0.,100000.]) ^ tail).is_empty()
    else:
        new_negative=mapped_x_region(negative,lambda x:x*active_um/200.)
    new = ((old & rect([0.,-100000.,750.,100000.]))+new_negative+new_external).merged()
    mux_body_checks[str(info)] = ((old ^ new) & unchanged_zone).is_empty()
    mux.shapes(li).clear(); mux.shapes(li).insert(new)
    for text_item in text_items:
        x = text_item.trans.disp.x
        if x > 750000:
            text_item.trans = pya.Trans(round(750000+(x-750000)*external_um/200.)-x, 0)*text_item.trans
        elif x<0 and revised:
            text_item.trans=pya.Trans(round(x*active_um/200.)-x,0)*text_item.trans
        mux.shapes(li).insert(text_item)

# 按PDK默认波导模板修正MUX公共窗口：20/1和10/2的固定总宽
# 分别是7.2和17.2 µm。先从两条x单调核心多边形提取中心线，再取
# 两条固定宽度窗口的并集；不能把3/8 µm误当成任意芯宽的边缘外扩。
mux_window_normalization = {}
if normalize_mux_windows:
    centerlines=[]
    for shape in mux.shapes(layout.layer(20,0)).each():
        if not shape.is_polygon():continue
        core_piece=pya.Region(shape.polygon);bounds=shape.polygon.bbox();step=500
        xs=list(range(bounds.left,bounds.right+1,step))
        if xs[-1]!=bounds.right:xs.append(bounds.right)
        points=[]
        for x in xs:
            probe=min(max(x,bounds.left+1),bounds.right-1)
            section=core_piece & pya.Region(pya.Box(probe,-100000,probe+1,100000))
            if section.is_empty():raise ValueError('MUX核心截面为空，不能提取中心线')
            points.append(pya.DPoint(x*layout.dbu,section.bbox().center().y*layout.dbu))
        if len(points)<2:raise ValueError('MUX核心中心线提取点不足')
        centerlines.append(points)
    if len(centerlines)!=2:raise ValueError('MUX应有两条直接LN1核心')
    clip=rect([-100.,-1000.,800.,1000.])
    for key,total_width_um in [((20,1),7.2), ((10,2),17.2)]:
        li = layout.layer(*key)
        before = pya.Region(mux.shapes(li)).merged()
        after=pya.Region()
        for points in centerlines:after+=pya.Region(pya.DPath(points,total_width_um).to_itype(layout.dbu))
        after=(after&clip).merged()
        mux.shapes(li).clear(); mux.shapes(li).insert(after)
        mux_window_normalization[str(key)] = {
            'fixed_total_width_um_per_centerline': total_width_um,
            'before_area_um2': before.area()*layout.dbu**2,
            'after_area_um2': after.area()*layout.dbu**2,
            'changed_area_um2': (before^after).area()*layout.dbu**2,
        }
selected = []
for inst in source_top.each_inst():
    name = inst.cell.name
    if name == mux_source_name or name in (
        "T_GSG_15MM_RECTANGULAR_BASELINE",
        "DUAL_RAIL_ACTIVE_LN_PARALLEL_WINDOWS_UNVERIFIED",
        "GSG150_LAUNCH_50UM_DRAFT",
    ):
        cell = mux if name == mux_source_name else copy_cell(name)
        tr=inst.dcplx_trans
        if name==mux_source_name and revised:
            contact=tr*pya.DPoint(-200.,-1.95)
            moved=tr*pya.DPoint(-active_um,-1.95)
            tr=pya.DCplxTrans(tr)
            tr.disp=pya.DVector(tr.disp.x+contact.x-moved.x,tr.disp.y+contact.y-moved.y)
        top.insert(pya.DCellInstArray(cell.cell_index(), shift*tr))
        selected.append(name)
if len(selected)!=8:
    raise ValueError("有源/MUX/GSG实例数量不符")
body_core = reg(layout, top, (20,0))

# 所有直接掩膜来自已回标中心线；有意交叉区域已留空。
for layer in p["route_geometry"]:
    shapes = top.shapes(layout.layer(*layer["layer"]))
    for record in layer["polygons"]:
        shapes.insert(polygon(record))
route_core_direct = pya.Region(top.shapes(layout.layer(20,0))).merged()
attachment = pya.Region()
for x,y in p["attachment_points_um"]:
    attachment += rect([x-.01,y-1.,x+.01,y+1.])
route_body_overlap_outside_ports = ((route_core_direct & body_core)-attachment).area()*.001**2

cross = copy_cell("crossing_ln")
if str(cross.dbbox()) != "(0,-22.5;45,22.5)":
    raise ValueError("交叉器黑盒边界不符")
for x,y in p["crossing_centers_um"]:
    top.insert(pya.DCellInstArray(cross.cell_index(), pya.DCplxTrans(1,0,False,x-22.5,y)))
edge = copy_cell("edge_coupler_9um_y_1550_ln")
for x,y in p["edge_coupler_internal_ports_um"]:
    top.insert(pya.DCellInstArray(edge.cell_index(), pya.DCplxTrans(1,0,False,x,y)))
taper = copy_cell("LN_TAPER_070_120_L100_DRAFT")
for item in p["tapers"]:
    top.insert(pya.DCellInstArray(taper.cell_index(),pya.DCplxTrans(1,item["angle_deg"],False,*item["narrow_um"])))
core = reg(layout, top, (20,0)); clad = reg(layout,top,(20,1)); trench = reg(layout,top,(10,2))
metal = reg(layout,top,(42,0))
taper_checks=[]
for item in p["tapers"]:
    tr = pya.DCplxTrans(1,item["angle_deg"],False,*item["narrow_um"])
    flags={}
    for key,coords in [("wide_1p2",[99.999,-.6,100,.6]),("narrow_0p7",[0,-.35,.001,.35]),
                       ("narrow_join",[-.001,-.35,0,.35])]:
        patch=pya.Region(pya.DPolygon(pya.DBox(*coords)).transformed(tr).to_itype(layout.dbu))
        flags[key]=(patch-core).is_empty()
    taper_checks.append({"name":item["name"],"checks":flags})
blackbox_regions=pya.Region()
for x,y in p["crossing_centers_um"]:
    blackbox_regions+=rect([x-22.5,y-22.5,x+22.5,y+22.5])
for x,y in p["edge_coupler_internal_ports_um"]:
    bounds=edge.dbbox(); blackbox_regions+=rect([bounds.left+x,bounds.bottom+y,bounds.right+x,bounds.top+y])
# 只填平左侧新路线窗口的亚微米尖缝；不改核心、MUX本体、黑盒或YSJ。
window_cleanup=[]
left_roi=rect([-10900,-1800,float(p.get('left_cleanup_xmax_um',-10500.)),1900])
for bounds in p.get('additional_cleanup_rois',[]):left_roi+=rect(bounds)
left_roi=left_roi.merged()
for key,expected_space in [((20,1),300),((10,2),200)]:
    li=layout.layer(*key)
    original_direct=pya.Region(top.shapes(li)).merged()
    # 去除离散短边后仅圆滑凹角；在完整区域中替换，不写纳米薄差集。
    smoothing_nm=100 if key==(20,1) else 20
    rounded=original_direct.smoothed(smoothing_nm,True).rounded_corners(1000,0,64)
    replacement=((original_direct-left_roi)+(rounded&left_roi))-blackbox_regions
    replacement=replacement.merged()
    # 多条窗口汇合会留下亚微米封闭针孔；只填左侧新窗口中的小孔。
    rebuilt=pya.Region();filled_holes=[]
    for poly in replacement.each():
        clean=pya.Polygon(list(poly.each_point_hull()))
        for h in range(poly.holes()):
            hp=pya.Polygon(list(poly.each_point_hole(h)))
            hole=pya.Region(hp)
            if hp.area()<=4000000 and (hole-left_roi).is_empty() and (hole&blackbox_regions).is_empty() and (hole&core).is_empty():
                filled_holes.append({'area_um2':hp.area()*1e-6,'bbox_dbu':str(hp.bbox())})
            else:
                clean.insert_hole(list(poly.each_point_hole(h)))
        rebuilt.insert(clean)
    replacement=rebuilt.merged()
    top.shapes(li).clear();top.shapes(li).insert(replacement)
    bay_closing=[]
    merge_caps=[]
    for cap in p.get('window_merge_caps',[]):
        if tuple(cap['layer'])!=key:continue
        patch=rect(cap['bounds_um'])-blackbox_regions
        if not (patch-left_roi).is_empty():raise ValueError('窗口端帽超出路由整理范围')
        added=patch-replacement
        replacement=(replacement+patch).merged()
        merge_caps.append({'bounds_um':cap['bounds_um'],'added_area_um2':added.area()*1e-6,'reason':cap['reason']})
    if merge_caps:
        before_cap_rounding=replacement
        replacement=(replacement.rounded_corners(300,0,64)-blackbox_regions).merged()
        merge_caps.append({'inner_corner_radius_um':.3,
                           'added_area_um2':(replacement-before_cap_rounding).area()*1e-6,
                           'removed_area_um2':(before_cap_rounding-replacement).area()*1e-6})
    top.shapes(li).clear();top.shapes(li).insert(replacement)
    if revised:
        # 仅填新路由窗口的亚微米窄凹口；不移动核心，不改MUX/黑盒。
        for closing_nm in (200,400,800):
            complete=reg(layout,top,key)
            test_region=complete-core if key==(20,1) else complete
            pairs=test_region.space_check(expected_space)
            local=pya.Region()
            for ep in pairs.each():
                patch=pya.Region(ep.bbox().enlarged(1000))
                if (patch-left_roi).is_empty():local+=patch
            if local.is_empty():break
            closed=replacement.sized(closing_nm).sized(-closing_nm)
            added=((closed-replacement)&local)-blackbox_regions
            if added.is_empty():continue
            replacement=(replacement+added).merged()
            top.shapes(li).clear();top.shapes(li).insert(replacement)
            bay_closing.append({'radius_um':closing_nm*.001,'added_area_um2':added.area()*1e-6})
    updated=reg(layout,top,key)
    relevant=(updated-core).space_check(expected_space) if key==(20,1) else updated.space_check(expected_space)
    local_count=sum(1 for edge_pair in relevant.each() if edge_pair.bbox().right < -10500000)
    window_cleanup.append({'layer':list(key),'maximum_smoothing_error_um':smoothing_nm*.001,'concave_radius_um':1.,
                           'added_area_um2':(replacement-original_direct).area()*1e-6,
                           'removed_area_um2':(original_direct-replacement).area()*1e-6,
                           'filled_small_holes':filled_holes,
                           'local_bay_closing':bay_closing,
                           'flat_merge_caps':merge_caps,
                           'remaining_left_space_markers':local_count})
clad=reg(layout,top,(20,1));trench=reg(layout,top,(10,2))
bus_flags=[]
for bounds in p["bus_link_core_boxes_um"]:
    bus_flags.append((rect(bounds)-core).is_empty())
rules={"core_width_030":core.width_check(300).size(),"core_space_030":core.space_check(300).size(),
       "etch_width_030":(clad-core).width_check(300).size(),"etch_space_030":(clad-core).space_check(300).size(),
       "trench_width_020":trench.width_check(200).size(),"trench_space_020":trench.space_check(200).size(),
       "metal_width_2":metal.width_check(2000).size(),"metal_space_3":metal.space_check(3000).size(),
       "core_not_in_clad":(core-clad).size(),"core_not_in_trench":(core-trench).size()}
frame=rect([-10900,-1900,10900,1900])
checks={"MUX_body_0_to_750_unchanged":all(mux_body_checks.values()),
        "MUX_original_auxiliary_tail_unchanged":all(aux_tail_checks.values()),
        "MUX_body_windows_normalized_to_pdk":bool(normalize_mux_windows),
        "route_body_overlap_outside_ports_zero":route_body_overlap_outside_ports<1e-12,
        "taper_endpoints_joined":all(all(t["checks"].values()) for t in taper_checks),
        "serial_crossing_bus_continuous":all(bus_flags),
        "blackbox_interiors_not_drawn":all((r & blackbox_regions).is_empty() for r in [core,clad,trench]),
        "physical_masks_inside_block":all((r-frame).is_empty() for r in [core,clad,trench,metal]),
        "one_top":len(layout.top_cells())==1,
        "new_rule_markers_clear":rules=={'core_width_030':40,'core_space_030':0,'etch_width_030':0,
            'etch_space_030':40,'trench_width_020':0,'trench_space_020':0,'metal_width_2':0,
            'metal_space_3':0,'core_not_in_clad':0,'core_not_in_trench':0}}
annotations=top.shapes(layout.layer(101,0))
for text_value, x,y in [(f"15MM_GSG150_A{active_um:g}_E{external_um:g}_CYCLE_MATCHED_DRAFT",-9500,-1800),
                       ("LN2_DEFERRED_BB_LOSS_UNKNOWN",-9500,-1830)]:
    annotations.insert(pya.DText(text_value,x,y).to_itype(layout.dbu))
report={"status":"new_GDS_geometry_check_not_optical_signoff","checks":checks,"rules":rules,
        "active_taper_um":active_um,"external_taper_um":external_um,
        "window_marker_boxes_um":{name:[[v.bbox().left*.001,v.bbox().bottom*.001,v.bbox().right*.001,v.bbox().top*.001] for v in ep.each()]
            for name,ep in [('LN1_etch_space',(clad-core).space_check(300)),('SiN_trench_space',trench.space_check(200))]},
        "copy_normalizations":copy_normalizations,"mux_window_normalization":mux_window_normalization,
        "left_window_cleanup":window_cleanup,
        "taper_checks":taper_checks,"route_body_overlap_outside_ports_um2":route_body_overlap_outside_ports,
        "crossing_count":len(p["crossing_centers_um"]),"taper_count":len(p["tapers"]),
        "source_selected_instances":selected,"source_cell_signatures":{
            name:digest_cell(source,source.cell(name)) for name in copied if name!=mux_source_name},
        "no_new_LN2":all(layout.get_info(li).layer!=21 for li in layout.layer_indices()),
        "bbox_um":str(top.dbbox()),"source_gds_unchanged":hashlib.sha256(Path(p['source_gds']).read_bytes()).hexdigest()==p['source_sha256']}
options=pya.SaveLayoutOptions();options.gds2_max_vertex_count=4000;options.gds2_multi_xy_records=False
layout.write(p["ours_gds"],options)
report["new_gds_written"]=p["ours_gds"]
report["sha256"]=hashlib.sha256(Path(p["ours_gds"]).read_bytes()).hexdigest()
preview=[]
for key in [(20,0),(20,1),(10,2),(42,0),(100,0)]:
    region=reg(layout,top,key)
    preview.append({"layer":list(key),"polygons":[{"exterior":[[pt.x*.001,pt.y*.001] for pt in poly.each_point_hull()],
        "holes":[[[pt.x*.001,pt.y*.001] for pt in poly.each_point_hole(h)] for h in range(poly.holes())]} for poly in region.each()]})
with open(p["preview_json"],"x",encoding="utf-8") as f:json.dump(preview,f,ensure_ascii=False)
with open(p["report_json"],"x",encoding="utf-8") as f:json.dump(report,f,ensure_ascii=False,indent=2)
# 单独生成父单元，稍后按GDS结构记录拼接，确保YSJ所有单元记录逐字节保留。
parent_layout=pya.Layout();parent_layout.dbu=.001
parent=parent_layout.create_cell(p["merged_top"])
for name in [p["ysj_top"],p["top_name"]]:
    dummy=parent_layout.create_cell(name)
    parent.insert(pya.CellInstArray(dummy.cell_index(),pya.Trans()))
parent_layout.write(p["parent_gds"])
print(json.dumps(report,ensure_ascii=False,indent=2))
if not all(checks.values()):
    raise RuntimeError("功能图形连接检查存在未通过项；只写检查稿，不合并")
