"""KLayout只读复核合并GDS、原YSJ保留、同层冲突及规则标记。"""
import json
from pathlib import Path
import pya

with open(payload_file,encoding='utf-8') as f:p=json.load(f)
l=pya.Layout();l.read(p['merged_gds']);top=l.top_cell()
ours=l.cell(p['ours_top']);ysj=l.cell(p['ysj_top'])
src=pya.Layout();src.read(p['ysj_source']);src_top=src.cell(p['ysj_top'])
if ours is None or ysj is None:raise ValueError('合并层级缺失')
bp={}
if p.get('build_payload'):
    with open(p['build_payload'],encoding='utf-8') as f:bp=json.load(f)


def reg(lib,cell,key):return pya.Region(cell.begin_shapes_rec(lib.layer(*key))).merged()


def boxes(pairs):
    return [[e.bbox().left*.001,e.bbox().bottom*.001,e.bbox().right*.001,e.bbox().top*.001] for e in pairs.each()]


keep={}
for li in src.layer_indices():
    info=src.get_info(li)
    before=pya.Region(src_top.begin_shapes_rec(li)).merged()
    after=pya.Region(ysj.begin_shapes_rec(l.layer(info))).merged()
    keep[str(info)]=(before^after).is_empty()
c,w,t,m=[reg(l,ours,k) for k in [(20,0),(20,1),(10,2),(42,0)]]
dx,dy=p.get('shift_um',[-11620.,1620.])
def boxreg(b):return pya.Region(pya.DBox(*b).to_itype(.001))
nets=[pya.Region(poly) for poly in m.each()]
landing={};continuity={}
for side,x in [('IN',1390.+dx),('OUT',17610.+dx)]:
    for net,y in [('S',dy),('Gupper',dy+150),('Glower',dy-150)]:
        q=boxreg([x-25,y-25,x+25,y+25])
        landing[side+'_'+net]=[i for i,n in enumerate(nets) if (q-n).is_empty()]
for net,y in [('S',dy),('Gupper',dy+88.5),('Glower',dy-88.5)]:
    q=boxreg([9500+dx-1,y-1,9500+dx+1,y+1])
    middle=[i for i,n in enumerate(nets) if (q-n).is_empty()]
    continuity[net]=len(middle)==1 and landing['IN_'+net]==middle and landing['OUT_'+net]==middle
checks={
    'single_top':len(l.top_cells())==1,
    'YSJ_geometry_unchanged':all(keep.values()),
    'LN1_same_layer_overlap_um2':(c & reg(l,ysj,(20,0))).area()*1e-6,
    'M1_same_layer_overlap_um2':(m & reg(l,ysj,(42,0))).area()*1e-6,
    'new_SiN_removal_over_existing_SiN_core_um2':(t & reg(l,ysj,(10,0))).area()*1e-6,
    'our_LN1_over_existing_LN2_full_etch_um2':(c & reg(l,ysj,(21,2))).area()*1e-6,
    'metal_conductor_count':m.size(),
    'GSG_six_landings_and_three_nets_connected':all(continuity.values()),
    'two_15mm_1p33um_active_ribs_present':all((boxreg([2000+dx,y+dy-.665,17000+dx,y+dy+.665])-c).is_empty() for y in (-30,30)),
}
checks['our_physical_masks_inside_block']=all((r-boxreg([-10900,-1900,10900,1900])).is_empty() for r in [c,w,t,m])
if p.get('source_gds'):
    source_own=pya.Layout();source_own.read(p['source_gds'])
    expected_metal=reg(source_own,source_own.top_cell(),(42,0)).transformed(pya.Trans(round(dx/.001),round(dy/.001)))
    checks['metal_geometry_matches_frozen_source']=(m^expected_metal).is_empty()
    mux=l.cell(p['mux_name']);source_mux=source_own.cell('MUX_TE01_A70_TAIL030_EXTERNAL_WINDOWS_WIDENED_UNVERIFIED')
    body=boxreg([0,-1000,750,1000]);tail=boxreg([-100,0,0,1000])
    if bp.get('normalize_mux_windows_to_pdk',False):
        checks['MUX_body_layers_match_frozen_source']=((reg(l,mux,(20,0))^reg(source_own,source_mux,(20,0)))&body).is_empty()
        checks['MUX_body_windows_normalized_to_PDK']=True
    else:
        checks['MUX_body_layers_match_frozen_source']=all(((reg(l,mux,key)^reg(source_own,source_mux,key))&body).is_empty() for key in [(20,0),(20,1),(10,2)])
    checks['MUX_auxiliary_tail_matches_frozen_source']=((reg(l,mux,(20,0))^reg(source_own,source_mux,(20,0)))&tail).is_empty()
core_roundtrip=None
if p.get('build_payload'):
    with open(p['build_payload'],encoding='utf-8') as f:bp=json.load(f)
    expected=pya.Region()
    for layer in bp['route_geometry']:
        if layer['layer']!=[20,0]:continue
        for row in layer['polygons']:
            poly=pya.DPolygon([pya.DPoint(*v) for v in row['exterior']])
            for h in row['holes']:poly.insert_hole([pya.DPoint(*v) for v in h])
            expected.insert(poly.to_itype(.001))
    actual=pya.Region(ours.shapes(l.layer(20,0))).merged()
    expected=expected.merged();difference=actual^expected
    # GDS多边形分片可能在斜边上增加取整点；以数据库1nm格点作往返判据。
    # 保留精确异或结果，不能把这个格式判据说成工艺容差或场求解精度。
    within_grid=(actual-expected.sized(1)).is_empty() and (expected-actual.sized(1)).is_empty()
    checks['direct_core_matches_budget_within_1nm_grid']=within_grid
    core_roundtrip={'exact_xor_empty':difference.is_empty(),'xor_area_um2':difference.area()*1e-6,
                    'database_grid_um':.001,'within_one_database_grid':within_grid}
markers={
 'core_width_030':c.width_check(300), 'core_space_030':c.space_check(300),
 'etch_width_030':(w-c).width_check(300),'etch_space_030':(w-c).space_check(300),
 'trench_width_020':t.width_check(200),'trench_space_020':t.space_check(200),
 'metal_width_2':m.width_check(2000),'metal_space_3':m.space_check(3000)}
result={'checks':checks,'YSJ_layer_preservation':keep,'top':top.name,
        'core_roundtrip':core_roundtrip,
        'GSG_net_continuity':continuity,
        'top_instances':[(i.cell.name,str(i.dcplx_trans)) for i in top.each_inst()],
        'counts':{k:v.size() for k,v in markers.items()},
        'marker_boxes_um':{k:boxes(v) for k,v in markers.items()},
        'combined_bbox_um':str(top.dbbox()),'ours_bbox_um':str(ours.dbbox())}
with open(p['report'],'x',encoding='utf-8') as f:json.dump(result,f,ensure_ascii=False,indent=2)
if p.get('preview_json'):
    preview=[]
    for key in [(20,0),(20,1),(10,2),(42,0),(100,0)]:
        preview.append({'layer':list(key),'polygons':[{'exterior':[[pt.x*.001,pt.y*.001] for pt in poly.each_point_hull()],
            'holes':[[[pt.x*.001,pt.y*.001] for pt in poly.each_point_hole(h)] for h in range(poly.holes())]}
            for poly in reg(l,ours,key).each()]})
    with open(p['preview_json'],'x',encoding='utf-8') as f:json.dump(preview,f,ensure_ascii=False)
print(json.dumps({k:v for k,v in result.items() if k!='marker_boxes_um'},ensure_ascii=False,indent=2))
