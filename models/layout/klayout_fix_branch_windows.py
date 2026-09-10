"""KLayout原生布尔检查局部分叉窗口修正，不改光学核心或电极。"""
import json
import hashlib
import pya
with open(payload_file,encoding='utf-8') as f:p=json.load(f)
l=pya.Layout();l.read(p['source']);top=l.top_cell();dbu=l.dbu


def reg(layer):return pya.Region(top.begin_shapes_rec(l.layer(*layer))).merged()


def checks():
    c=reg((20,0));e=(reg((20,1))-c).merged();t=reg((10,2));m=reg((42,0))
    return {'core_width_030':c.width_check(300).size(),'core_space_030':c.space_check(300).size(),
        'etch1_width_030':e.width_check(300).size(),'etch1_space_030':e.space_check(300).size(),
        'sin_width_020':t.width_check(200).size(),'sin_space_020':t.space_check(200).size(),
        'metal_width_2':m.width_check(2000).size(),'metal_space_3':m.space_check(3000).size(),
        'core_not_in_clad':(c-reg((20,1))).size(),'core_not_in_trench':(c-t).size()}


original={(20,0):reg((20,0)),(42,0):reg((42,0)),(20,1):reg((20,1)),(10,2):reg((10,2))}
before=checks();source_instances=sorted((i.cell.name,str(i.dcplx_trans)) for i in top.each_inst())
changed={}
for layer,boxes in p['patches'].items():
    key=tuple(map(int,layer.split('/')));patch=pya.Region()
    for coords in boxes:
        b=pya.DBox(*coords).to_itype(dbu);top.shapes(l.layer(*key)).insert(b);patch.insert(b)
    changed[layer]=(reg(key)-original[key]).area()*dbu**2
    # 只填充分叉窗口局部的凹角，原生整数几何圆角避免尖缝。
    current=pya.Region(top.shapes(l.layer(*key))).merged()
    # 先合并曲线离散中的短边，避免圆角半径被0.5um采样小边限制。
    rounded=current.smoothed(20,True).rounded_corners(1000,0,64)
    local=pya.Region()
    windows=[[690,20,860,230],[690,-230,860,-20]]
    if key==(10,2):windows += [[18140,16,18320,240],[18140,-240,18320,-16]]
    for coords in windows:
        local.insert(pya.DBox(*coords).to_itype(dbu))
    # 直接替换局部完整边界，避免把圆角差集拆成接近1nm的独立碎片。
    replacement=((current-local)+(rounded&local)).merged()
    top.shapes(l.layer(*key)).clear()
    top.shapes(l.layer(*key)).insert(replacement)
    changed[layer]=(reg(key)-original[key]).area()*dbu**2
after=checks()
flags={'LN1_core_unchanged':(reg((20,0))^original[(20,0)]).is_empty(),
       'M1_unchanged':(reg((42,0))^original[(42,0)]).is_empty(),
       'all_instances_unchanged':source_instances==sorted((i.cell.name,str(i.dcplx_trans)) for i in top.each_inst()),
       'original_MUX_40_width_and_40_etch_space_only':after['core_width_030']==40 and after['etch1_space_030']==40,
       'no_new_width_or_space':all(after[k]==0 for k in ['core_space_030','etch1_width_030','sin_width_020','sin_space_020','metal_width_2','metal_space_3']),
       'core_covered':after['core_not_in_clad']==0 and after['core_not_in_trench']==0}
result={'status':'局部公共刻蚀窗口修正；名义平面几何检查','source_report':p['source_report'],
        'source_gds':p['source'],'source_sha256':p['source_sha256'],
        'output_gds':p['output'],'top_cell':p['top'],'patch_boxes_um':p['patches'],
        'added_area_um2':changed,'before':before,'after':after,'checks':flags,'geometry_pass':all(flags.values()),
        'delay_note':'光学核心中心线长度未改变，仍沿用300/350/300ps名义预算；窗口介质变化未独立求模',
        'remaining':['MUX本体40个线宽及40个刻蚀间距标记','21层最终映射','HFSS接口S参数','弯曲损耗与多模验证']}
result['boundary_smoothing']={'maximum_simplification_error_um':.02,'keep_horizontal_vertical':True,
 'concave_radius_um':1.,'convex_radius_um':0.,'corner_points_per_circle':64,
 'scope':'仅20/1与10/2局部端口窗口；不改变20/0脊'}
e=(reg((20,1))-reg((20,0))).merged();t=reg((10,2))
result['remaining_branch_edge_pairs_dbu']={
 'etch1':[str(ep) for ep in e.space_check(300).each() if ep.bbox().left<1000000 or ep.bbox().left>18100000],
 'sin':[str(ep) for ep in t.space_check(200).each()]}
if all(flags.values()):
    top.name=p['top']
    options=pya.SaveLayoutOptions();options.gds2_max_vertex_count=4000;options.gds2_multi_xy_records=False
    l.write(p['output'],options)
    reread=pya.Layout();reread.read(p['output']);reread_top=reread.top_cell()
    equality={};within_grid={};xor_area={}
    for key in [(20,0),(20,1),(10,2),(42,0)]:
        expected=reg(key);read=pya.Region(reread_top.begin_shapes_rec(reread.layer(*key))).merged()
        equality[str(key)]=(expected^read).is_empty()
        within_grid[str(key)]=(expected-read.sized(1)).is_empty() and (read-expected.sized(1)).is_empty()
        xor_area[str(key)]=(expected^read).area()*dbu**2
    result['gds_roundtrip_exact']=equality
    result['gds_roundtrip_within_1nm_grid']=within_grid
    result['gds_roundtrip_xor_area_um2']=xor_area
    result['gds_max_vertex_count']=4000
    previous_layout,previous_top=l,top
    l,top=reread,reread_top
    result['actual_saved_gds_markers']=checks()
    l,top=previous_layout,previous_top
    result['saved_markers_equal_design']=result['actual_saved_gds_markers']==after
    if not all(within_grid.values()) or not result['saved_markers_equal_design']:
        with open(p['report'],'w',encoding='utf-8') as f:json.dump(result,f,ensure_ascii=False,indent=2)
        raise RuntimeError('GDS回读超出1nm格点误差，或产生新规则标记')
    with open(p['output'],'rb') as f:result['output_sha256']=hashlib.sha256(f.read()).hexdigest()
    preview=[]
    for key in [(20,0),(20,1),(10,2),(42,0)]:
        preview.append({'layer':list(key),'polygons':[{'exterior':[[pt.x*dbu,pt.y*dbu] for pt in poly.each_point_hull()],
            'holes':[[[pt.x*dbu,pt.y*dbu] for pt in poly.each_point_hole(h)] for h in range(poly.holes())]} for poly in reg(key).each()]})
    with open(p['preview'],'w',encoding='utf-8') as f:json.dump(preview,f,ensure_ascii=False)
with open(p['report'],'w',encoding='utf-8') as f:json.dump(result,f,ensure_ascii=False,indent=2)
if not all(flags.values()):raise RuntimeError(str(after))
