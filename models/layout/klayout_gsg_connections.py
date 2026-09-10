"""KLayout批处理：添加GSG连接，检查短路/断路/波导覆盖及手册几何规则。"""
import json
import hashlib
import pya


with open(payload_file,'r',encoding='utf-8') as stream:payload=json.load(stream)
layout=pya.Layout();layout.read(payload['source']);top=layout.top_cell();dbu=layout.dbu
if len(layout.top_cells())!=1 or dbu!=.001:raise RuntimeError('源工程顶层或格点不符')
if top.name!='EO4P_10G_15MM_EULER_DELAY_CLOSED_DRAFT':raise RuntimeError('不是预期v16')


def signature(cell):
    rows=[(str(layout.get_info(li)),sorted(str(s) for s in cell.shapes(li).each()))
          for li in layout.layer_indices() if not cell.shapes(li).is_empty()]
    instances=sorted((i.cell.name,str(i.dcplx_trans)) for i in cell.each_inst())
    return hashlib.sha256(json.dumps([rows,instances]).encode()).hexdigest()


def region(layer,cell=top):return pya.Region(cell.begin_shapes_rec(layout.layer(*layer))).merged()


def box(coords):return pya.Region(pya.DBox(*coords).to_itype(dbu))


def markers():
    core=region((20,0));clad=region((20,1));trench=region((10,2));metal=region((42,0))
    etch=(clad-core).merged()
    return {'LN核心线宽030':core.width_check(300).size(),'LN核心间距030':core.space_check(300).size(),
            'LN1刻蚀区线宽030':etch.width_check(300).size(),'LN1刻蚀区间距030':etch.space_check(300).size(),
            'SiN挖除线宽020':trench.width_check(200).size(),'SiN挖除间距020':trench.space_check(200).size(),
            'M1线宽2um':metal.width_check(2000).size(),'M1间距3um':metal.space_check(3000).size(),
            'LN核心未覆盖SiN挖除':(core-trench).size()}


before_signatures={c.name:signature(c) for c in layout.each_cell() if c!=top}
old_instances=sorted((i.cell.name,str(i.dcplx_trans)) for i in top.each_inst())
optical_layers=[(20,0),(20,1),(10,2),(21,0),(21,1),(21,2)]
old_optics={key:region(key) for key in optical_layers}
old_metal=region((42,0));old_markers=markers()
cell=layout.create_cell(payload['cell']);config=payload['connector'];metal_layer=layout.layer(42,0)
annotation=layout.layer(101,0)
for name,points in config['polygons'].items():
    cell.shapes(metal_layer).insert(pya.DPolygon([pya.DPoint(*p) for p in points]).to_itype(dbu))
for name,(x,y) in config['contact_centers_um'].items():
    # 非工艺层上标记50x50um有效落针区及接触中心，不伪造额外开窗层。
    cell.shapes(annotation).insert(pya.DBox(x-25,y-25,x+25,y+25).to_itype(dbu))
    cell.shapes(annotation).insert(pya.DText(name+'_CONTACT',x,y).to_itype(dbu))
cell.shapes(annotation).insert(pya.DText(f'PITCH_{config["pitch_um"]:g}_GEOMETRY_ONLY',0,-config['pitch_um']-60).to_itype(dbu))

length=config['total_length_um']
transforms={'RF_IN':pya.DCplxTrans(1,0,False,2000-length,0),
            'RF_OUT':pya.DCplxTrans(1,180,True,17000+length,0)}
for name,tr in transforms.items():
    top.insert(pya.DCellInstArray(cell.cell_index(),tr))
    top.shapes(annotation).insert(pya.DText(name+'_EXTERNAL_GSG',tr.disp.x,-225).to_itype(dbu))
top.name=payload['top']
# 去掉旧的两个RF占位框，保留库黑盒实例及其所有边界。
removed_placeholders=0
expected=[pya.DBox(1840,-180,2160,180).to_itype(dbu),pya.DBox(16840,-180,17160,180).to_itype(dbu)]
for shape in list(top.shapes(layout.layer(100,0)).each()):
    if shape.is_box() and shape.box in expected:
        shape.delete();removed_placeholders+=1

metal=region((42,0));core=region((20,0));added=(metal-old_metal).merged()
conductors=list(metal.each());new_markers=markers()
new_regions=region((42,0),cell)
min_contact_clearance=None
low,high=0,20000
while high-low>1:
    mid=(low+high)//2
    if added.separation_check(core,mid).is_empty():low=mid
    else:high=mid
min_contact_clearance=low*dbu
connectivity={};contact_centers={}
for label,tr in transforms.items():
    contact_centers[label]={}
    for net,(x,y) in config['contact_centers_um'].items():
        pt=tr*pya.DPoint(x,y);contact_centers[label][net]=[pt.x,pt.y]
        landing=box([pt.x-25,pt.y-25,pt.x+25,pt.y+25])
        match=[i for i,p in enumerate(conductors) if (landing-pya.Region(p)).is_empty()]
        connectivity[label+'_'+net]={'conductor_indices':match,'landing_fully_on_metal':len(match)==1}
continuity={}
for net,y in [('S',0),('Gupper',88.5),('Glower',-88.5)]:
    i=[i for i,p in enumerate(conductors) if (box([9500-1,y-1,9500+1,y+1])-pya.Region(p)).is_empty()]
    continuity[net]=len(i)==1 and all(connectivity[label+'_'+net]['conductor_indices']==i for label in transforms)
checks={'optical_masks_unchanged':all((old_optics[k]^region(k)).is_empty() for k in optical_layers),
        'all_inherited_cells_unchanged':all(signature(layout.cell(n))==v for n,v in before_signatures.items()),
        'inherited_instances_preserved':all(item in sorted((i.cell.name,str(i.dcplx_trans)) for i in top.each_inst()) for item in old_instances),
        'active_metal_unchanged':((metal^old_metal)&box([2000,-500,17000,500])).is_empty(),
        'no_existing_metal_removed':(old_metal-metal).is_empty(),
        'three_connected_conductors':len(conductors)==3,
        'G_S_G_all_six_landings_connected':all(continuity.values()),
        'no_metal_waveguide_overlap':(added&core).is_empty(),
        'new_metal_clearance_at_least_2um_engineering_rule':min_contact_clearance>=2.,
        'M1_min_width_and_space_pass':new_markers['M1线宽2um']==0 and new_markers['M1间距3um']==0,
        'no_new_optical_rule_markers':all(new_markers[k]==v for k,v in old_markers.items() if not k.startswith('M1')),
        'both_old_RF_placeholder_boxes_replaced':removed_placeholders==2}
preview=[]
for key in [(20,0),(20,1),(10,2),(42,0)]:
    preview.append({'layer':list(key),'polygons':[{'exterior':[[p.x*dbu,p.y*dbu] for p in poly.each_point_hull()],
                    'holes':[[[p.x*dbu,p.y*dbu] for p in poly.each_point_hole(h)] for h in range(poly.holes())]}
                    for poly in region(key).each()]})
result={'status':'双端外接GSG几何工作稿；不是50ohm射频性能验证',
        'source_gds':payload['source'],'source_sha256':payload['source_sha256'],
        'output_gds':payload['output'],'top_cell':top.name,'pitch_um':config['pitch_um'],
        'connection_parameters':{k:v for k,v in config.items() if k!='polygons'},
        'probe_contact_centers_um':contact_centers,'connectivity':connectivity,
        'net_continuity':continuity,'new_metal_to_LN1_min_clearance_um':min_contact_clearance,
        'checks':checks,'gsg_geometry_pass':all(checks.values()),
        'manual_markers_before':old_markers,'manual_markers_after':new_markers,
        'inherited_cell_signatures':before_signatures,
        'optical_centerline_delay_geometry_unchanged':checks['optical_masks_unchanged'],
        'remaining':['接口S参数与阻抗未求解','实际探针型号及触针尺寸未确认',
                     '外接50ohm终端/测试链路尚未接入','原8个分支窗口规则标记及MUX已有标记',
                     '21层最终映射按用户要求暂缓','未进行制造容差或官方完整DRC']}
if all(checks.values()):
    layout.write(payload['output'])
    verify=pya.Layout();verify.read(payload['output'])
    result['gds_roundtrip_one_top']=len(verify.top_cells())==1
    result['output_sha256']=hashlib.sha256(open(payload['output'],'rb').read()).hexdigest()
else:result['output_gds_written']=False
with open(payload['report'],'w',encoding='utf-8') as stream:json.dump(result,stream,ensure_ascii=False,indent=2)
with open(payload['preview_output'],'w',encoding='utf-8') as stream:json.dump(preview,stream,ensure_ascii=False)
if not all(checks.values()):raise RuntimeError('GSG检查未通过：'+str([k for k,v in checks.items() if not v]))
