"""从v6重排四个MUX连接：包络与有源轨道平行，LN1脊在电极外使用缓慢S弯。"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from klink import KLinkClient
from shapely.geometry import LineString, Polygon

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "results/layout/尾端线宽修正_工作稿_v6/四程10GHz_15mm_尾端线宽修正_待验证工作稿.gds"


def smootherstep(u: np.ndarray) -> np.ndarray:
    return u**3 * (10.0 + u * (-15.0 + 6.0*u))


def s_taper_vertices(x0: float, x1: float, center0: float, center1: float,
                     width0: float, width1: float, samples: int = 201) -> list[list[float]]:
    if not x1 > x0 or min(width0,width1) <= 0 or samples < 3:
        raise ValueError("S形拉锥参数无效")
    x = np.linspace(x0,x1,samples)
    p = smootherstep((x-x0)/(x1-x0))
    center = center0+(center1-center0)*p
    width = width0+(width1-width0)*p
    lower = np.column_stack((x,center-width/2))
    upper = np.column_stack((x[::-1],(center+width/2)[::-1]))
    return np.vstack((lower,upper)).tolist()


def centerline_length(x0: float, x1: float, center0: float, center1: float) -> float:
    x=np.linspace(x0,x1,4001)
    p=smootherstep((x-x0)/(x1-x0))
    y=center0+(center1-center0)*p
    return float(np.sum(np.hypot(np.diff(x),np.diff(y))))


def minimum_radius(x0: float, x1: float, center0: float, center1: float) -> float:
    x=np.linspace(x0,x1,4001)
    p=smootherstep((x-x0)/(x1-x0))
    y=center0+(center1-center0)*p
    first=np.gradient(y,x)
    second=np.gradient(first,x)
    curvature=np.abs(second)/(1+first**2)**1.5
    return float(1/np.max(curvature))


def output_envelope_vertices(width_um: float, samples: int = 201) -> list[list[float]]:
    """共同包络在外部双支路分离段展开，端面匹配两条路由包络的并集。"""
    x=np.linspace(750.0,950.0,samples)
    p=smootherstep((x-750.0)/200.0)
    half=width_um/2
    lower=-half+0.15*p
    upper=half+3.40*p
    return (np.vstack(([[-200.0,-half],[750.0,-half]],
                       np.column_stack((x[1:],lower[1:])),
                       np.column_stack((x[::-1],upper[::-1])),
                       [[-200.0,half]]))).tolist()


ACTIVE = s_taper_vertices(-200,0,0,-1.95,1.33,1.40)
EXT_MAIN = s_taper_vertices(750,950,-1.80,.15,1.10,.70)
EXT_AUX = s_taper_vertices(750,950,1.45,3.40,.40,.70)
CLAD_ENVELOPE = output_envelope_vertices(7.2)
TRENCH_ENVELOPE = output_envelope_vertices(17.2)

KLAYOUT_CODE = r'''
import hashlib
import json

def signature(layout,cell):
    shapes=[(str(layout.get_info(li)),sorted(str(s) for s in cell.shapes(li).each()))
            for li in layout.layer_indices() if not cell.shapes(li).is_empty()]
    insts=sorted((i.cell.name,str(i.dcplx_trans)) for i in cell.each_inst())
    return hashlib.sha256(json.dumps([sorted(shapes),insts]).encode()).hexdigest()

layout=pya.Layout()
layout.read(payload['input'])
top=layout.top_cell()
if top is None or layout.dbu!=.001 or len(layout.top_cells())!=1:
    raise RuntimeError('v6顶层或格点不匹配')
muxes=[c for c in layout.each_cell() if c.name.startswith('MUX_TE01_A70_TAIL030')]
if len(muxes)!=1:
    raise RuntimeError('未唯一找到v6 MUX')
mux=muxes[0]
instances=[i for i in top.each_inst() if i.cell==mux]
if len(instances)!=4:
    raise RuntimeError('MUX实例数不是四个')
li_core=layout.layer(20,0)
original_mux_layers={str(layout.get_info(li)):pya.Region(mux.shapes(li)).merged() for li in layout.layer_indices()}
other_cells={c.name:signature(layout,c) for c in layout.each_cell() if c not in (top,mux)}
other_top_layers={str(layout.get_info(li)):pya.Region(top.begin_shapes_rec(li)).merged()
                  for li in layout.layer_indices() if (layout.get_info(li).layer,layout.get_info(li).datatype)
                  not in ((20,0),(20,1),(10,2),(101,0))}
core_shapes=list(mux.shapes(li_core).each())
bboxes={'active':pya.Box(-200000,-2650,0,-1250),
        'main':pya.Box(750000,-2350,950000,-1250),
        'aux':pya.Box(750000,1100,950000,1800)}
selected={name:[s for s in core_shapes if s.bbox()==bbox] for name,bbox in bboxes.items()}
if any(len(v)!=1 for v in selected.values()):
    raise RuntimeError('无法唯一识别三段原直拉锥')
for shape in selected.values():
    shape[0].delete()
for name in ('active','main','aux'):
    polygon=pya.DPolygon([pya.DPoint(*xy) for xy in payload['vertices'][name]])
    mux.shapes(li_core).insert(polygon.to_itype(layout.dbu))
for layer_key in ((20,1),(10,2)):
    li=layout.layer(*layer_key)
    shapes=list(mux.shapes(li).each())
    width=payload['envelope_widths'][str(layer_key)]
    expected=pya.Box(-200000,-round(width*500),950000,round(width*500))
    matching=[s for s in shapes if s.bbox()==expected]
    if len(matching)!=1 or len(shapes)!=1:
        raise RuntimeError('无法唯一识别原共同包络：'+str(layer_key))
    matching[0].delete()
    polygon=pya.DPolygon([pya.DPoint(*xy) for xy in payload['envelopes'][str(layer_key)]])
    mux.shapes(li).insert(polygon.to_itype(layout.dbu))

old_transforms=[]
new_transforms=[]
for inst in instances:
    tr=inst.dcplx_trans
    old_transforms.append(str(tr))
    target_y=30.0 if tr.disp.y>0 else -30.0
    tr.disp=pya.DVector(tr.disp.x,target_y)
    inst.dcplx_trans=tr
    new_transforms.append(str(tr))
mux.name='MUX_TE01_A70_TAIL030_CENTERED_SPORTS_UNVERIFIED'
top.name='EO4P_10G_15MM_MUX_REPLANNED_DRAFT'
mux.shapes(layout.layer(101,0)).insert(pya.Text('MUX_CENTERED_ENVELOPE_S_BEND_PORTS_UNVERIFIED',pya.Trans(-190000,-6800)))

checks={
 'four_mux_instances':len(instances)==4,
 'other_cells_unchanged':all(signature(layout,layout.cell(n))==h for n,h in other_cells.items()),
 'unrelated_process_layers_unchanged_inside_mux':all((pya.Region(mux.shapes(li)).merged()^old).is_empty()
      for li in layout.layer_indices() if str(layout.get_info(li)) in original_mux_layers
      and (layout.get_info(li).layer,layout.get_info(li).datatype) not in ((20,0),(20,1),(10,2),(101,0))
      for old in [original_mux_layers[str(layout.get_info(li))]]),
 'noncore_top_layers_unchanged':all((pya.Region(top.begin_shapes_rec(li)).merged()^old).is_empty()
      for li in layout.layer_indices() if str(layout.get_info(li)) in other_top_layers
      for old in [other_top_layers[str(layout.get_info(li))]]),
}
core=pya.Region(top.begin_shapes_rec(li_core)).merged()
clad=pya.Region(top.begin_shapes_rec(layout.layer(20,1))).merged()
trench=pya.Region(top.begin_shapes_rec(layout.layer(10,2))).merged()
metal=pya.Region(top.begin_shapes_rec(layout.layer(42,0))).merged()
joint_rows=[]
for index,inst in enumerate(instances,1):
    row={'index':index,'transform':str(inst.dcplx_trans)}
    for name,x,y,width in [('active',-200,0,1.33),('external_main',950,.15,.7),('external_aux',950,3.4,.7)]:
        test=pya.Region(pya.DPolygon(pya.DBox(x-.001,y-width/2,x+.001,y+width/2)).transformed(inst.dcplx_trans).to_itype(layout.dbu))
        ok=(test-core).is_empty()
        point=inst.dcplx_trans*pya.DPoint(x,y)
        row[name]={'global_um':[point.x,point.y],'core_connected':ok}
        checks['join_'+str(index)+'_'+name]=ok
    for name,layer,width in [('clad',layout.layer(20,1),7.2),('trench',layout.layer(10,2),17.2)]:
        test=pya.Region(pya.DPolygon(pya.DBox(-200.001,-width/2,-199.999,width/2)).transformed(inst.dcplx_trans).to_itype(layout.dbu))
        ok=(test-pya.Region(top.begin_shapes_rec(layer)).merged()).is_empty()
        row['active_'+name+'_parallel_connected']=ok
        checks['active_'+name+'_'+str(index)]=ok
        for branch,y in [('main',.15),('aux',3.4)]:
            output_test=pya.Region(pya.DPolygon(pya.DBox(949.999,y-width/2,950.001,y+width/2)).transformed(inst.dcplx_trans).to_itype(layout.dbu))
            output_ok=(output_test-pya.Region(top.begin_shapes_rec(layer)).merged()).is_empty()
            row['external_'+branch+'_'+name+'_full_width_connected']=output_ok
            checks['external_'+branch+'_'+name+'_'+str(index)]=output_ok
    joint_rows.append(row)
checks['no_metal_core_overlap']=(metal&core).is_empty()
checks['core_space_030']=core.space_check(300).is_empty()
mux_core=pya.Region(mux.shapes(li_core)).merged()
mux_clad=pya.Region(mux.shapes(layout.layer(20,1))).merged()
mux_trench=pya.Region(mux.shapes(layout.layer(10,2))).merged()
checks['mux_clad_covers_mux_core']=(mux_core-mux_clad).is_empty()
checks['mux_trench_covers_mux_core']=(mux_core-mux_trench).is_empty()
if not all(checks.values()):
    raise RuntimeError('重新规划连接检查失败：'+json.dumps(checks))
etch1=clad-core
markers={'core_width_030':core.width_check(300).size(),'core_space_030':core.space_check(300).size(),
         'etch1_width_030':etch1.width_check(300).size(),'etch1_space_030':etch1.space_check(300).size(),
         'sin_trench_width_020':trench.width_check(200).size(),'sin_trench_space_020':trench.space_check(200).size(),
         'metal_width_2':metal.width_check(2000).size(),'metal_space_3':metal.space_check(3000).size()}
layout.write(payload['output'])
{'top_cell':top.name,'mux_cell':mux.name,'old_transforms':old_transforms,'new_transforms':new_transforms,
 'instances_repositioned':4,'joint_checks':joint_rows,'checks':checks,'manual_rule_markers':markers,
 'envelope_design':'有源侧至MUX本体保持平行共中心；仅外部双支路分离段展开并匹配两条路由包络',
 'core_design':'有源主脊及两个外部端口在各自200um无电极段内采用五次平滑S形拉锥',
 'optical_validation_pass':False}
'''


def main() -> None:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',type=Path,default=BASE)
    p.add_argument('--output-dir',type=Path,required=True)
    a=p.parse_args()
    vertices={'active':ACTIVE,'main':EXT_MAIN,'aux':EXT_AUX}
    envelopes={'(20, 1)':CLAD_ENVELOPE,'(10, 2)':TRENCH_ENVELOPE}
    specifications={'active':(-200,0,0,-1.95,1.33,1.40),
                    'main':(750,950,-1.80,.15,1.10,.70),
                    'aux':(750,950,1.45,3.40,.40,.70)}
    numerical={}
    for name,(x0,x1,c0,c1,w0,w1) in specifications.items():
        shape=Polygon(vertices[name])
        if not shape.is_valid:
            raise ValueError('S形拉锥无效：'+name)
        for x,c,w in ((x0,c0,w0),(x1,c1,w1)):
            section=shape.intersection(LineString([(x,-10),(x,10)]))
            if abs(section.length-w)>1e-8 or abs((section.bounds[1]+section.bounds[3])/2-c)>1e-8:
                raise ValueError('端口尺寸不正确：'+name)
        length=centerline_length(x0,x1,c0,c1)
        radius=minimum_radius(x0,x1,c0,c1)
        numerical[name]={'centerline_length_um':length,'extra_length_vs_200um':length-(x1-x0),
                         'minimum_centerline_radius_um':radius,'pdk_radius_80um_pass':radius>=80}
        if radius<80:
            raise ValueError('S形中心线最小曲率半径小于80um：'+name)
    for width,key in ((7.2,'(20, 1)'),(17.2,'(10, 2)')):
        shape=Polygon(envelopes[key])
        if not shape.is_valid:
            raise ValueError('外部展开包络无效：'+key)
        expected={-200:(-width/2,width/2),750:(-width/2,width/2),
                  950:(.15-width/2,3.4+width/2)}
        for x,(low,high) in expected.items():
            cut=shape.intersection(LineString([(x,-20),(x,20)]))
            if abs(cut.bounds[1]-low)>1e-8 or abs(cut.bounds[3]-high)>1e-8:
                raise ValueError('外部包络端面不匹配：'+key)
    a.output_dir.mkdir(parents=True,exist_ok=False)
    output=(a.output_dir/'四程10GHz_15mm_MUX重新规划平行接续_待验证工作稿.gds').resolve()
    source_hash=hashlib.sha256(a.input.read_bytes()).hexdigest()
    payload={'input':str(a.input.resolve()),'output':str(output),'vertices':vertices,
             'envelopes':envelopes,'envelope_widths':{'(20, 1)':7.2,'(10, 2)':17.2}}
    c=KLinkClient()
    c.connect()
    try:
        response=c.call('exec.python',{'code':'payload = '+repr(payload)+'\n'+KLAYOUT_CODE})
        if response.get('exception'):
            raise RuntimeError(str(response['exception']))
        report=response['return_value']
        report['file_info']=c.call('layout.file_info',{'path':str(output),'detail':'counts'})
        c.call('layout.show_file',{'path':str(output),'mode':'new'})
        c.call('view.show_cell',{'cell':report['top_cell'],'zoom_fit':True})
        for layer,color in [('20/0','#C62858'),('20/1','#339A90'),('10/2','#65A7CF'),('42/0','#D5A12C')]:
            c.call('layer.set_style',{'layer':layer,'fill_color':color,'frame_color':color,
                   'dither_pattern':0 if layer in ('20/0','42/0') else 1,'line_width':1})
        c.call('layer.set_visible',{'layers':['101/0','100/0','100/30','56/30'],'visible':False})
        c.call('view.hier_levels',{'min':0,'max':10})
        for name,bbox in [('左侧平行接续与S弯.png',[1725,-43,2075,43]),
                          ('右侧平行接续与S弯.png',[16925,-43,17275,43])]:
            c.call('view.screenshot',{'mode':'path','path':str((a.output_dir/name).resolve()),
                  'bbox_um':bbox,'width_px':1900,'height_px':620})
        c.call('view.zoom_box',{'bbox_um':[1725,20,2075,42]})
    finally:
        c.close()
    if hashlib.sha256(a.input.read_bytes()).hexdigest()!=source_hash:
        raise RuntimeError('v6源文件发生变化')
    report.update({'working_gds':str(output),'source_gds':str(a.input.resolve()),
                   'source_sha256':source_hash,'source_unchanged':True,'numerical_geometry':numerical,
                   'centerline_long_routes_unchanged':True,'foundry_signoff':False,
                   'remaining':['三段新S形拉锥的模式传输、反射和群时延尚未仿真',
                                '外部双支路的共同刻蚀包络接续仍需结合最终截面核对',
                                'MUX本体辅助脊仍有0.3um法向宽度/格点标记',
                                '未处理21层、圆弯、RF焊盘和终端；继承v6其余限制']})
    (a.output_dir/'MUX重新规划接续检查.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k!='file_info'},ensure_ascii=False,indent=2))


if __name__=='__main__':
    main()
