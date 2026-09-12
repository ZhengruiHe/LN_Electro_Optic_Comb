"""从v6生成直LN1脊连接审阅版；不改变工艺几何，只验证、另存并显示。"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from klink import KLinkClient

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "results/layout/_历史归档/迭代版本_20260912/尾端线宽修正_工作稿_v6/四程10GHz_15mm_尾端线宽修正_待验证工作稿.gds"

CODE = r'''
import json
source=pya.Layout()
source.read(payload['input'])
layout=pya.Layout()
layout.read(payload['input'])
source_top=source.top_cell()
top=layout.top_cell()
if top is None or layout.dbu!=.001 or len(layout.top_cells())!=1:
    raise RuntimeError('v6顶层或1nm格点不匹配')
muxes=[c for c in layout.each_cell() if c.name.startswith('MUX_TE01_A70_TAIL030')]
if len(muxes)!=1:
    raise RuntimeError('未唯一找到v6直连接MUX')
mux=muxes[0]
instances=[i for i in top.each_inst() if i.cell==mux]
if len(instances)!=4:
    raise RuntimeError('MUX实例不是四个')
li_core=layout.layer(20,0)
local_core=pya.Region(mux.shapes(li_core)).merged()
segments=[('active',-200,0,-1.95,1.33,1.40,pya.Box(-200000,-2650,0,-1250)),
          ('external_main',750,950,-1.80,1.10,.70,pya.Box(750000,-2350,950000,-1250)),
          ('external_aux',750,950,1.45,.40,.70,pya.Box(750000,1100,950000,1800))]
local_checks={}
core_shapes=list(mux.shapes(li_core).each())
for name,x0,x1,center,w0,w1,bbox in segments:
    candidates=[s for s in core_shapes if s.bbox()==bbox]
    if len(candidates)!=1:
        raise RuntimeError('未唯一识别直拉锥：'+name)
    polygon=candidates[0].polygon
    sections=[]
    for x,w in ((x0,w0),(x1,w1)):
        x_dbu=round(x/layout.dbu)
        values=[point.y*layout.dbu for point in polygon.each_point_hull() if point.x==x_dbu]
        if len(values)!=2:
            raise RuntimeError('拉锥端点顶点数量不为2：'+name)
        measured_center=(min(values)+max(values))/2
        measured_width=max(values)-min(values)
        sections.append({'x_um':x,'expected_width_um':w,'measured_center_um':measured_center,
                         'measured_width_um':measured_width})
    local_checks[name]={'straight_center_um':center,'end_sections':sections,
                        'centers_match':all(abs(v['measured_center_um']-center)<=.001 for v in sections),
                        'widths_match':all(abs(v['measured_width_um']-v['expected_width_um'])<=.002 for v in sections)}
checks={'four_mux_instances':len(instances)==4,
        'all_three_connections_have_constant_center':all(v['centers_match'] for v in local_checks.values()),
        'all_endpoint_widths_match':all(v['widths_match'] for v in local_checks.values())}
core=pya.Region(top.begin_shapes_rec(li_core)).merged()
joints=[]
for index,inst in enumerate(instances,1):
    row={'index':index,'transform':str(inst.dcplx_trans)}
    for name,x,y,w in [('active',-200,-1.95,1.33),('external_main',950,-1.8,.7),('external_aux',950,1.45,.7)]:
        test=pya.Region(pya.DPolygon(pya.DBox(x-.001,y-w/2,x+.001,y+w/2)).transformed(inst.dcplx_trans).to_itype(layout.dbu))
        ok=(test-core).is_empty()
        point=inst.dcplx_trans*pya.DPoint(x,y)
        row[name]={'global_um':[point.x,point.y],'connected':ok}
        checks['joint_'+str(index)+'_'+name]=ok
    joints.append(row)
process_layers=[]
for li in layout.layer_indices():
    info=layout.get_info(li)
    if (info.layer,info.datatype)==(101,0):
        continue
    source_li=source.layer(info.layer,info.datatype)
    unchanged=(pya.Region(top.begin_shapes_rec(li)).merged()^
               pya.Region(source_top.begin_shapes_rec(source_li)).merged()).is_empty()
    process_layers.append({'layer':str(info),'unchanged':unchanged})
checks['all_process_geometries_unchanged']=all(row['unchanged'] for row in process_layers)
checks['all_instance_transforms_unchanged']=(
    sorted((i.cell.name,str(i.dcplx_trans)) for i in source_top.each_inst())==
    sorted((i.cell.name,str(i.dcplx_trans)) for i in top.each_inst()))
if not all(checks.values()):
    raise RuntimeError('直脊规划核查失败：'+json.dumps(checks))
old_top=top.name
old_mux=mux.name
top.name='EO4P_10G_15MM_STRAIGHT_RIB_PLAN_DRAFT'
mux.name='MUX_TE01_A70_TAIL030_STRAIGHT_PORTS_UNVERIFIED'
mux.shapes(layout.layer(101,0)).insert(pya.Text('STRAIGHT_RIB_PORTS_MASK_ENVELOPE_PENDING',pya.Trans(-190000,-6800)))
layout.write(payload['output'])
{'status':'直LN1脊连接审阅版；未优化圆弯、时延或工艺包络','top_cell':top.name,'mux_cell':mux.name,
 'source_top_cell':old_top,'source_mux_cell':old_mux,'local_straight_checks':local_checks,
 'joint_checks':joints,'checks':checks,'process_geometry_changed':False,'optical_validation_pass':False,
 'deferred':['20/1与10/2端面台阶和最终截面映射','所有圆弯与10GHz时延回标','RF焊盘和终端']}
'''


def main() -> None:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,default=BASE)
    parser.add_argument('--output-dir',type=Path,required=True)
    args=parser.parse_args()
    args.output_dir.mkdir(parents=True,exist_ok=False)
    output=(args.output_dir/'四程10GHz_15mm_直LN1脊连接规划_待验证工作稿.gds').resolve()
    before=hashlib.sha256(args.input.read_bytes()).hexdigest()
    client=KLinkClient()
    client.connect()
    try:
        response=client.call('exec.python',{'code':'payload = '+repr({'input':str(args.input.resolve()),'output':str(output)})+'\n'+CODE})
        if response.get('exception'):
            raise RuntimeError(str(response['exception']))
        report=response['return_value']
        report['file_info']=client.call('layout.file_info',{'path':str(output),'detail':'counts'})
        client.call('layout.show_file',{'path':str(output),'mode':'new'})
        client.call('view.show_cell',{'cell':report['top_cell'],'zoom_fit':True})
        client.call('layer.set_style',{'layer':'20/0','fill_color':'#C62858','frame_color':'#C62858',
                                      'dither_pattern':0,'line_width':2})
        client.call('layer.set_visible',{'layers':['10/2','20/1','42/0','100/0','100/30','56/30','101/0'],'visible':False})
        client.call('view.hier_levels',{'min':0,'max':10})
        for name,bbox,w,h in [('左侧直脊连接.png',[1725,-43,2075,43],1900,620),
                              ('右侧直脊连接.png',[16925,-43,17275,43],1900,620),
                              ('直LN1脊全图.png',[0,-530,18850,905],2400,600)]:
            client.call('view.screenshot',{'mode':'path','path':str((args.output_dir/name).resolve()),
                                          'bbox_um':bbox,'width_px':w,'height_px':h})
        client.call('view.zoom_box',{'bbox_um':[1725,20,2075,42]})
    finally:
        client.close()
    if hashlib.sha256(args.input.read_bytes()).hexdigest()!=before:
        raise RuntimeError('v6源文件发生变化')
    report.update({'working_gds':str(output),'source_gds':str(args.input.resolve()),
                   'source_sha256':before,'source_unchanged':True,'foundry_signoff':False,
                   'note':'v8的S形端口方案已否决；本文件恢复v6直脊，只作为连接与端口规划基线。'})
    (args.output_dir/'直脊连接规划检查.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k!='file_info'},ensure_ascii=False,indent=2))


if __name__=='__main__':
    main()
