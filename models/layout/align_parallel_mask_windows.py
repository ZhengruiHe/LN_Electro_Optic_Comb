"""在v9直脊基线上平行对齐四个有源端20/1与10/2窗口，不改LN1脊、圆弯或时延。"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from klink import KLinkClient

ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/'results/layout/直LN1脊连接规划_工作稿_v9_01/四程10GHz_15mm_直LN1脊连接规划_待验证工作稿.gds'

CODE=r'''
import hashlib
import json
def signature(layout,cell,skip_text=False):
    rows=[]
    for li in layout.layer_indices():
        if skip_text and (layout.get_info(li).layer,layout.get_info(li).datatype)==(101,0):
            continue
        values=sorted(str(s) for s in cell.shapes(li).each())
        if values:
            rows.append((str(layout.get_info(li)),values))
    return hashlib.sha256(json.dumps(rows).encode()).hexdigest()

layout=pya.Layout()
layout.read(payload['input'])
top=layout.top_cell()
active=layout.cell('DUAL_RAIL_ACTIVE_LN')
muxes=[c for c in layout.each_cell() if c.name.startswith('MUX_TE01_A70_TAIL030_STRAIGHT_PORTS')]
if top is None or active is None or len(muxes)!=1 or layout.dbu!=.001 or len(layout.top_cells())!=1:
    raise RuntimeError('v9结构、格点或顶层不匹配')
mux=muxes[0]
instances=[i for i in top.each_inst() if i.cell==mux]
if len(instances)!=4:
    raise RuntimeError('MUX实例数不是四个')
other_cells={c.name:signature(layout,c) for c in layout.each_cell() if c not in (top,active,mux)}
before_core=pya.Region(top.begin_shapes_rec(layout.layer(20,0))).merged()
before_metal=pya.Region(top.begin_shapes_rec(layout.layer(42,0))).merged()
before_instances=sorted(str(i.dcplx_trans) for i in top.each_inst())

# 有源双轨20/1由7.2um扩大为8.1um，中心线和长度不变。
active_li=layout.layer(20,1)
active_shapes=list(active.shapes(active_li).each())
expected_active=[pya.Box(2000000,26400,17000000,33600),pya.Box(2000000,-33600,17000000,-26400)]
if sorted(str(s.box) for s in active_shapes)!=sorted(str(b) for b in expected_active):
    raise RuntimeError('有源20/1不是两条预期7.2um矩形')
active.shapes(active_li).clear()
for center in (30000,-30000):
    active.shapes(active_li).insert(pya.Box(2000000,center-4050,17000000,center+4050))

# MUX窗口整体保持平行直线；以主脊中心-1.95um为共同中心。
for layer_key,old_box,new_box in [
    ((20,1),pya.Box(-200000,-3600,950000,3600),pya.Box(-200000,-6000,950000,2100)),
    ((10,2),pya.Box(-200000,-8600,950000,8600),pya.Box(-200000,-10550,950000,6650))]:
    li=layout.layer(*layer_key)
    shapes=list(mux.shapes(li).each())
    if len(shapes)!=1 or not shapes[0].is_box() or shapes[0].box!=old_box:
        raise RuntimeError('MUX原窗口不匹配：'+str(layer_key))
    mux.shapes(li).clear()
    mux.shapes(li).insert(new_box)

mux.name='MUX_TE01_A70_TAIL030_PARALLEL_WINDOWS_UNVERIFIED'
active.name='DUAL_RAIL_ACTIVE_LN_PARALLEL_WINDOWS_UNVERIFIED'
top.name='EO4P_10G_15MM_PARALLEL_WINDOWS_DRAFT'
mux.shapes(layout.layer(101,0)).insert(pya.Text('PARALLEL_MASK_WINDOWS_20_1_W8P1_10_2_W17P2',pya.Trans(-190000,-6800)))

core=pya.Region(top.begin_shapes_rec(layout.layer(20,0))).merged()
clad=pya.Region(top.begin_shapes_rec(layout.layer(20,1))).merged()
trench=pya.Region(top.begin_shapes_rec(layout.layer(10,2))).merged()
metal=pya.Region(top.begin_shapes_rec(layout.layer(42,0))).merged()
mux_core=pya.Region(mux.shapes(layout.layer(20,0))).merged()
mux_clad=pya.Region(mux.shapes(layout.layer(20,1))).merged()
mux_trench=pya.Region(mux.shapes(layout.layer(10,2))).merged()
checks={
 'four_mux_instances':len(instances)==4,
 'all_other_cells_unchanged':all(signature(layout,layout.cell(n))==h for n,h in other_cells.items()),
 'all_20_0_core_unchanged':(core^before_core).is_empty(),
 'all_M1_metal_unchanged':(metal^before_metal).is_empty(),
 'all_instance_transforms_unchanged':before_instances==sorted(str(i.dcplx_trans) for i in top.each_inst()),
 'mux_core_inside_20_1':(mux_core-mux_clad).is_empty(),
 'mux_core_inside_10_2':(mux_core-mux_trench).is_empty(),
 'no_metal_core_overlap':(metal&core).is_empty(),
}
joint_rows=[]
for index,inst in enumerate(instances,1):
    row={'index':index,'transform':str(inst.dcplx_trans)}
    for name,layer,width in [('20/1',layout.layer(20,1),8.1),('10/2',layout.layer(10,2),17.2)]:
        test=pya.Region(pya.DPolygon(pya.DBox(-200.001,-1.95-width/2,-199.999,-1.95+width/2)).transformed(inst.dcplx_trans).to_itype(layout.dbu))
        region=pya.Region(top.begin_shapes_rec(layer)).merged()
        ok=(test-region).is_empty()
        row[name]={'full_width_connected':ok,'width_um':width}
        checks['active_join_'+str(index)+'_'+name]=ok
    joint_rows.append(row)
if not all(checks.values()):
    raise RuntimeError('平行窗口检查失败：'+json.dumps(checks))
etch1=clad-core
markers={'core_width_030':core.width_check(300).size(),'core_space_030':core.space_check(300).size(),
         'etch1_width_030':etch1.width_check(300).size(),'etch1_space_030':etch1.space_check(300).size(),
         'sin_trench_width_020':trench.width_check(200).size(),'sin_trench_space_020':trench.space_check(200).size(),
         'metal_width_2':metal.width_check(2000).size(),'metal_space_3':metal.space_check(3000).size()}
layout.write(payload['output'])
{'status':'直LN1脊+平行公共刻蚀窗口规划；非光学或流片签核','top_cell':top.name,
 'mux_cell':mux.name,'active_cell':active.name,'instances_updated':4,'joint_checks':joint_rows,
 'checks':checks,'manual_rule_markers':markers,
 'window_parameters_um':{'20/1_width':8.1,'10/2_width':17.2,'local_center':-1.95,
                         'mux_core_y_min':mux_core.bbox().bottom*.001,'mux_core_y_max':mux_core.bbox().top*.001,
                         'minimum_nominal_transverse_margin':.3},
 'process_geometry_changed':True,'optical_validation_pass':False}
'''


def main() -> None:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',type=Path,default=BASE)
    p.add_argument('--output-dir',type=Path,required=True)
    a=p.parse_args()
    a.output_dir.mkdir(parents=True,exist_ok=False)
    output=(a.output_dir/'四程10GHz_15mm_直脊平行公共窗口_待验证工作稿.gds').resolve()
    source_hash=hashlib.sha256(a.input.read_bytes()).hexdigest()
    client=KLinkClient()
    client.connect()
    try:
        response=client.call('exec.python',{'code':'payload = '+repr({'input':str(a.input.resolve()),'output':str(output)})+'\n'+CODE})
        if response.get('exception'):
            raise RuntimeError(str(response['exception']))
        report=response['return_value']
        report['file_info']=client.call('layout.file_info',{'path':str(output),'detail':'counts'})
        client.call('layout.show_file',{'path':str(output),'mode':'new'})
        client.call('view.show_cell',{'cell':report['top_cell'],'zoom_fit':True})
        for layer,color in [('20/0','#C62858'),('20/1','#339A90'),('10/2','#65A7CF'),('42/0','#D5A12C')]:
            client.call('layer.set_style',{'layer':layer,'fill_color':color,'frame_color':color,
                        'dither_pattern':0 if layer in ('20/0','42/0') else 1,'line_width':1})
        client.call('layer.set_visible',{'layers':['101/0','100/0','100/30','56/30'],'visible':False})
        client.call('view.hier_levels',{'min':0,'max':10})
        for name,bbox in [('左侧四个平行窗口接头.png',[1725,-45,2075,45]),
                          ('右侧四个平行窗口接头.png',[16925,-45,17275,45])]:
            client.call('view.screenshot',{'mode':'path','path':str((a.output_dir/name).resolve()),
                                          'bbox_um':bbox,'width_px':1900,'height_px':650})
        client.call('view.zoom_box',{'bbox_um':[1725,20,2075,42]})
    finally:
        client.close()
    if hashlib.sha256(a.input.read_bytes()).hexdigest()!=source_hash:
        raise RuntimeError('v9源文件发生变化')
    report.update({'working_gds':str(output),'source_gds':str(a.input.resolve()),
                   'source_sha256':source_hash,'source_unchanged':True,'foundry_signoff':False,
                   'limitations':['8.1um的20/1窗口改变了原EME/HFSS截面，必须重新验证光学模式和射频速度',
                                  '外部主/辅助双支路的包络分离端仍需在最终截面下复核',
                                  'MUX本体辅助脊仍有0.3um法向宽度/格点标记',
                                  '21层、圆弯、10GHz时延、RF焊盘和终端仍按用户要求暂缓']})
    (a.output_dir/'平行公共窗口检查.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k!='file_info'},ensure_ascii=False,indent=2))


if __name__=='__main__':
    main()
