"""仅修正MUX辅助尾端的0.15um违规，另存GDS，不改本体或宣称低反射通过。"""
import argparse
import hashlib
import json
from pathlib import Path

from klink import KLinkClient

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT/'results/layout/_历史归档/迭代版本_20260912/黑盒交叉_拉锥工作稿_v5/四程10GHz_15mm_PDK黑盒交叉_拉锥接续_待验证工作稿.gds'

CODE = r'''
import hashlib
import json
def signature(layout,cell):
    shapes = [(str(layout.get_info(li)), sorted(str(s) for s in cell.shapes(li).each()))
              for li in layout.layer_indices() if not cell.shapes(li).is_empty()]
    return hashlib.sha256(json.dumps([sorted(shapes),sorted((i.cell.name,str(i.dcplx_trans)) for i in cell.each_inst())]).encode()).hexdigest()
layout=pya.Layout()
layout.read(payload['input'])
top=layout.top_cell()
mux=layout.cell('MUX_TE01_A70_EME_NOMINAL')
if mux is None or layout.dbu!=.001 or len(layout.top_cells())!=1:
    raise RuntimeError('基线结构与预期不符')
li=layout.layer(20,0)
old_core=pya.Region(top.begin_shapes_rec(li)).merged()
old_mux=pya.Region(mux.shapes(li)).merged()
other_cells={c.name:signature(layout,c) for c in layout.each_cell() if c!=top and c!=mux}
other_layers={str(layout.get_info(i)):pya.Region(top.begin_shapes_rec(i)).merged()
              for i in layout.layer_indices() if layout.get_info(i).layer not in (20,101)}
old_tail=pya.Region(pya.Polygon([pya.Point(-100000,1325),pya.Point(0,1250),
                               pya.Point(0,1550),pya.Point(-100000,1475)]))
candidates=[s for s in mux.shapes(li).each() if s.bbox().left==-100000 and s.bbox().right==0]
if len(candidates)!=1 or not (pya.Region(candidates[0].polygon)^old_tail).is_empty():
    raise RuntimeError('无法唯一匹配原0.15->0.30um尾端，拒绝盲改')
instances=[i for i in top.each_inst() if i.cell==mux]
if len(instances)!=4:
    raise RuntimeError('MUX实例不是四个')
old_markers={'core_width':old_core.width_check(300).size()}
old_etch=pya.Region(top.begin_shapes_rec(layout.layer(20,1))).merged()-old_core
old_markers['etch1_space']=old_etch.space_check(300).size()
new_tail=pya.Region(pya.Box(-100000,1250,0,1550))
candidates[0].delete()
mux.shapes(li).insert(new_tail)
text_shapes=mux.shapes(layout.layer(101,0))
for s in list(text_shapes.each()):
    if s.is_text() and s.text.string=='AUX_TAIL_0P15_PENDING_EME':
        s.delete()
text_shapes.insert(pya.Text('AUX_TAIL_0P30_REFLECTION_UNVERIFIED',pya.Trans(-100000,1400)))
mux.name='MUX_TE01_A70_TAIL030_UNVERIFIED'
top.name='EO4P_10G_15MM_TAIL030_DRAFT'
core=pya.Region(top.begin_shapes_rec(li)).merged()
allowed=pya.Region()
locations=[]
for inst in instances:
    allowed+=new_tail.transformed(inst.cplx_trans)
    point=inst.dcplx_trans*pya.DPoint(-100,1.4)
    joint=inst.dcplx_trans*pya.DPoint(0,1.4)
    locations.append({'tip_um':[point.x,point.y],'mux_joint_um':[joint.x,joint.y]})
checks={
 'other_cells_unchanged':all(signature(layout,layout.cell(n))==h for n,h in other_cells.items()),
 'core_unchanged_outside_four_tails':((core^old_core)-allowed).is_empty(),
 'mux_core_unchanged_outside_tail':((pya.Region(mux.shapes(li)).merged()^old_mux)-new_tail).is_empty(),
 'tail_width_030_pass':new_tail.width_check(300).is_empty(),
 'tail_to_body_connected':(pya.Region(pya.Box(-1,1250,1,1550))-pya.Region(mux.shapes(li)).merged()).is_empty(),
 'other_process_layers_unchanged':all((pya.Region(top.begin_shapes_rec(i)).merged()^other_layers[str(layout.get_info(i))]).is_empty()
                                    for i in layout.layer_indices() if str(layout.get_info(i)) in other_layers),
}
# 20/1包络不得随本次尾端线宽修改。
source=pya.Layout()
source.read(payload['input'])
checks['clad20_1_unchanged']=(pya.Region(source.top_cell().begin_shapes_rec(source.layer(20,1)))^
                            pya.Region(top.begin_shapes_rec(layout.layer(20,1)))).is_empty()
if not all(checks.values()):
    raise RuntimeError('尾端局部检查失败：'+json.dumps(checks))
etch=pya.Region(top.begin_shapes_rec(layout.layer(20,1))).merged()-core
new_markers={'core_width':core.width_check(300).size(),'etch1_space':etch.space_check(300).size()}
remaining=[]
for pair in core.width_check(300).each():
    p=pair.bbox().center()
    for inst in instances:
        q=inst.dcplx_trans.inverted()*pya.DPoint(p.x*.001,p.y*.001)
        if -200<=q.x<=950 and -8.6<=q.y<=8.6:
            remaining.append([q.x,q.y])
layout.write(payload['output'])
{'top_cell':top.name,'tail_length_um':100,'old_tip_width_um':.15,'new_tail_width_um':.30,
 'locations':locations,'checks':checks,'markers_before':old_markers,'markers_after':new_markers,
 'remaining_marker_centers_mux_local_um':remaining,'optical_validation_pass':False,
 'note':'仅修正末端100um违规尺寸；平头终止尚未验证反射。MUX本体不改变。'}
'''

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',type=Path,default=BASE)
    p.add_argument('--output-dir',type=Path,required=True)
    a=p.parse_args()
    before=hashlib.sha256(a.input.read_bytes()).hexdigest()
    a.output_dir.mkdir(parents=True,exist_ok=False)
    output=(a.output_dir/'四程10GHz_15mm_尾端线宽修正_待验证工作稿.gds').resolve()
    c=KLinkClient()
    c.connect()
    try:
        r=c.call('exec.python',{'code':'payload = '+repr({'input':str(a.input.resolve()),'output':str(output)})+'\n'+CODE})
        if r.get('exception'):
            raise RuntimeError(str(r['exception']))
        report=r['return_value']
        report['file_info']=c.call('layout.file_info',{'path':str(output),'detail':'counts'})
        c.call('layout.show_file',{'path':str(output),'mode':'new'})
        c.call('view.show_cell',{'cell':report['top_cell'],'zoom_fit':True})
        for layer,color in [('20/0','#C62858'),('20/1','#339A90'),('10/2','#65A7CF'),('42/0','#D5A12C')]:
            c.call('layer.set_style',{'layer':layer,'fill_color':color,'frame_color':color,
                                      'dither_pattern':0 if layer in ('20/0','42/0') else 1,'line_width':1})
        c.call('layer.set_visible',{'layers':['101/0','100/0','100/30','56/30'],'visible':False})
        c.call('view.hier_levels',{'min':0,'max':10})
        bbox=[1775,27,1925,39]
        c.call('view.screenshot',{'mode':'path','path':str((a.output_dir/'左上MUX尾端_实际GDS.png').resolve()),
                                 'bbox_um':bbox,'width_px':1800,'height_px':340})
        c.call('view.zoom_box',{'bbox_um':bbox})
    finally:
        c.close()
    if hashlib.sha256(a.input.read_bytes()).hexdigest()!=before:
        raise RuntimeError('原文件哈希发生变化')
    report.update({'working_gds':str(output),'source_gds':str(a.input.resolve()),'source_sha256':before,
                   'original_unchanged':True,'remaining':['MUX邻近本体最小线宽标记尚未关闭',
                   '平头末端的反射和与主路回耦未验证','未修改21层、欧拉弯或射频接口；继承v5其余限制']})
    (a.output_dir/'尾端修正检查.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k!='file_info'},ensure_ascii=False,indent=2))

if __name__=='__main__':
    main()
