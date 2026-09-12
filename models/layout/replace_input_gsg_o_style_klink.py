"""用PDK o型直线CPW风格替换输入GSG150引出，保留输出50R终端。"""
import argparse, hashlib, json
from pathlib import Path
from klink import KLinkClient

p=argparse.ArgumentParser()
p.add_argument('--source',type=Path,required=True)
p.add_argument('--output-dir',type=Path,required=True)
p.add_argument('--port',type=int,default=8765)
a=p.parse_args(); a.source=a.source.resolve(); a.output_dir=a.output_dir.resolve()
a.output_dir.mkdir(parents=True,exist_ok=False)
output=(a.output_dir/'YSJ与四程15mm_M04_MUX_750um_50R终端_O型输入_功能候选.gds').resolve()
report=(a.output_dir/'O型输入版图检查.json').resolve()
before=hashlib.sha256(a.source.read_bytes()).hexdigest()
code=r'''
import json, pya
layout=pya.Layout(); layout.read(payload['source'])
top=layout.top_cell(); eo=layout.cell('EO4P_15MM_GSG150_HORIZONTAL_DRAFT')
if top is None or eo is None or len(layout.top_cells())!=1 or abs(layout.dbu-.001)>1e-12:
    raise RuntimeError('输入顶层或格点不符合要求')
def region(cell,key): return pya.Region(cell.begin_shapes_rec(layout.layer(*key))).merged()
opt_keys=[(10,0),(10,1),(10,2),(20,0),(20,1),(20,2),(21,0),(21,1),(21,2)]
opt_before={k:region(top,k) for k in opt_keys}
launches=[i for i in eo.each_inst() if i.cell.name=='GSG150_LAUNCH_50UM_DRAFT']
if len(launches)!=1 or str(launches[0].dcplx_trans)!='r0 *1 -10255,1620':
    raise RuntimeError('输入GSG150实例不符合v5h基线')
old_launch_cell=launches[0].cell
launches[0].delete()
if all(i.cell!=old_launch_cell for c in layout.each_cell() for i in c.each_inst()):
    layout.delete_cell(old_launch_cell.cell_index())
cell=layout.create_cell('GSG150_LAUNCH_O_STYLE_50UM_DRAFT')
m1=layout.layer(42,0); ann=layout.layer(101,0)
# O型直线式输入：50um落针段 -> 100um直线渐变 -> 43um电极主干；总长321um。
x=[0.,120.,220.,321.]
def poly(lower,upper):
    return pya.DPolygon([pya.DPoint(x[k],lower[k]) for k in range(4)]+[pya.DPoint(x[k],upper[k]) for k in range(3,-1,-1)]).to_itype(layout.dbu)
sw=[25.,25.,21.5,21.5]; si=[38.5,38.5,38.5,38.5]; so=[190.,190.,138.5,138.5]
cell.shapes(m1).insert(poly([-v for v in sw],sw))
cell.shapes(m1).insert(poly(si,so))
cell.shapes(m1).insert(poly([-v for v in so],[-v for v in si]))
for net,y in [('Gupper',150.),('S',0.),('Glower',-150.)]:
    cell.shapes(ann).insert(pya.DBox(25-25,y-25,25+25,y+25).to_itype(layout.dbu))
    cell.shapes(ann).insert(pya.Text(net+'_CONTACT',pya.Trans(25,y)))
cell.shapes(ann).insert(pya.Text('O_STYLE_GSG150_INPUT',pya.Trans(0,-220)))
new_t=pya.DCplxTrans(1,0,False,-9940,1620)
new_inst=pya.DCellInstArray(cell.cell_index(),new_t)
eo.insert(new_inst)
top.name='YSJ_EO4P_15MM_M04_50R_TERM_O_INPUT_HORIZONTAL_DRAFT'
# 连接与规则检查
checks={
 'one_top':len(layout.top_cells())==1,
 'one_o_input':len([i for i in eo.each_inst() if i.cell==cell])==1,
 'old_input_removed':not any(i.cell.name=='GSG150_LAUNCH_50UM_DRAFT' for i in eo.each_inst()),
 'output_gsg_already_removed':not any(i.cell.name=='GSG150_LAUNCH_50UM_DRAFT' and str(i.dcplx_trans)=='m90 *1 6015,1620' for i in eo.each_inst()),
 'straight_input_length_321um':cell.bbox().right*layout.dbu==321.,
 'm1_min_width_2um':pya.Region(cell.shapes(m1)).width_check(2000).is_empty(),
 'm1_min_space_3um':pya.Region(cell.shapes(m1)).space_check(3000).is_empty(),
 'optical_layers_unchanged':all((opt_before[k]^region(top,k)).is_empty() for k in opt_keys),
 'no_new_ln2':cell.shapes(layout.layer(21,0)).is_empty() and cell.shapes(layout.layer(21,1)).is_empty() and cell.shapes(layout.layer(21,2)).is_empty()
}
einst=[i for i in eo.each_inst() if i.cell.name=='T_GSG_15MM_RECTANGULAR_BASELINE']
if len(einst)!=1: raise RuntimeError('电极实例缺失')
electrode_m1=pya.Region(einst[0].cell.shapes(m1)).transformed_icplx(einst[0].dcplx_trans)
new_m1=pya.Region(cell.shapes(m1)).transformed_icplx(new_t)
# 信号及两地都应与15mm电极端部产生面积重叠，不仅是边界相切。
sig=pya.Region(pya.Box(-10000,-10000,10000,10000)).transformed_icplx(new_t)
gu=pya.Region(pya.Box(0,38500,41000,138500)).transformed_icplx(new_t)
gl=pya.Region(pya.Box(0,-138500,41000,-38500)).transformed_icplx(new_t)
checks['signal_hits_15mm_electrode']=not (sig & electrode_m1).is_empty()
checks['upper_ground_hits_15mm_electrode']=not (gu & electrode_m1).is_empty()
checks['lower_ground_hits_15mm_electrode']=not (gl & electrode_m1).is_empty()
tb=top.bbox(); checks['within_21p8mm_block']=(tb.left*layout.dbu>=-10900 and tb.right*layout.dbu<=10900.002 and tb.bottom*layout.dbu>=-1900 and tb.top*layout.dbu<=1900)
if not all(checks.values()): raise RuntimeError('O型输入版图检查失败:'+json.dumps(checks))
layout.write(payload['output'])
result={'status':'O_style_GSG150_input_geometry_candidate','output_gds':payload['output'],'source_gds':payload['source'],'source_sha256':payload['source_sha256'],'top_cell':top.name,'input_cell':cell.name,'removed_right_gsg150':True,'retained_input_pitch_um':150.0,'input_geometry_length_um':321.0,'input_pad_width_um':50.0,'input_trunk_width_um':43.0,'changed_layers':['42/0','101/0'],'checks':checks,'optical_validation_pass':False,'rf_validation_pass':False,'note':'参考PDK o型直线CPW风格；未复制或修改黑盒。GSG150到43um主干采用两段直线过渡，RF阻抗仍需求解。'}
open(payload['report'],'w',encoding='utf-8').write(json.dumps(result,ensure_ascii=False,indent=2)); result
'''
payload={'source':str(a.source),'output':str(output),'report':str(report),'source_sha256':before}
c=KLinkClient(host='127.0.0.1',port=a.port); c.connect()
try:
    r=c.call('exec.python',{'code':'payload='+repr(payload)+'\n'+code})
    if r.get('exception'): raise RuntimeError(str(r['exception']))
    result=r['return_value']; result['file_info']=c.call('layout.file_info',{'path':str(output),'detail':'counts'})
    c.call('layout.show_file',{'path':str(output),'mode':'new'})
    for layer,color in [('42/0','#D5A12C'),('41/0','#E64A19'),('20/0','#C62858'),('20/1','#55A995')]:
        c.call('layer.set_style',{'layer':layer,'fill_color':color,'frame_color':color,'dither_pattern':0 if layer in ('42/0','41/0','20/0') else 1,'line_width':1})
    c.call('layer.set_visible',{'layers':['10/0','10/1','10/2','20/2','21/2','56/30','70/30','100/0','100/30','101/0'],'visible':False})
    c.call('view.screenshot',{'mode':'path','path':str((a.output_dir/'O型RF输入局部.png').resolve()),'bbox_um':[-10000,1400,-9500,1840],'width_px':2200,'height_px':1100})
finally: c.close()
if hashlib.sha256(a.source.read_bytes()).hexdigest()!=before: raise RuntimeError('源GDS发生变化')
print(json.dumps(result,ensure_ascii=False,indent=2))
