"""在15 mm电极输出端加入估算50 ohm的HTR终端，删除右侧GSG输出探针。"""
import argparse, hashlib, json
from pathlib import Path
from klink import KLinkClient

ROOT=Path(__file__).resolve().parents[2]
p=argparse.ArgumentParser()
p.add_argument('--source',type=Path,required=True)
p.add_argument('--output-dir',type=Path,required=True)
p.add_argument('--port',type=int,default=8765)
a=p.parse_args(); a.source=a.source.resolve(); a.output_dir=a.output_dir.resolve()
a.output_dir.mkdir(parents=True,exist_ok=False)
output=(a.output_dir/'YSJ与四程15mm_M04_MUX_750um_50R终端_功能候选.gds').resolve()
report=(a.output_dir/'50R终端版图检查.json').resolve()
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
if len(launches)!=2: raise RuntimeError('未找到两个GSG150输入输出单元')
out_launch=[i for i in launches if str(i.dcplx_trans)=='m90 *1 6015,1620']
in_launch=[i for i in launches if str(i.dcplx_trans)=='r0 *1 -10255,1620']
if len(out_launch)!=1 or len(in_launch)!=1: raise RuntimeError('GSG150输入输出变换不符合冻结版')
out_launch[0].delete()
term=layout.create_cell('RF_TERMINATION_50R_RS20_EST_DRAFT')
htr=layout.layer(41,0); m1=layout.layer(42,0); ann=layout.layer(101,0)
term.shapes(htr).insert(pya.Box(0,-6250,31250,6250))
term.shapes(m1).insert(pya.Box(-10000,-10000,10000,10000))
term.shapes(m1).insert(pya.Box(25000,-10000,50000,10000))
term.shapes(m1).insert(pya.Box(-1000,38500,41000,138500))
term.shapes(m1).insert(pya.Box(-1000,-138500,41000,-38500))
term.shapes(m1).insert(pya.Box(40000,-138500,60000,138500))
term.shapes(ann).insert(pya.Text('RF_TERMINATION_50R_EST_RS20',pya.Trans(0,-170000)))
# EO单元中的电极实例位于 x=-9620..5380 µm，因此终端放在EO坐标 x=5380，
# 而不是把电极单元内部的15000 µm误当成EO顶层坐标。
eo.insert(pya.DCellInstArray(term.cell_index(),pya.DCplxTrans(1,0,False,5380,1620)))
top.name='YSJ_EO4P_15MM_M04_50R_TERM_HORIZONTAL_DRAFT'
term_htr=pya.Region(term.shapes(htr)).merged(); term_m1=pya.Region(term.shapes(m1)).merged()
checks={
 'one_top':len(layout.top_cells())==1,
 'one_input_gsg150':len([i for i in eo.each_inst() if i.cell.name=='GSG150_LAUNCH_50UM_DRAFT'])==1,
 'right_gsg_removed':not any(i.cell.name=='GSG150_LAUNCH_50UM_DRAFT' and str(i.dcplx_trans)=='m90 *1 6015,1620' for i in eo.each_inst()),
 'one_termination':len([i for i in eo.each_inst() if i.cell==term])==1,
 'htr_2p5_squares':term_htr.bbox()==pya.Box(0,-6250,31250,6250),
 'signal_ground_m1_not_short':(pya.Region(pya.Box(-10000,-10000,10000,10000)) & pya.Region(pya.Box(25000,-10000,50000,10000))).is_empty(),
 'ground_bridge_hits_upper':not (pya.Region(pya.Box(40000,-138500,60000,138500)) & pya.Region(pya.Box(-1000,38500,41000,138500))).is_empty(),
 'ground_bridge_hits_lower':not (pya.Region(pya.Box(40000,-138500,60000,138500)) & pya.Region(pya.Box(-1000,-138500,41000,-38500))).is_empty(),
 'metal_min_width_2um':term_m1.width_check(2000).is_empty(),
 'metal_min_space_3um':term_m1.space_check(3000).is_empty(),
 'optical_layers_unchanged':all((opt_before[k]^region(top,k)).is_empty() for k in opt_keys),
'no_new_ln2':term.shapes(layout.layer(21,0)).is_empty() and term.shapes(layout.layer(21,1)).is_empty() and term.shapes(layout.layer(21,2)).is_empty(),
}
einst=[i for i in eo.each_inst() if i.cell.name=='T_GSG_15MM_RECTANGULAR_BASELINE']
tinst=[i for i in eo.each_inst() if i.cell==term]
if len(einst)!=1 or len(tinst)!=1: raise RuntimeError('电极或终端实例数量不符')
electrode_m1=pya.Region(einst[0].cell.shapes(m1)).transformed_icplx(einst[0].dcplx_trans)
term_m1_eo=term_m1.transformed_icplx(tinst[0].dcplx_trans)
term_htr_eo=term_htr.transformed_icplx(tinst[0].dcplx_trans)
sig_pad_eo=pya.Region(pya.Box(-10000,-10000,10000,10000)).transformed_icplx(tinst[0].dcplx_trans)
gup_eo=pya.Region(pya.Box(-1000,38500,41000,138500)).transformed_icplx(tinst[0].dcplx_trans)
gdn_eo=pya.Region(pya.Box(-1000,-138500,41000,-38500)).transformed_icplx(tinst[0].dcplx_trans)
checks['signal_pad_hits_electrode']=not (sig_pad_eo & electrode_m1).is_empty()
checks['upper_ground_hits_electrode']=not (gup_eo & electrode_m1).is_empty()
checks['lower_ground_hits_electrode']=not (gdn_eo & electrode_m1).is_empty()
checks['htr_hits_signal_pad']=not (term_htr_eo & sig_pad_eo).is_empty()
checks['htr_hits_ground_pad']=not (term_htr_eo & term_m1_eo).is_empty()
tb=top.bbox(); checks['within_21p8mm_block']=(tb.left*layout.dbu>=-10900 and tb.right*layout.dbu<=10900.002 and tb.bottom*layout.dbu>=-1900 and tb.top*layout.dbu<=1900)
if not all(checks.values()): raise RuntimeError('50R终端几何检查失败:'+json.dumps(checks))
layout.write(payload['output'])
result={'status':'50R_HTR_termination_geometry_candidate','output_gds':payload['output'],'source_gds':payload['source'],'source_sha256':payload['source_sha256'],'top_cell':top.name,'termination_cell':term.name,'removed_output_gsg150':True,'retained_input_gsg150':True,'sheet_resistance_ohm_per_square':20.0,'target_resistance_ohm':50.0,'resistor_width_um':12.5,'effective_resistor_length_um':31.25,'square_count':2.5,'changed_layers':['41/0','42/0','101/0'],'checks':checks,'optical_validation_pass':False,'rf_50ohm_validation_pass':False,'note':'按HTR方阻20 ohm/square几何估算；未包含接触电阻、寄生和S参数。'}
open(payload['report'],'w',encoding='utf-8').write(json.dumps(result,ensure_ascii=False,indent=2)); result
'''
payload={'source':str(a.source),'output':str(output),'report':str(report),'source_sha256':before}
c=KLinkClient(host='127.0.0.1',port=a.port); c.connect()
try:
    r=c.call('exec.python',{'code':'payload='+repr(payload)+'\n'+code})
    if r.get('exception'): raise RuntimeError(str(r['exception']))
    result=r['return_value']; result['file_info']=c.call('layout.file_info',{'path':str(output),'detail':'counts'})
    c.call('layout.show_file',{'path':str(output),'mode':'new'})
    for layer,color in [('41/0','#B87333'),('42/0','#D5A12C'),('20/0','#C62858'),('20/1','#55A995')]:
        c.call('layer.set_style',{'layer':layer,'fill_color':color,'frame_color':color,'dither_pattern':0 if layer in ('41/0','42/0','20/0') else 1,'line_width':1})
    c.call('layer.set_visible',{'layers':['10/0','10/1','10/2','20/2','21/2','56/30','70/30','100/0','100/30','101/0'],'visible':False})
    c.call('view.screenshot',{'mode':'path','path':str((a.output_dir/'50R终端局部.png').resolve()),'bbox_um':[5150,1450,5500,1800],'width_px':2200,'height_px':900})
finally: c.close()
if hashlib.sha256(a.source.read_bytes()).hexdigest()!=before: raise RuntimeError('源GDS发生变化')
print(json.dumps(result,ensure_ascii=False,indent=2))
