"""仅修正自定义HTR终端的有效方块数；不改光路、不启动求解、不处理地桥跨光芯。"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
from klink import KLinkClient

PAYLOAD_CODE = r"""
import hashlib, json, pya
l=pya.Layout(); l.read(payload['source'])
if len(l.top_cells())!=1 or l.dbu!=.001:
    raise RuntimeError('源文件顶层或格点不符')
top=l.top_cell(); eo=l.cell('EO4P_15MM_GSG150_HORIZONTAL_DRAFT')
term=l.cell('RF_TERMINATION_50R_RS20_EST_DRAFT')
if eo is None or term is None: raise RuntimeError('未找到当前电极与终端单元')
def reg(layout,cell,key):
    li=layout.find_layer(*key)
    return pya.Region() if li is None else pya.Region(cell.begin_shapes_rec(li)).merged()
def signature(layout,c):
    return hashlib.sha256(json.dumps({
        'shapes':sorted((str(layout.get_info(li)),sorted(str(s) for s in c.shapes(li).each()))
                        for li in layout.layer_indices() if not c.shapes(li).is_empty()),
        'instances':sorted((i.cell.name,str(i.cplx_trans)) for i in c.each_inst())
    },sort_keys=True).encode()).hexdigest()
def bb(b):
    return [v*l.dbu for v in [b.left,b.bottom,b.right,b.top]]
unchanged={c.name:signature(l,c) for c in l.each_cell() if c not in (term,top)}
old_source_bbox=top.bbox()
htr_li=l.layer(41,0); metal_li=l.layer(42,0)
old_htr=reg(l,term,(41,0)); old_metal=reg(l,term,(42,0))
expected_old_htr=pya.Region(pya.Box(0,-6250,31250,6250))
expected_old_metal=pya.Region()
for b in [(-10000,-10000,10000,10000),(25000,-10000,50000,10000),
          (-1000,38500,41000,138500),(-1000,-138500,41000,-38500),
          (40000,-138500,60000,138500)]:
    expected_old_metal.insert(pya.Box(*b))
if not (old_htr^expected_old_htr).is_empty() or not (old_metal^expected_old_metal).is_empty():
    raise RuntimeError('终端不是已核查的旧版几何，拒绝盲改')
ti=[i for i in eo.each_inst() if i.cell==term]
if len(ti)!=1 or str(ti[0].dcplx_trans)!='r0 *1 5380,1620':
    raise RuntimeError('终端摆位不符')
core=reg(l,eo,(20,0))
old_overlap=(old_metal.transformed(ti[0].cplx_trans)&core).merged()
rs=payload['sheet_resistance']; target=payload['target']; width=payload['width']
contact_overlap=10.
leff=target/rs*width
signal_edge=10.
ground_edge=signal_edge+leff
draw_end=ground_edge+contact_overlap
bridge_start=ground_edge+15.
bridge_end=bridge_start+20.
def B(x0,y0,x1,y1):
    return pya.DBox(x0,y0,x1,y1).to_itype(l.dbu)
term.shapes(htr_li).clear(); term.shapes(metal_li).clear()
term.shapes(htr_li).insert(B(0,-width/2,draw_end,width/2))
new_boxes=[
    (-10.,-10.,signal_edge,10.),
    (ground_edge,-10.,ground_edge+25.,10.),
    (-1.,38.5,bridge_start+1.,138.5),
    (-1.,-138.5,bridge_start+1.,-38.5),
    (bridge_start,-138.5,bridge_end,138.5)]
for b in new_boxes: term.shapes(metal_li).insert(B(*b))
top.name='YSJ_EO4P_M04_OINPUT_HTR50_LEFF_V7'

def audit(layout):
    t=layout.cell(term.name); e=layout.cell(eo.name); tp=layout.top_cell()
    inst=[i for i in e.each_inst() if i.cell==t][0]
    h=reg(layout,t,(41,0)); m=reg(layout,t,(42,0))
    mm=reg(layout,e,(42,0))
    expected=pya.Region(B(signal_edge,-width/2,ground_edge,width/2))
    exposed=(h-m).merged()
    all_exposed=(h.transformed(inst.cplx_trans)-mm).merged()
    pieces=list(mm.each())
    def component(x,y):
        anchor=pya.Region(B(x-.1,y-.1,x+.1,y+.1))
        return [i for i,p in enumerate(pieces) if (anchor-pya.Region(p)).is_empty()]
    ids={'signal':component(-9915,1620),'upper_ground':component(-9915,1770),
         'lower_ground':component(-9915,1470)}
    distinct=all(len(v)==1 for v in ids.values()) and ids['upper_ground']==ids['lower_ground'] and ids['signal']!=ids['upper_ground']
    contacts={}
    if distinct:
        hworld=h.transformed(inst.cplx_trans)
        contacts={'signal_contact_area_um2':(hworld&pya.Region(pieces[ids['signal'][0]])).area()*layout.dbu**2,
                  'ground_contact_area_um2':(hworld&pya.Region(pieces[ids['upper_ground'][0]])).area()*layout.dbu**2}
    new_overlap=(m.transformed(inst.cplx_trans)&reg(layout,e,(20,0))).merged()
    actual_l=exposed.bbox().width()*layout.dbu
    actual_w=exposed.bbox().height()*layout.dbu
    checks={
        'one_top':len(layout.top_cells())==1,
        'all_other_cells_unchanged':all(signature(layout,layout.cell(name))==digest for name,digest in unchanged.items()),
        'exposed_region_is_exact_rectangle':(exposed^expected).is_empty(),
        'no_other_metal_covers_resistor_gap':(all_exposed^expected.transformed(inst.cplx_trans)).is_empty(),
        'resistance_50R_under_ideal_contact_model':abs(rs*actual_l/actual_w-target)<1e-9,
        'signal_and_common_ground_are_distinct_metal_components':distinct,
        'both_resistor_contacts_overlap_metal':len(contacts)==2 and all(v>0 for v in contacts.values()),
        'local_M1_width_2um':m.width_check(2000).is_empty(),
        'local_M1_space_3um':m.space_check(3000).is_empty(),
        'local_HTR_width_2um':h.width_check(2000).is_empty(),
        'new_terminal_inside_physical_block':((m+h).transformed(inst.cplx_trans)-pya.Region(B(-10900,-1900,10900,1900))).is_empty(),
        'full_bbox_unchanged':tp.bbox()==old_source_bbox
    }
    return {
        'checks':checks,
        'readback_effective_length_um':actual_l,
        'readback_effective_width_um':actual_w,
        'readback_square_count':actual_l/actual_w,
        'nominal_resistance_ohm':rs*actual_l/actual_w,
        'total_TiN_length_um':h.bbox().width()*layout.dbu,
        'metal_contact_edges_local_x_um':[signal_edge,ground_edge],
        'terminal_origin_in_ours_um':[inst.dcplx_trans.disp.x,inst.dcplx_trans.disp.y],
        'contact_areas':contacts,
        'metal_components':ids,
        'ground_bridge_optical_overlap':{
            'before_area_um2':old_overlap.area()*layout.dbu**2,
            'after_area_um2':new_overlap.area()*layout.dbu**2,
            'after_polygon_count':new_overlap.size(),
            'after_bboxes_um':[bb(p.bbox()) for p in new_overlap.each()],
            'closed':False,
            'note':'仍存在金属覆盖光芯；本次不改变光路或虚构光学验证。'
        }
    }
built=audit(l)
if not all(built['checks'].values()): raise RuntimeError(json.dumps(built,ensure_ascii=False))
l.write(payload['output'])
verify=pya.Layout(); verify.read(payload['output'])
result=audit(verify)
if not all(result['checks'].values()): raise RuntimeError('写出后的GDS回读检查失败')
result.update({
    'status':'resistor_dimensions_corrected_not_tapeout_ready',
    'source_gds':payload['source'],'output_gds':payload['output'],
    'sheet_resistance_ohm_per_square':rs,
    'target_ohm':target,'top_cell':verify.top_cell().name,
    'changed_layers':['41/0','42/0'],'resistor_geometry_pass':True,
    'sheet_resistance_is_user_authorized_estimate':True,
    'contact_and_temperature_model':'ideal contacts; no TCR/heating calibration',
    'rf_signoff':False,'tapeout_ready':False,
    'notes':['未改YSJ、光路、MUX、O型输入、15mm电极本体和LN2。',
             '地桥及地侧接触右移，以保留金属接触边缘间31.25um电阻区。',
             '地桥跨光芯仍未关闭；没有启动光学、HFSS或容差求解。']
})
json.dumps(result,ensure_ascii=False)
"""

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',required=True,type=Path)
    p.add_argument('--output-dir',required=True,type=Path)
    p.add_argument('--sheet-resistance',type=float,default=20.)
    p.add_argument('--target',type=float,default=50.)
    p.add_argument('--width',type=float,default=12.5)
    a=p.parse_args()
    source=a.source.resolve(); dest=a.output_dir.resolve()
    if not source.is_file(): raise FileNotFoundError(source)
    if min(a.sheet_resistance,a.target,a.width)<=0: raise ValueError('参数必须为正')
    if a.width<2: raise ValueError('HTR线宽低于本轮采用的2um下限')
    if dest.exists(): raise FileExistsError(dest)
    before=hashlib.sha256(source.read_bytes()).hexdigest()
    c=KLinkClient(); c.connect()
    try:
        dest.mkdir(parents=True)
        output=dest/'YSJ与四程15mm_M04_O型输入_HTR有效50R_候选.gds'
        payload={'source':str(source),'output':str(output),'sheet_resistance':a.sheet_resistance,
                 'target':a.target,'width':a.width}
        r=c.call('exec.python',{'code':'payload='+repr(payload)+'\n'+PAYLOAD_CODE})
        if r.get('exception'): raise RuntimeError(r['exception'])
        result=json.loads(r['return_value'])
        if hashlib.sha256(source.read_bytes()).hexdigest()!=before: raise RuntimeError('源GDS发生变化')
        result['source_sha256']=before
        result['output_sha256']=hashlib.sha256(output.read_bytes()).hexdigest()
        mismatch=[]
        for load in (24.,45.,50.,55.):
            gamma=(load-50.)/(load+50.)
            mismatch.append({'R_ohm':load,'Z0_ohm':50.,'Gamma':gamma,
                             'reflected_power_percent':100*gamma**2,
                             'VSWR':(1+abs(gamma))/(1-abs(gamma))})
        result['ideal_load_mismatch_budget']=mismatch
        (dest/'电阻有效长度与接触回读检查.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
        c.call('layout.show_file',{'path':str(output),'mode':'new'})
        c.call('view.show_cell',{'cell':result['top_cell'],'zoom_fit':True})
        c.call('view.hier_levels',{'min':0,'max':20})
        c.call('layer.set_visible',{'layers':['41/0','42/0','20/0'],'visible':True})
        for layer,color in [('41/0','#E64A19'),('42/0','#D5A12C'),('20/0','#C62858')]:
            c.call('layer.set_style',{'layer':layer,'fill_color':color,'frame_color':color,
                                     'dither_pattern':0,'line_width':2})
        c.call('layer.set_visible',{'layers':['10/0','10/1','10/2','20/1','20/2','21/2',
                                            '56/30','70/30','100/0','100/30','101/0'],'visible':False})
        c.call('view.screenshot',{'mode':'path','path':str(dest/'终端有效电阻修正_地桥跨光芯仍待处理.png'),
                                'bbox_um':[5340,1460,5500,1780],'width_px':1400,'height_px':1000})
        c.call('view.zoom_box',{'bbox_um':[5360,1580,5470,1660]})
        print(json.dumps(result,ensure_ascii=False,indent=2))
    finally: c.close()

if __name__=='__main__': main()
