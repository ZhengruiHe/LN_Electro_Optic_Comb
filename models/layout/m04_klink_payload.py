import json, math, pya

source = pya.Layout(); source.read(payload['source_gds'])
layout = pya.Layout(); layout.read(payload['source_gds'])
top = layout.top_cell(); mux = layout.cell(payload['mux_cell'])
if top is None or mux is None or len(layout.top_cells()) != 1 or abs(layout.dbu - .001) > 1e-12:
    raise RuntimeError('顶层、MUX或格点不符合要求')

def smooth(u):
    return u*u*(3.0-2.0*u)

def progress(u):
    a = u*u; b = a + 2.3333333333*(1-u)*(1-u)
    return a/b if b else 0.0

def make_points(section, count=240):
    out=[]
    for i in range(count+1):
        limits={'active':(-100.0,0.0),'body':(0.0,750.0),'external':(750.0,800.0)}
        x0,x1=limits[section]; x=x0+(x1-x0)*i/count
        if section=='active':
            u=(x+100)/100; q=smooth(u); mw=1.33+.12*q; aw=.35; mc=-1.95-.025*q; ac=1.40+.025*q
        elif section=='body':
            u=x/750; q=progress(u); mw=1.45-.33*q; aw=.35+.07*q; gap=.8+1.7*(.5+.5*math.cos(2*math.pi*u)); mc=-gap/2-mw/2; ac=gap/2+aw/2
        else:
            # 外部两支路必须回到现有路由的端口中心：主脊约 -1.80 µm、辅助脊约 1.45 µm。
            # 之前把它们错误地平移到 0.15/3.40 µm，造成 MUX 边界后的实际断口。
            u=(x-750)/50; q=smooth(u); mw=1.12+(.70-1.12)*q; aw=.42+(.70-.42)*q; mc=-1.81+.01*q; ac=1.46-.01*q
        out.append((x,mw,mc,aw,ac))
    return out

path=make_points('active')[:-1]+make_points('body')[1:-1]+make_points('external')[1:]
main=[(x,mw,mc) for x,mw,mc,aw,ac in path]; aux=[(x,aw,ac) for x,mw,mc,aw,ac in path]

def poly(items):
    lo=[pya.DPoint(x,c-w/2) for x,w,c in items]; hi=[pya.DPoint(x,c+w/2) for x,w,c in reversed(items)]
    return pya.DPolygon(lo+hi).to_itype(layout.dbu)

core_layer=layout.layer(20,0); clad_layer=layout.layer(20,1); trench_layer=layout.layer(10,2)
for li in (core_layer,clad_layer,trench_layer): mux.shapes(li).clear()
mux.shapes(core_layer).insert(poly(main)); mux.shapes(core_layer).insert(poly(aux))
for li,width in ((clad_layer,7.2),(trench_layer,17.2)):
    region=pya.Region()
    region += pya.Region(pya.DPath([pya.DPoint(x,c) for x,_,c in main],width).to_itype(layout.dbu))
    region += pya.Region(pya.DPath([pya.DPoint(x,c) for x,_,c in aux],width).to_itype(layout.dbu))
    mux.shapes(li).insert(region)
mux.name='MUX_TE01_A70_M04_A100_E50_BODY750_DRAFT'; top.name='YSJ_EO4P_15MM_M04_HORIZONTAL_DRAFT'
core=pya.Region(mux.shapes(core_layer)).merged(); clad=pya.Region(mux.shapes(clad_layer)).merged(); trench=pya.Region(mux.shapes(trench_layer)).merged()
samples=[{'x_um':x,'main_width_um':mw,'aux_width_um':aw,'edge_gap_um':ac-aw/2-(mc+mw/2)} for x,mw,mc,aw,ac in path]
checks={'one_top':len(layout.top_cells())==1,'core_valid':not core.is_empty(),'core_inside_clad':(core-clad).is_empty(),'core_inside_trench':(core-trench).is_empty(),'minimum_gap_gt_0p6':min(s['edge_gap_um'] for s in samples)>.6,'core_width_markers':core.width_check(300).size(),'core_space_markers':core.space_check(300).size(),'etch_space_markers':(clad-core).space_check(300).size()}
if not all(checks[k] for k in ('one_top','core_valid','core_inside_clad','core_inside_trench','minimum_gap_gt_0p6')): raise RuntimeError('M04几何检查失败:'+json.dumps(checks))
layout.write(payload['output_gds'])
result={'status':'M04_KLink_layout_candidate','output_gds':payload['output_gds'],'source_gds':payload['source_gds'],'top_cell':top.name,'mux_cell':mux.name,'checks':checks,'sampled_geometry':samples,'changed_layers':['20/0','20/1','10/2'],'new_LN2_layer':False,'fixed_length_um':750.0,'coupling_length_um':550.0,'full_chain_pass':False,'tapeout_signoff':False}
open(payload['report'],'w',encoding='utf-8').write(json.dumps(result,ensure_ascii=False,indent=2)); result
