"""Local optical metal absorption estimate for an LN rib under a TiN bridge."""
from __future__ import annotations
import argparse, json, math, sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'models/optical'))
from mode_pdk_sweep import load_config, load_lumapi, ln_index_string, scalar, mode_localization_metrics

def rect(mode,name,x0,x1,y0,y1,material,mesh=3,index=None):
    mode.addrect(); mode.set('name',name); mode.set('x min',x0*1e-6); mode.set('x max',x1*1e-6)
    mode.set('y min',y0*1e-6); mode.set('y max',y1*1e-6); mode.set('z span',1e-6)
    mode.set('material',material); mode.set('override mesh order from material database',True); mode.set('mesh order',mesh)
    if index is not None: mode.set('index',index)

def poly(mode,name,pts,material,mesh=1,index=None):
    mode.addpoly(); mode.set('name',name); mode.set('vertices',np.asarray(pts)*1e-6); mode.set('z span',1e-6)
    mode.set('material',material); mode.set('override mesh order from material database',True); mode.set('mesh order',mesh)
    if index is not None: mode.set('index',index)

def build(mode,cfg,with_metal,metal_material):
    mode.switchtolayout(); mode.deleteall()
    s=cfg['stack_um']; w=cfg['waveguide_um']; f=cfg['fde']
    top_width=1.33; etch=.2; ang=math.radians(70.)
    rib_bottom=top_width+2*etch/math.tan(ang)
    slab_bottom=7.2+2*(s['ln']-etch)/math.tan(ang)
    ln_bottom=s['sin']+s['ln_sin_interlayer_oxide']; ln_top=ln_bottom+s['ln']; slab_top=ln_bottom+s['ln']-etch
    half=12.
    rect(mode,'Bottom_Oxide',-half,half,-2.5,0.,cfg['material_models']['sio2'],4)
    rect(mode,'SiN_Removed_Oxide_Fill',-half,half,0.,s['sin'],cfg['material_models']['sio2'],3)
    rect(mode,'LN_SiN_Interlayer_Oxide',-half,half,s['sin'],ln_bottom,cfg['material_models']['sio2'],4)
    poly(mode,'LN_Residual_Slab',[[-slab_bottom/2,ln_bottom],[slab_bottom/2,ln_bottom],[7.2/2,slab_top],[-7.2/2,slab_top]],'<Object defined dielectric>',1,ln_index_string(1.55))
    poly(mode,'LN_Rib',[[-rib_bottom/2,slab_top],[rib_bottom/2,slab_top],[top_width/2,ln_top],[-top_width/2,ln_top]],'<Object defined dielectric>',1,ln_index_string(1.55))
    rect(mode,'Top_Oxide',-half,half,ln_top,ln_top+1.,cfg['material_models']['sio2'],4)
    if with_metal:
        rect(mode,'M1_Bridge_20um',-10.,10.,ln_top+1.,ln_top+2.,metal_material,1)
    mode.addfde(); mode.set('solver type','2D Z normal'); mode.set('x min',-half*1e-6); mode.set('x max',half*1e-6)
    mode.set('y min',-2.5e-6); mode.set('y max',3.2e-6); mode.set('x min bc','PML'); mode.set('x max bc','PML'); mode.set('y min bc','PML'); mode.set('y max bc','PML')
    mode.set('mesh cells x',360); mode.set('mesh cells y',240); mode.set('wavelength',1.55e-6); mode.set('number of trial modes',12)
    mode.set('use max index',False); mode.set('n',1.8); mode.set('calculate group index',False)
    return {'top_width_um':top_width,'etch_depth_um':etch,'sidewall_deg_from_horizontal':70.,'top_cladding_um':1.,'metal_thickness_um':1.,'metal_width_um':20. if with_metal else 0.}

def solve(mode,cfg,with_metal,metal_material):
    geom=build(mode,cfg,with_metal,metal_material); n=int(mode.findmodes()); rows=[]
    for i in range(1,n+1):
        p=f'FDE::data::mode{i}'
        te=float(scalar(mode.getdata(p,'TE polarization fraction')).real)
        ne=scalar(mode.getdata(p,'neff'))
        loss=float(scalar(mode.getdata(p,'loss')).real)
        loc=mode_localization_metrics(mode,p,1.5)
        rows.append({'solver_mode':i,'neff_real':ne.real,'neff_imag':ne.imag,'loss_db_per_m':loss,'loss_db_per_cm':loss/100.,'te_fraction':te,**loc})
    te_rows=[r for r in rows if r['te_fraction']>=.55]
    te_rows.sort(key=lambda r:(-r['central_energy_fraction'],-r['te_fraction']))
    if not te_rows: raise RuntimeError('未找到TE候选')
    return {'with_tin_metal':with_metal,'geometry':geom,'modes':rows,'selected_TE0':te_rows[0]}

def main():
    p=argparse.ArgumentParser(); p.add_argument('--output-dir',type=Path,required=True); p.add_argument('--metal-material',default='TiN - Palik'); a=p.parse_args()
    a.output_dir=a.output_dir.resolve(); a.output_dir.mkdir(parents=True,exist_ok=False)
    cfg=load_config(ROOT/'models/optical/mode_pdk_config.json'); api=load_lumapi(cfg)
    result={'status':'running','wavelength_nm':1550.,'material_file':r'C:\Program Files\Netease\GameViewer\Download\LN材料参数.mdf','metal_material':a.metal_material,'source_config':str((ROOT/'models/optical/mode_pdk_config.json').resolve()),'cases':[]}
    with api.MODE(hide=True) as mode:
        result['solver_version']=mode.version(); result['material_exists_metal']=bool(mode.materialexists(a.metal_material)); result['material_exists_sio2']=bool(mode.materialexists('SiO2 (Glass) - Palik'))
        if not result['material_exists_metal']: raise RuntimeError('MODE材料库中没有指定金属: '+a.metal_material)
        for flag in (False,True):
            case=solve(mode,cfg,flag,a.metal_material); result['cases'].append(case)
            mode.save(str((a.output_dir/('metal_present.lms' if flag else 'metal_absent.lms')).resolve()))
    bare=result['cases'][0]['selected_TE0']; metal=result['cases'][1]['selected_TE0']
    result['incremental_loss_db_per_cm']=metal['loss_db_per_cm']-bare['loss_db_per_cm']
    result['incremental_loss_db_per_um']=result['incremental_loss_db_per_cm']/1e4
    result['metal_bridge_length_um']=20.
    result['estimated_bridge_extra_loss_db']=result['incremental_loss_db_per_um']*20.
    result['status']='completed_local_mode_absorption_estimate'
    result['interpretation']='局部连续TiN覆盖的单位长度吸收；有限20um桥边缘散射、真实金属位置和工艺参数未包含'
    (a.output_dir/'金属吸收损耗结果.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))
if __name__=='__main__': main()
