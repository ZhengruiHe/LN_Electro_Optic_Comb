"""计算0.7um无源LN rib沿晶体Y/Z传播的TE0群折射率，供欧拉回路时延预算。"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from mode_pdk_sweep import (DEFAULT_CONFIG, add_polygon, add_rect, load_config,
                            load_lumapi, mode_localization_metrics, scalar,
                            zelmon_ln_indices)

ROOT=Path(__file__).resolve().parents[2]


def index_string(wavelength_um: float, propagation_axis: str) -> str:
    no,ne=zelmon_ln_indices(np.asarray([wavelength_um]))
    if propagation_axis=='Y':
        # 仿真横向x=晶体Z，竖直y=晶体X，传播z=晶体Y。
        values=(ne[0],no[0],no[0])
    elif propagation_axis=='Z':
        # 仿真横向x=晶体Y，竖直y=晶体X，传播z=晶体Z。
        values=(no[0],no[0],ne[0])
    else:
        raise ValueError('传播晶轴必须为Y或Z')
    return ';'.join(f'{value:.12g}' for value in values)


def build(mode,cfg,wavelength_um,axis):
    mode.switchtolayout(); mode.deleteall()
    stack,wg,materials,fde=cfg['stack_um'],cfg['waveguide_um'],cfg['material_models'],cfg['fde']
    half=fde['x_span_um']/2
    sin_top=stack['sin']; ln_bottom=sin_top+stack['ln_sin_interlayer_oxide']
    ln_top=ln_bottom+stack['ln']; residual=.2; angle=70.; width=.7; platform=7.2
    top_oxide=ln_top+stack['top_cladding']
    add_rect(mode,'Bottom_Oxide',-half,half,fde['y_min_um'],0,materials['sio2'],4)
    add_rect(mode,'SiN_Removed_Oxide_Fill',-half,half,0,sin_top,materials['sio2'],3)
    add_rect(mode,'LN_SiN_Interlayer_Oxide',-half,half,sin_top,ln_bottom,materials['sio2'],4)
    offset=2*residual/math.tan(math.radians(angle))
    idx=index_string(wavelength_um,axis)
    add_polygon(mode,'LN_Residual_Slab',np.asarray([[-(platform+offset)/2,ln_bottom],
        [(platform+offset)/2,ln_bottom],[platform/2,ln_bottom+residual],[-platform/2,ln_bottom+residual]]),
        '<Object defined dielectric>',1,idx)
    slab_top=ln_bottom+residual
    add_polygon(mode,'LN_Rib',np.asarray([[-(width+offset)/2,slab_top],[(width+offset)/2,slab_top],
        [width/2,ln_top],[-width/2,ln_top]]),'<Object defined dielectric>',1,idx)
    add_rect(mode,'Top_Oxide',-half,half,ln_bottom,top_oxide,materials['sio2'],4)
    mode.addfde(); mode.set('solver type','2D Z normal')
    mode.set('x min',-half*1e-6); mode.set('x max',half*1e-6)
    mode.set('y min',fde['y_min_um']*1e-6); mode.set('y max',fde['y_max_um']*1e-6)
    for side in ('x min','x max','y min','y max'):
        mode.set(side+' bc','PML')
    mode.set('mesh cells x',fde['mesh_cells_x']); mode.set('mesh cells y',fde['mesh_cells_y'])
    mode.set('wavelength',wavelength_um*1e-6); mode.set('number of trial modes',12)
    mode.set('calculate group index',False); mode.set('detailed dispersion calculation',False)
    mode.set('simulation temperature',cfg['temperature_c'])


def solve(mode,cfg,wavelength_um,axis):
    build(mode,cfg,wavelength_um,axis)
    count=int(mode.findmodes())
    candidates=[]
    for number in range(1,count+1):
        path=f'FDE::data::mode{number}'
        te=float(scalar(mode.getdata(path,'TE polarization fraction')).real)
        if te<.55:
            continue
        loc=mode_localization_metrics(mode,path,1.5)
        candidates.append({'solver_mode':number,'neff':scalar(mode.getdata(path,'neff')).real,
                           'te_fraction':te,**loc})
    candidates.sort(key=lambda row:row['neff'],reverse=True)
    eligible=[row for row in candidates if row['central_energy_fraction']>=.7 and abs(row['peak_x_um'])<=.5]
    if not eligible:
        raise RuntimeError(f'{axis}轴{wavelength_um*1000:g}nm未找到中心局域TE0')
    return eligible[0],candidates


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--step-nm',type=float,default=5.)
    a=p.parse_args()
    if a.step_nm<=0:
        raise ValueError('波长步长必须为正')
    a.output_dir.mkdir(parents=True,exist_ok=False)
    cfg=load_config(DEFAULT_CONFIG); api=load_lumapi(cfg)
    result={'status':'本机MODE名义方向群折射率；不是弯曲模式或工艺容差结果',
            'width_um':.7,'ln_total_um':.4,'rib_height_um':.2,'residual_slab_um':.2,
            'platform_top_width_um':7.2,'sidewall_angle_deg_from_horizontal':70,
            'wavelength_nm':1550.,'step_nm':a.step_nm,'axes':{}}
    with api.MODE(hide=True) as mode:
        result['software_version']=mode.version()
        for axis in ('Y','Z'):
            samples=[]
            for wavelength_nm in (1550-a.step_nm,1550.,1550+a.step_nm):
                selected,candidates=solve(mode,cfg,wavelength_nm/1000,axis)
                samples.append({'wavelength_nm':wavelength_nm,'selected':selected,
                                'candidate_count':len(candidates)})
            lower,center,upper=[row['selected'] for row in samples]
            derivative=(upper['neff']-lower['neff'])/(2*a.step_nm/1000)
            ng=center['neff']-1.55*derivative
            result['axes'][axis]={'propagation_crystal_axis':axis,'group_index':ng,
                'neff_1550':center['neff'],'te_fraction_1550':center['te_fraction'],
                'central_energy_fraction_1550':center['central_energy_fraction'],
                'peak_x_um_1550':center['peak_x_um'],'samples':samples}
    result['directional_difference']=result['axes']['Z']['group_index']-result['axes']['Y']['group_index']
    result['limitations']=['只求晶体Y/Z两个主方向；任意角度暂用cos²插值做工程时延积分',
                           '无金属、无弯曲、无粗糙度；不等于欧拉弯损耗验证',
                           '20/1按7.2um名义平台模型；v10/v11局部窗口变化未回填']
    path=a.output_dir/'0p7um_rib_两个晶向群折射率.json'
    path.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))


if __name__=='__main__':
    main()
