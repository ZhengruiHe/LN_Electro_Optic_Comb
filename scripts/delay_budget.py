"""先求10GHz方向群时延和欧拉回路长度；只写预算，不生成或修改GDS。"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'models/layout'))

from euler_routes import (  # noqa: E402
    euler_compact_finish_return_route,
    euler_input_route,
    euler_loop2,
)
from route_geometry import check_route_network  # noqa: E402

C_UM_PER_PS=299.792458


def route_delay(points,ng_y,ng_z):
    length=weighted_y=weighted_z=0.0
    for first,second in zip(points,points[1:]):
        dx=second[0]-first[0]; dy=second[1]-first[1]
        ds=math.hypot(dx,dy)
        if ds==0:
            continue
        cos2=(dx/ds)**2; sin2=(dy/ds)**2
        length+=ds; weighted_y+=ds*cos2; weighted_z+=ds*sin2
    delay=(ng_y*weighted_y+ng_z*weighted_z)/C_UM_PER_PS
    return {'centerline_length_um':length,'crystal_Y_weighted_length_um':weighted_y,
            'crystal_Z_weighted_length_um':weighted_z,'directional_delay_ps':delay,
            'equivalent_average_group_index':delay*C_UM_PER_PS/length}


def builders():
    return [
        lambda length:euler_compact_finish_return_route(
            x_right=18150,x_left=850,y_start=29.85,y_end=33.4,
            target_length=length,outward_sign=1),
        euler_loop2,
        lambda length:euler_compact_finish_return_route(
            x_right=18150,x_left=850,y_start=-29.85,y_end=-33.4,
            target_length=length,outward_sign=-1),
    ]


def solve_length(builder,initial_length,target_delay,ng_y,ng_z):
    length=initial_length
    history=[]
    for _ in range(4):
        points,bends=builder(length)
        delay=route_delay(points,ng_y,ng_z)
        error=target_delay-delay['directional_delay_ps']
        history.append({'trial_length_um':length,'delay_ps':delay['directional_delay_ps'],'error_ps':error})
        length+=error*C_UM_PER_PS/ng_y
    points,bends=builder(length)
    delay=route_delay(points,ng_y,ng_z)
    return length,points,bends,delay,history


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--ng',type=Path,default=ROOT/'results/optical/欧拉时延输入_0p7um方向群折射率_v1/0p7um_rib_两个晶向群折射率.json')
    parser.add_argument('--baseline',type=Path,default=ROOT/'results/layout/四程电光梳_10GHz_15mm_3T_3.5T_3T_名义_检查.json')
    parser.add_argument('--output-dir',type=Path,required=True)
    args=parser.parse_args()
    ng_data=json.loads(args.ng.read_text(encoding='utf-8'))
    baseline=json.loads(args.baseline.read_text(encoding='utf-8'))
    ng_y=ng_data['axes']['Y']['group_index']; ng_z=ng_data['axes']['Z']['group_index']
    if baseline['loop_periods']!=[3.0,3.5,3.0] or baseline['electrode']['active_length_um']!=15000:
        raise ValueError('预算仅适用于当前10GHz、15mm、3T/3.5T/3T基线')
    args.output_dir.mkdir(parents=True,exist_ok=False)
    rows=[]; route_points={}; route_meta={}
    for index,(base,builder) in enumerate(zip(baseline['delay_budget'],builders()),1):
        passive_target=(base['target_delay_ps']-base['active_delay_ps']-base['two_mux_bodies_delay_ps']
                        -base['two_active_tapers_delay_ps']-base['two_external_tapers_delay_ps'])
        old_length=base['passive_route_target_um']
        old_points,old_bends=builder(old_length)
        old_directional=route_delay(old_points,ng_y,ng_z)
        new_length,points,bends,delay,history=solve_length(builder,old_length,passive_target,ng_y,ng_z)
        total=(base['active_delay_ps']+base['two_mux_bodies_delay_ps']+base['two_active_tapers_delay_ps']
               +base['two_external_tapers_delay_ps']+delay['directional_delay_ps'])
        residual=total-base['target_delay_ps']
        old_total=(base['active_delay_ps']+base['two_mux_bodies_delay_ps']+base['two_active_tapers_delay_ps']
                   +base['two_external_tapers_delay_ps']+old_directional['directional_delay_ps'])
        row={'loop':index,'periods':base['periods'],'target_total_delay_ps':base['target_delay_ps'],
             'fixed_active_delay_ps':base['active_delay_ps'],'fixed_two_mux_delay_ps':base['two_mux_bodies_delay_ps'],
             'fixed_two_active_taper_delay_ps':base['two_active_tapers_delay_ps'],
             'fixed_two_external_taper_delay_ps':base['two_external_tapers_delay_ps'],
             'required_passive_route_delay_ps':passive_target,'old_circular_budget_length_um':old_length,
             'old_length_with_euler_directional_delay_ps':old_directional['directional_delay_ps'],
             'old_length_total_delay_ps':old_total,'old_length_phase_error_deg_at_10GHz':(old_total-base['target_delay_ps'])*3.6,
             'solved_euler_route_length_um':new_length,'length_correction_um':new_length-old_length,
             **delay,'total_delay_ps':total,'residual_delay_ps':residual,
             'residual_phase_deg_at_10GHz':residual*3.6,'euler_bend_count':len(bends['bends']),
             'minimum_radius_um':bends['minimum_radius_um'],'euler_fraction':bends['euler_fraction'],
             'solver_history':history}
        rows.append(row); route_points['loop'+str(index)]=points; route_meta['loop'+str(index)]=bends
    input_points,input_meta=euler_input_route()
    route_points['input']=input_points
    route_points['output']=[[18150.,-26.6],[18700.,-26.6]]
    route_meta['input']=input_meta
    topology=check_route_network(
        route_points,
        [{'routes':['input','loop2'],'center_um':[300.,250.],'angle_deg':90.}],
    )
    delay_closed=all(abs(row['residual_phase_deg_at_10GHz'])<1e-4 for row in rows)
    blocking=[]
    if not topology['topology_screen_pass']:
        blocking.append('输入、输出或三条回路的中心线拓扑筛查未通过')
    result={'status':'欧拉版图生成前的10GHz名义方向时延预算；非完整级联签核',
            'frequency_ghz':10.,'period_ps':100.,'target_periods':[3.,3.5,3.],
            'passive_rib':{'width_um':.7,'ng_crystal_Y':ng_y,'ng_crystal_Z':ng_z,
                           'angular_interpolation':'ng(theta)=ngY*cos(theta)^2+ngZ*sin(theta)^2'},
            'euler':{'fraction':.4,'minimum_radius_um':80.,'geometry_only':True},
            'loops':rows,'route_points_um':route_points,'route_metadata':route_meta,
            'topology_screen':topology,
            'delay_closed':delay_closed,
            'layout_generation_allowed':delay_closed and topology['topology_screen_pass'],
            'blocking_before_gds':blocking,
            'limitations':['任意角度群折射率采用Y/Z两主轴cos²插值，是工程模型，不是弯曲本征模求解',
                           'MUX群时延沿用旧7.2um EME本体；v10/v11窗口变化尚未回填',
                           'loop2把黑盒交叉器及0.7↔1.2um拉锥按0.7um中心线积分；真实复传输群时延未知',
                           '没有计算欧拉弯损耗、模式串扰、工艺容差或21层影响']}
    (args.output_dir/'欧拉回路_10GHz方向时延预算.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    fields=[key for key in rows[0] if key not in ('solver_history',)]
    with (args.output_dir/'欧拉回路_10GHz方向时延预算.csv').open('w',encoding='utf-8-sig',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=fields); writer.writeheader()
        writer.writerows([{k:v for k,v in row.items() if k in fields} for row in rows])
    print(json.dumps({"layout_generation_allowed":result['layout_generation_allowed'],
        "loops":[{k:r[k] for k in ('loop','solved_euler_route_length_um','length_correction_um','total_delay_ps',
                  'residual_phase_deg_at_10GHz','minimum_radius_um')} for r in rows]},ensure_ascii=False,indent=2))


if __name__=='__main__':
    main()
