"""固定器件的周期匹配布线：整体移动左侧路由束，回标三条回路。

默认10GHz，3T/3.5T/3T。MUX耦合核心、原尾端、金属和晶向不改；
有源侧接续100um，外侧50um。只写新预算目录，不写GDS。
"""
import argparse
import json
import math
from pathlib import Path
from shapely.geometry import LineString, Point, box
from shapely.ops import substring
from route_15mm_left_delays import PathBuilder, build_left
from route_geometry import check_route_network
from euler_routes import bend_displacement, radius_for_semicircle_displacement
from route_compact_254_trial import inputs, fixed_delays
from delay_budget import route_delay
from search_fixed_block_placement import geometry

ROOT=Path(__file__).resolve().parents[2]


def build(active=100.,external=50.,periods=(3.,3.5,3.),spread_access=False):
    if not 100<=active<=200 or not 50<=external<=100:
        raise ValueError('本轮只允许100..200um有源接续、50..100um外接续，不改原辅助尾端')
    if periods[0]%1 or periods[2]%1 or abs(periods[1]%1-.5)>1e-9:
        raise ValueError('当前极性序列要求整数、半整数、整数周期')
    old,_=inputs();ny=old['passive_rib']['ng_crystal_Y'];nz=old['passive_rib']['ng_crystal_Z']
    fixed=fixed_delays(active,external)
    base=ROOT/'results/layout/_历史归档/迭代版本_20260912/合并YSJ_版图确认_20260909_v1'
    study=json.loads((base/'15mm重新布局可行性_v9_四出口完整局部检查/15mm上方模块可行性.json').read_text(encoding='utf-8'))
    dx,dy=-11620.,1620.
    shrink=(200-active)+(100-external)
    left_shift=shrink
    xl=950+dx+shrink
    configs=[
        dict(down_col=-10845.+left_shift,entry_y=650.,return_col=-10200.+left_shift,u_count=1,
             fold_left=None,final_rectangular=True,final_start_y=1050.,final_radius=83.,final_x=-10728.638053930058+left_shift),
        dict(down_col=-10865.+left_shift,entry_y=300.,return_col=-8000.,u_count=3,
             fold_left=-7600.,final_rectangular=True,final_start_y=1180.,final_radius=80.,final_x=-10720.+left_shift),
        dict(down_col=-10890.+left_shift,entry_y=-400.,return_col=-10530.+left_shift,u_count=1,
             fold_left=None,final_rectangular=False,final_start_y=None,final_radius=80.,final_x=-10671.+left_shift),
    ]
    if spread_access:
        configs[0].update(down_col=-10730.,final_radius=130.,final_x=xl-10.)
        configs[1].update(down_col=-10760.,final_radius=100.,final_x=xl-10.)
        configs[2].update(down_col=-10820.)
    right_prefixes={}
    if spread_access:
        # 起始端口不动，只拉开回程长直段。顺序沿原嵌套关系保持。
        first_lane=dy+29.85+bend_displacement(math.pi,80.)[1]
        for name,y,lane in [('loop1',dy+29.85,first_lane),('loop2',dy+26.6,first_lane+20.),
                             ('output',dy-26.6,first_lane+40.),('loop3',dy-29.85,first_lane+60.)]:
            p=PathBuilder([[18050.+dx-shrink,y]],0.);p.straight(x=p.points[-1][0]+5.)
            p.bend(math.pi,radius_for_semicircle_displacement(lane-y),'拉开间距的右端欧拉折返')
            right_prefixes[name]=p
    routes={};metas={}
    ends=[dy+33.4,dy-30.15,dy-33.4]
    for i,(cfg,end) in enumerate(zip(configs,ends),1):
        prefix=(right_prefixes[f'loop{i}'].points if spread_access else
                [[x+dx-shrink,y+dy] for x,y in study['right_escape_local_points_um'][f'loop{i}']])
        fold=[-8800.,-6500.,-10100.][i-1]+left_shift
        target=periods[i-1]*100.-fixed[i-1]
        history=[]
        for _ in range(4):
            pts,bends=build_left(prefix,i,fold,xl,end,cfg)
            delay=route_delay(pts,ny,nz);error=target-delay['directional_delay_ps']
            history.append(dict(fold_x_um=fold,error_ps=error))
            fold+=error*299.792458/ny/(cfg['u_count']+1)
        pts,bends=build_left(prefix,i,fold,xl,end,cfg)
        routes[f'loop{i}']=pts
        metas[f'loop{i}']={**route_delay(pts,ny,nz),'bends':bends,'config':cfg,'right_fold_x_um':fold,
            'target_total_ps':periods[i-1]*100.,'fixed_delay_ps':fixed[i-1],'target_passive_delay_ps':target,'history':history}
    q=bend_displacement(math.pi/2,80)[0]
    input_col=-8340.
    incoming=PathBuilder([[-10380.,-1700.]],0.)
    incoming.straight(x=-10280.);incoming.straight(x=input_col+500.)
    incoming.bend(math.pi,80.,'输入底部折返');incoming.straight(x=input_col+q)
    incoming.bend(-math.pi/2,80.,'输入上行');incoming.straight(y=1120.-q)
    input_radius=115. if spread_access else 80.
    input_q=bend_displacement(math.pi/2,input_radius)[0]
    incoming.bend(math.pi/2,80.,'输入向西');incoming.straight(x=xl-10. if spread_access else -10729.638053930058+left_shift)
    incoming.bend(-math.pi/2,input_radius,'输入MUX前上转');incoming.straight(y=dy+30.15-input_q)
    incoming.bend(-math.pi/2,input_radius,'输入接MUX');incoming.straight(x=xl)
    routes['input']=incoming.points;metas['input']={**route_delay(incoming.points,ny,nz),'bends':incoming.bends}
    prefix=(right_prefixes['output'].points if spread_access else
            [[x+dx-shrink,y+dy] for x,y in study['right_escape_local_points_um']['output']])
    out=PathBuilder(prefix,math.pi);out.straight(x=(-10790. if spread_access else -10877.+left_shift)+q)
    out.bend(math.pi/2,80.,'输出向下');out.straight(y=100.+q);out.bend(math.pi/2,80.,'输出向东')
    out.straight(x=-9500.+left_shift);out.bend(-math.pi,radius_for_semicircle_displacement(600.),'输出折返')
    out.straight(x=-10280.);out.straight(x=-10380.)
    routes['output']=out.points;metas['output']={**route_delay(out.points,ny,nz),'bends':out.bends}
    topology=check_route_network(routes);crossings=[]
    for c in topology['unexpected_intersections']:
        r=dict(c)
        if 'x_um' in r:
            point=Point(r['x_um'],r['y_um']);r['straight_245um_each_route']={}
            for name in r['routes']:
                line=LineString(routes[name]);s=line.project(point);sub=substring(line,max(0,s-122.5),min(line.length,s+122.5))
                r['straight_245um_each_route'][name]=(s>=122.5 and s+122.5<=line.length and abs(sub.length-math.dist(sub.coords[0],sub.coords[-1]))<1e-5)
        crossings.append(r)
    source=json.loads((base/'YSJ源版图几何.json').read_text(encoding='utf-8'))
    old_core=geometry(source,(20,0));old_sin=geometry(source,(10,0));old_etch=geometry(source,(21,2))
    checks={}
    for name,pts in routes.items():
        core=LineString(pts).buffer(.35,cap_style=2);window=LineString(pts).buffer(8.6,cap_style=2)
        checks[name]={'YSJ_core_overlap_um2':core.intersection(old_core).area,
                     'SiN_removal_over_YSJ_SiN_um2':window.intersection(old_sin).area,
                     'LN1_over_YSJ_LN2_etch_um2':core.intersection(old_etch).area,
                     'window_inside_block':box(-10900,-1900,10900,1900).covers(window)}
    result=dict(status='周期匹配路由预算；须回读GDS',electrode_length_um=15000,shift_um=[dx,dy],
        active_taper_um=active,external_taper_um=external,loop_periods=list(periods),
        left_route_shift_um=left_shift,crossing_bus_x_um=-10530.+left_shift,
        spread_access=spread_access,right_return_lane_pitch_um=20. if spread_access else None,
        window_merge_caps=[] if spread_access else [{'layer':[20,1],
            'bounds_um':[-10832.+left_shift,1872.,-10822.+left_shift,1879.],
            'reason':'把两条回程公共窗口的收尖凹口截成宽端，不改核心'}],
        routes=routes,delays=metas,topology=topology,crossing_details=crossings,shape_checks=checks,
        input_and_final_output_included=True,full_layout_pass=False)
    return result


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output-dir',type=Path,required=True)
    ap.add_argument('--active',type=float,default=100.);ap.add_argument('--external',type=float,default=50.)
    ap.add_argument('--periods',type=float,nargs=3,default=[3.,3.5,3.])
    ap.add_argument('--spread-access',action='store_true',help='拉开接入段和回程长直段间距；不改MUX、电极')
    a=ap.parse_args()
    r=build(a.active,a.external,a.periods,a.spread_access);a.output_dir.mkdir(parents=True,exist_ok=False)
    with (a.output_dir/'左侧延时路线试算.json').open('x',encoding='utf-8') as f:json.dump(r,f,ensure_ascii=False,indent=2)
    print(json.dumps({'intervals':[{k:v[k] for k in ['target_total_ps','fixed_delay_ps','directional_delay_ps','right_fold_x_um']} for n,v in r['delays'].items() if n.startswith('loop')],
        'self_crossing':r['topology']['self_crossing_routes'],'crossings':r['crossing_details'],'shape_checks':r['shape_checks']},ensure_ascii=False,indent=2))


if __name__=='__main__':main()
