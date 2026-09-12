"""短接续与约254 ps第二回路的布线试排；只写新预算/预览，不写GDS。

保留15 mm电极和750 um MUX核心。未通过交叉/边界检查的候选不能
用于合并GDS；这里的群时延沿用已保存的方向ng和MUX相位斜率。
"""
import argparse
import json
import math
import sys
from pathlib import Path
import numpy as np
from shapely.geometry import LineString, box, Point
from shapely.ops import substring
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from euler_routes import bend_displacement, radius_for_semicircle_displacement, partial_euler_local, append_euler
from route_15mm_left_delays import PathBuilder, build_left
from route_geometry import check_route_network

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
from delay_budget import route_delay


def inputs():
    d = ROOT/'results/layout/_历史归档/迭代版本_20260912'
    old = json.loads((d/'欧拉回路时延预算_v3_无斜直线/欧拉回路_10GHz方向时延预算.json').read_text(encoding='utf-8'))
    merged = json.loads((d/'YSJ合并15mm_功能候选_v10/合并版时延预算.json').read_text(encoding='utf-8'))
    return old, merged


def fixed_delays(active, external):
    old, merged = inputs()
    active_pair = sum(row['fixed_two_active_taper_delay_ps'] for row in old['loops'][:2])/2
    external_pair = old['loops'][0]['fixed_two_external_taper_delay_ps']/2
    return [r['total_ps']-r['passive_route_ps']-active_pair*(1-active/200.)
            -external_pair*(1-external/100.) for r in merged['intervals']]


def trial(active=100., external=50., extension=5., dx=-11620., dy=1620.):
    old, _ = inputs()
    ny = old['passive_rib']['ng_crystal_Y']; nz = old['passive_rib']['ng_crystal_Z']
    fixed = fixed_delays(active, external)
    xl = 2000-active-750-external+dx
    xr = 17000+active+750+external+dx
    q = bend_displacement(math.pi/2, 80.)[0]
    common_center = 29.85+(bend_displacement(math.pi,80.)[1]+3.55)/2
    start_ys = dict(loop1=29.85,loop2=26.6,loop3=-29.85,output=-26.6)
    prefixes = {}
    for name,y in start_ys.items():
        p = PathBuilder([[xr,y+dy]],0.)
        p.straight(x=xr+5)
        radius = radius_for_semicircle_displacement(2*(common_center-y))
        p.bend(math.pi,radius,'右端嵌套欧拉折返')
        prefixes[name] = p
    ends = [dy+33.4,dy-30.15,dy-33.4]
    short = prefixes['loop2']
    short.straight(x=xl-extension)
    short.bend(math.pi/2,80.,'左端转下')
    short.straight(y=ends[1]+q)
    short.bend(math.pi/2,80.,'左端转向MUX')
    short.straight(x=xl)
    routes = {'loop2':short.points}
    bends = {'loop2':short.bends}
    configs = {
        1:dict(down_col=-10845.,entry_y=650.,return_col=-10200.,u_count=1,
               fold_left=None,final_rectangular=True,final_start_y=1050.,
               final_radius=83.,final_x=-10728.638053930058),
        3:dict(down_col=-10890.,entry_y=-400.,return_col=-10530.,u_count=1,
               fold_left=None,final_rectangular=False,final_start_y=None,
               final_radius=80.,final_x=-10671.),
    }
    for i in (1,3):
        name=f'loop{i}'; cfg=configs[i]; fold=-8800. if i==1 else -10100.
        for _ in range(4):
            pts,meta=build_left(prefixes[name].points,i,fold,xl,ends[i-1],cfg)
            error = 300.-fixed[i-1]-route_delay(pts,ny,nz)['directional_delay_ps']
            fold += error*299.792458/ny/2
        routes[name]=pts; bends[name]=meta
    incoming=PathBuilder([[-10380.,-1700.]],0.)
    incoming.straight(x=-8000.)
    incoming.bend(math.pi,80.,'输入底部折返')
    incoming.straight(x=-8400.+q)
    incoming.bend(-math.pi/2,80.,'输入上行')
    incoming.straight(y=1120.-q)
    incoming.bend(math.pi/2,80.,'输入向西')
    incoming.straight(x=-10729.638053930058)
    incoming.bend(-math.pi/2,80.,'输入MUX前上转')
    incoming.straight(y=dy+30.15-q)
    incoming.bend(-math.pi/2,80.,'输入接MUX')
    incoming.straight(x=xl)
    routes['input']=incoming.points
    out=prefixes['output']
    out.straight(x=-10877.+q);out.bend(math.pi/2,80.,'输出向下')
    out.straight(y=100.+q);out.bend(math.pi/2,80.,'输出向东')
    out.straight(x=-9500.)
    out.bend(-math.pi,radius_for_semicircle_displacement(600.),'输出折返')
    out.straight(x=-10380.)
    routes['output']=out.points
    topology=check_route_network(routes)
    intersections=[]
    for record in topology['unexpected_intersections']:
        r=dict(record)
        if 'x_um' in r:
            point=Point(r['x_um'],r['y_um']); usable={}
            for name in r['routes']:
                line=LineString(routes[name]);s=line.project(point)
                segment=substring(line,max(0,s-47.5),min(line.length,s+47.5))
                usable[name]=(s>=47.5 and s+47.5<=line.length and
                              abs(segment.length-math.dist(segment.coords[0],segment.coords[-1]))<1e-5)
            r['straight_95um_each_route']=usable
            r['can_reserve_crossing_L25']=(abs(r['angle_deg']-90)<1e-3 and all(usable.values()))
        intersections.append(r)
    delays={name:route_delay(points,ny,nz) for name,points in routes.items()}
    intervals=[dict(loop=i,total_ps=fixed[i-1]+delays[f'loop{i}']['directional_delay_ps'],fixed_ps=fixed[i-1]) for i in (1,2,3)]
    within={name:box(-10900,-1900,10900,1900).covers(LineString(pts).buffer(8.35)) for name,pts in routes.items()}
    return dict(status='未通过的试排也保留；不可直接作为功能GDS',active_taper_um=active,external_taper_um=external,
                left_extension_um=extension,shift_um=[dx,dy],routes=routes,bends=bends,delays=delays,
                intervals=intervals,intersections=intersections,self_crossing=topology['self_crossing_routes'],
                windows_inside_block=within,full_layout_pass=False,
                note='中心线与工程群时延试算；尚未检查器件体和YSJ冲突；未运行场求解。')


def finish_report(routes, bends, fixed, **params):
    old,_=inputs();ny=old['passive_rib']['ng_crystal_Y'];nz=old['passive_rib']['ng_crystal_Z']
    topology=check_route_network(routes); crossings=[]
    for rec in topology['unexpected_intersections']:
        r=dict(rec)
        if 'x_um' in r:
            point=Point(r['x_um'],r['y_um']);usable={}
            for name in r['routes']:
                line=LineString(routes[name]);s=line.project(point)
                seg=substring(line,max(0,s-47.5),min(line.length,s+47.5))
                usable[name]=s>=47.5 and s+47.5<=line.length and abs(seg.length-math.dist(seg.coords[0],seg.coords[-1]))<1e-5
            r['straight_95um_each_route']=usable
            r['can_reserve_crossing_L25']=abs(r['angle_deg']-90)<1e-3 and all(usable.values())
        crossings.append(r)
    delays={n:route_delay(p,ny,nz) for n,p in routes.items()}
    return dict(status='布线候选，未生成GDS',routes=routes,bends=bends,delays=delays,
        intervals=[dict(loop=i,total_ps=fixed[i-1]+delays[f'loop{i}']['directional_delay_ps'],fixed_ps=fixed[i-1]) for i in (1,2,3)],
        intersections=crossings,self_crossing=topology['self_crossing_routes'],
        windows_inside_block={n:box(-10900,-1900,10900,1900).covers(LineString(p).buffer(8.35)) for n,p in routes.items()},
        full_layout_pass=False,**params)


def add_s(p, length, displacement):
    if abs(displacement)<1e-9:
        p.straight(x=p.points[-1][0]-length);return
    # 当前只用于西行光路：局部正y指向全局南，故转角取负号。
    angle=-2*math.atan2(displacement,length)
    _,info=partial_euler_local(angle,80.)
    radius=80*length/(2*info['displacement_x_um'])
    for _ in range(3):
        _,info=partial_euler_local(angle,radius)
        radius*=length/(2*info['displacement_x_um'])
    if radius<80:raise ValueError('S接续半径小于80um')
    p.bend(angle,radius,'平滑移位1');p.bend(-angle,radius,'平滑移位2')


def full_trial(active=50.,external=50.,extension=160.,dx=-11620.,dy=1620.):
    old,_=inputs();ny=old['passive_rib']['ng_crystal_Y'];nz=old['passive_rib']['ng_crystal_Z'];fixed=fixed_delays(active,external)
    xl=2000-active-750-external+dx;xr=17000+active+750+external+dx
    q=bend_displacement(math.pi/2,80)[0]
    blue_end=dy-30.15;blue_lane=blue_end-2*q-150.
    common_center=(blue_lane+dy+26.6)/2
    starts=dict(loop2=dy+26.6,loop3=dy-29.85,output=dy-26.6)
    prefixes={}
    for name,y in starts.items():
        p=PathBuilder([[xr,y]],0.);p.straight(x=xr+5)
        p.bend(-math.pi,radius_for_semicircle_displacement(2*(y-common_center)),'右侧下翻嵌套')
        prefixes[name]=p
    p=prefixes['loop2'];p.straight(x=xl-extension);p.bend(-math.pi/2,80,'短回路左端上转')
    p.straight(y=blue_end-q);p.bend(-math.pi/2,80,'短回路进入下轨');p.straight(x=xl)
    routes={'loop2':p.points};bends={'loop2':p.bends}
    # 第一回路仍上翻；左侧补足300ps。
    p1=PathBuilder([[xr,dy+29.85]],0.);p1.straight(x=xr+5)
    p1.bend(math.pi,radius_for_semicircle_displacement(bend_displacement(math.pi,80)[1]+3.55),'第一回路上翻')
    cfg=dict(down_col=-10845.,entry_y=650.,return_col=-10200.,u_count=1,fold_left=None,
             final_rectangular=True,final_start_y=1050.,final_radius=83.,final_x=-10728.638053930058)
    fold=-8800.
    for _ in range(4):
        pts,meta=build_left(p1.points,1,fold,xl,dy+33.4,cfg)
        fold+=(300-fixed[0]-route_delay(pts,ny,nz)['directional_delay_ps'])*299.792458/ny/2
    routes['loop1']=pts;bends['loop1']=meta
    green_cross=blue_end-q-50.
    pg=prefixes['loop3'];pg.straight(x=-8000.);add_s(pg,1000.,green_cross-pg.points[-1][1])
    prefix=pg.points[:]
    def green(fold):
        g=PathBuilder(prefix,math.pi);g.straight(x=-10891.+q);g.bend(math.pi/2,80,'第三路左侧下转')
        g.straight(y=-400.+q);g.bend(math.pi/2,80,'第三路进入补偿直段');g.straight(x=fold)
        g.bend(math.pi,80,'第三路补偿折返')
        col=xl-q-5;g.straight(x=col+q);g.bend(-math.pi/2,80,'第三路返回上行')
        g.straight(y=dy-33.4-q);g.bend(-math.pi/2,80,'第三路进入下轨辅助端');g.straight(x=xl)
        return g
    fold=-8800.
    for _ in range(4):
        g=green(fold);fold+=(300-fixed[2]-route_delay(g.points,ny,nz)['directional_delay_ps'])*299.792458/ny/2
    routes['loop3']=g.points;bends['loop3']=pg.bends+g.bends
    out=prefixes['output'];out.straight(x=-9300.)
    output_height=blue_lane+q+52.;add_s(out,600.,output_height-out.points[-1][1])
    out.straight(x=-10000.);out.bend(math.pi/2,80,'输出下转')
    out.straight(y=-500.+q);out.bend(-math.pi/2,80,'输出转向左端面');out.straight(x=-10380.)
    routes['output']=out.points;bends['output']=out.bends
    inc=PathBuilder([[-10380.,-1700.]],0.);inc.straight(x=-7000.);inc.bend(math.pi,80,'输入下部折返')
    inc.straight(x=-7500.+q);inc.bend(-math.pi/2,80,'输入上行');inc.straight(y=1120.-q)
    inc.bend(math.pi/2,80,'输入向西');inc.straight(x=-10729.638053930058)
    inc.bend(-math.pi/2,80,'输入上转');inc.straight(y=dy+30.15-q);inc.bend(-math.pi/2,80,'输入接上轨');inc.straight(x=xl)
    routes['input']=inc.points;bends['input']=inc.bends
    return finish_report(routes,bends,fixed,active_taper_um=active,external_taper_um=external,left_extension_um=extension,
                         shift_um=[dx,dy],note='上下分组重排的比较候选，未获准作为254ps最终版。')


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output-dir',type=Path,required=True)
    ap.add_argument('--extension',type=float,default=5.)
    ap.add_argument('--active',type=float,default=100.)
    ap.add_argument('--external',type=float,default=50.)
    ap.add_argument('--full-network',action='store_true')
    a=ap.parse_args();a.output_dir.mkdir(parents=True,exist_ok=False)
    r=(full_trial if a.full_network else trial)(a.active,a.external,a.extension)
    with (a.output_dir/'紧凑回路试排.json').open('x',encoding='utf-8') as f:json.dump(r,f,ensure_ascii=False,indent=2)
    plt.rcParams.update({'font.sans-serif':['Microsoft YaHei'],'axes.unicode_minus':False})
    fig,ax=plt.subplots(figsize=(13,8),layout='constrained')
    for name,pts in r['routes'].items():
        axy=np.asarray(pts);ax.plot(axy[:,0],axy[:,1],label=name,lw=1)
    for i,c in enumerate(r['intersections']):
        if 'x_um' in c:
            ax.scatter(c['x_um'],c['y_um'],color='green' if c['can_reserve_crossing_L25'] else 'red',s=25)
            ax.text(c['x_um']+15,c['y_um'],str(i),fontsize=8)
    ax.set(xlim=(-10920,-9000),ylim=(-700,1920),aspect='equal',xlabel='x (µm)',ylabel='y (µm)',title='短回路试排：红点不能直接放置正交交叉器（不是最终版图）')
    ax.legend();ax.grid(alpha=.2);fig.savefig(a.output_dir/'左侧布线冲突诊断.png',dpi=150);plt.close(fig)
    print(json.dumps({k:r[k] for k in ['intervals','self_crossing','windows_inside_block','intersections']},ensure_ascii=False,indent=2))


if __name__=='__main__':main()
