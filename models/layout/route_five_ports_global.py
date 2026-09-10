"""五路联合重排：保留器件轮廓，先求路由和时延，再允许生成GDS。

输出只写新目录。A/E是器件外部接续长度，不压缩750um耦合轮廓。
右下MUX轴向镜像仅为摆位候选，不改变晶轴为90度；须回读确认。
"""
import argparse
import json
import math
import sys
from pathlib import Path
from functools import lru_cache
import numpy as np
from scipy.optimize import brentq
from shapely.geometry import LineString, box
from shapely.affinity import translate, scale
from shapely.ops import unary_union
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from route_15mm_left_delays import PathBuilder
from route_compact_254_trial import inputs, fixed_delays, finish_report
from euler_routes import bend_displacement, radius_for_semicircle_displacement, partial_euler_local
from search_fixed_block_placement import geometry, draw_layer
from delay_budget import route_delay

ROOT=Path(__file__).resolve().parents[2]
R=80.
Q=bend_displacement(math.pi/2,R)[0]
D=bend_displacement(math.pi,R)[1]
radius_for=lru_cache(maxsize=128)(radius_for_semicircle_displacement)


def s_offset(p,length,offset):
    """沿当前切向的前进长度及左法向偏移。"""
    angle=2*math.atan2(offset,length)
    _,info=partial_euler_local(angle,R)
    radius=R*length/(2*info['displacement_x_um'])
    for _ in range(3):
        _,info=partial_euler_local(angle,radius)
        radius*=length/(2*info['displacement_x_um'])
    if radius<R:raise ValueError('S形连接半径不足80um')
    p.bend(angle,radius,'路由S渐移1');p.bend(-angle,radius,'路由S渐移2')


def generate(active=20.,external=20.,dx=-11620.,dy=-1600.,target=254.,output_y=-1840.):
    old,_=inputs();ny=old['passive_rib']['ng_crystal_Y'];nz=old['passive_rib']['ng_crystal_Z']
    fixed=fixed_delays(active,external)
    xl=2000-active-750-external+dx;xr=17000+active+750+external+dx
    # 从两只左端口逆向构造并排扇出，之后反转即为实际入射光路。
    center=dy+33.4+D/2
    fan={}
    for name,y,lead in [('input',dy+30.15,0.),('loop1',dy+33.4,60.)]:
        p=PathBuilder([[xl,y]],math.pi);p.straight(x=xl-5)
        p.bend(-math.pi,radius_for(2*(center-y)),'左端口嵌套折返')
        p.straight(x=xl-5+lead)
        p.bend(math.pi/2,R,'左端口正交扇出')
        fan[name]=p
    blue_y=max(p.points[-1][1] for p in fan.values())+50.
    for p in fan.values():p.straight(y=blue_y+70.)
    def blue(extension):
        p=PathBuilder([[xr,dy+26.6]],0.);p.straight(x=xr+5)
        p.bend(math.pi,radius_for(blue_y-dy-26.6),'第二回路右端折返')
        p.straight(x=xl-extension);p.bend(math.pi/2,R,'第二回路左端向下')
        p.straight(y=dy-30.15+Q);p.bend(math.pi/2,R,'第二回路接下轨');p.straight(x=xl)
        return p
    bp=blue(50.)
    base=fixed[1]+route_delay(bp.points,ny,nz)['directional_delay_ps']
    extension=50.+(target-base)*299.792458/ny/2
    if extension<35.:
        extension=35.  # 留出左端扇出与短回路的几何间隔，记录实际预算而非强报目标。
    bp=blue(extension)
    routes={'loop2':bp.points};bends={'loop2':bp.bends}
    # 第一回路的小U位于大U内部；在左侧空区正交越过第二回路后补时延。
    p=PathBuilder([[xr,dy+29.85]],0.);p.straight(x=xr+5);p.bend(math.pi,R,'第一回路右端小折返')
    cross_x=-9000.
    p.straight(x=cross_x+Q);p.bend(-math.pi/2,R,'第一回路正交出线')
    p.straight(y=blue_y+60.);p.bend(math.pi/2,R,'第一回路进入左侧补偿区')
    prefix=p.points[:];prefix_bends=p.bends[:]
    fc=fan['loop1'].points[-1][0]
    def one(height):
        q=PathBuilder(prefix,math.pi);q.straight(x=-10000.+Q);q.bend(-math.pi/2,R,'第一回路上行')
        q.straight(y=height-Q);q.bend(math.pi/2,R,'第一回路向西')
        q.straight(x=-10680.);q.bend(-math.pi,R,'第一回路补偿折返')
        q.straight(x=fc-Q);q.bend(-math.pi/2,R,'第一回路返回下行')
        q.straight(y=fan['loop1'].points[-1][1])
        q.points.extend(fan['loop1'].points[-2::-1])
        return q
    hlow=max(blue_y+60+3*Q,blue_y+70-D+Q)+1
    def err1(h):return fixed[0]+route_delay(one(h).points,ny,nz)['directional_delay_ps']-300.
    h1=brentq(err1,hlow,1800.)
    p1=one(h1);routes['loop1']=p1.points;bends['loop1']=prefix_bends+p1.bends+fan['loop1'].bends
    # 右下MUX采用沿传播轴镜像的实例摆位；辅助输出在下侧直接去右端面。
    right_center=(dy+26.6+blue_y)/2
    pg=PathBuilder([[xr,dy-30.15]],0.);pg.straight(x=xr+5)
    pg.bend(math.pi,radius_for(2*(right_center-(dy-30.15))),'第三回路右端外层折返')
    pg.straight(x=-7500.);s_offset(pg,1000.,-(blue_y+240.-pg.points[-1][1]))
    green_prefix=pg.points[:]
    def three(height):
        q=PathBuilder(green_prefix,math.pi);q.straight(x=-9600.+Q);q.bend(-math.pi/2,R,'第三回路进入补偿区')
        q.straight(y=height-Q);q.bend(math.pi/2,R,'第三回路向西')
        q.straight(x=-10860.+Q);q.bend(math.pi/2,R,'第三回路左侧下行')
        q.straight(y=dy-33.4+Q);q.bend(math.pi/2,R,'第三回路从下侧接入辅助端');q.straight(x=xl)
        return q
    hlo=blue_y+240+2*Q+1
    def err3(h):return fixed[2]+route_delay(three(h).points,ny,nz)['directional_delay_ps']-300.
    hg=brentq(err3,hlo,1870.)
    p3=three(hg);routes['loop3']=p3.points;bends['loop3']=pg.bends+p3.bends
    # 输入与补偿回路分区，左端面位置避开模组核心；不参与三个周期预算。
    pi=PathBuilder([[-10380.,-1720.]],0.);pi.straight(x=-10280.)
    pi.bend(-math.pi,R,'输入端面后下折返');pi.straight(x=-10600.+Q);pi.bend(-math.pi/2,R,'输入左侧上行')
    hi=h1-70.
    pi.straight(y=hi-Q);pi.bend(-math.pi/2,R,'输入向东')
    fi=fan['input'].points[-1][0]
    pi.straight(x=fi-Q);pi.bend(-math.pi/2,R,'输入接扇出段')
    pi.straight(y=fan['input'].points[-1][1]);pi.points.extend(fan['input'].points[-2::-1])
    routes['input']=pi.points;bends['input']=pi.bends+fan['input'].bends
    out=PathBuilder([[xr,dy-33.4]],0.);out.straight(x=xr+100.)
    s_offset(out,1000.,output_y-out.points[-1][1]);out.straight(x=10380.)
    routes['output']=out.points;bends['output']=out.bends
    r=finish_report(routes,bends,fixed,active_taper_um=active,external_taper_um=external,
        shift_um=[dx,dy],mux_mirror_flags={'right_lower':True},blue_lane_y_um=blue_y,
        extension_um=extension,loop1_height_um=h1,loop3_height_um=hg,
        edge_ports_um=[[-10380.,-1720.],[10380.,output_y]],
        scope='仅端口摆位、外部接续与路由；未写GDS；镜像摆位/20um接续不是新场求解结果')
    # 原YSJ的跨层刻蚀安全也必须检查，不只检查同层线交点。
    source=json.loads((ROOT/'results/layout/合并YSJ_版图确认_20260909_v1/YSJ源版图几何.json').read_text(encoding='utf-8'))
    r['YSJ_route_checks']={}
    for name,pts in routes.items():
        line=LineString(pts);core=line.buffer(.35);trench=line.buffer(8.35)
        r['YSJ_route_checks'][name]={
            'LN1_overlap_um2':core.intersection(geometry(source,(20,0))).area,
            'SiN_removal_over_SiN_core_um2':trench.intersection(geometry(source,(10,0))).area,
            'LN1_over_old_LN2_etch_um2':core.intersection(geometry(source,(21,2))).area}
    return r


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output-dir',type=Path,required=True)
    ap.add_argument('--active',type=float,default=20.);ap.add_argument('--external',type=float,default=20.)
    ap.add_argument('--dx',type=float,default=-11620.);ap.add_argument('--dy',type=float,default=-1600.)
    a=ap.parse_args();a.output_dir.mkdir(parents=True,exist_ok=False)
    r=generate(a.active,a.external,a.dx,a.dy)
    with (a.output_dir/'五路联合布线.json').open('x',encoding='utf-8') as f:json.dump(r,f,ensure_ascii=False,indent=2)
    plt.rcParams.update({'font.sans-serif':['Microsoft YaHei'],'axes.unicode_minus':False})
    fig,axes=plt.subplots(1,2,figsize=(17,8),layout='constrained')
    source=json.loads((ROOT/'results/layout/合并YSJ_版图确认_20260909_v1/YSJ源版图几何.json').read_text(encoding='utf-8'))
    for ax in axes:
        draw_layer(ax,source,(20,0),0,0,'#bbbbbb')
        for name,pts in r['routes'].items():
            p=np.asarray(pts);ax.plot(p[:,0]/1000,p[:,1]/1000,label=name,lw=.9)
        for i,c in enumerate(r['intersections']):
            if 'x_um' in c:ax.scatter(c['x_um']/1000,c['y_um']/1000,c='green' if c['can_reserve_crossing_L25'] else 'red',s=20);ax.text(c['x_um']/1000,c['y_um']/1000,str(i),fontsize=8)
        ax.grid(alpha=.2);ax.set(aspect='equal',xlabel='x (mm)',ylabel='y (mm)')
    axes[0].set(xlim=(-10.92,-8.5),ylim=(-1.92,1.92),title='左侧联合布线（须回读GDS）');axes[0].legend()
    axes[1].set(xlim=(-11.1,11.1),ylim=(-2,2),title='全block位置：灰色为YSJ LN1')
    fig.savefig(a.output_dir/'五路联合布线预览.png',dpi=160);plt.close(fig)
    keys=['intervals','self_crossing','windows_inside_block','YSJ_route_checks','loop1_height_um','loop3_height_um','extension_um']
    print(json.dumps({k:r[k] for k in keys},ensure_ascii=False,indent=2))
    print(json.dumps({'crossings':len(r['intersections']),'bad':[c for c in r['intersections'] if not c.get('can_reserve_crossing_L25')]},ensure_ascii=False,indent=2))


if __name__=='__main__':main()
