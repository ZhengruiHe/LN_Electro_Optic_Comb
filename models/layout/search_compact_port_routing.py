"""将固定器件视作端口模块，搜索上下翻转/出线分组；不运行场求解。

镜像只作为平面端口规划候选，不能自动当作已确认版图。不会改源GDS。
保留15mm电极、750um MUX耦合轮廓、100/50um外接续假设。
"""
import argparse
import itertools
import json
import math
from pathlib import Path
from functools import lru_cache
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from shapely.geometry import LineString
from route_compact_254_trial import inputs, fixed_delays, finish_report
from route_15mm_left_delays import PathBuilder, build_left
from euler_routes import bend_displacement, radius_for_semicircle_displacement
from delay_budget import route_delay

radius_for_displacement=lru_cache(maxsize=256)(radius_for_semicircle_displacement)
Q=bend_displacement(math.pi/2,80.)[0]
D=bend_displacement(math.pi,80.)[1]


def build(flips,split,active=100.,external=50.,dx=-11620.,dy=1620.):
    lu,ru,ll,rl=flips
    lx=lambda main:dy+30+( (-.15 if main else -3.4) if lu else (.15 if main else 3.4))
    ly=lambda main:dy-30+( (.15 if main else 3.4) if ll else (-.15 if main else -3.4))
    starts={'loop1':dy+30+(.15 if ru else -.15),'loop2':dy+30+(3.4 if ru else -3.4),
            'loop3':dy-30+(-.15 if rl else .15),'output':dy-30+(-3.4 if rl else 3.4)}
    ends={'loop1':lx(False),'loop2':ly(True),'loop3':ly(False),'input':lx(True)}
    order=sorted(starts,key=starts.get,reverse=True)
    signs={n:1 if i<split else -1 for i,n in enumerate(order)}
    centers={}
    for sign in (1,-1):
        values=[starts[n] for n in order if signs[n]==sign]
        if values:centers[sign]=(max(values)+(D+3.55)/2 if sign==1 else min(values)-(D+3.55)/2)
    xl=2000-active-750-external+dx;xr=17000+active+750+external+dx
    prefixes={}
    for n,y in starts.items():
        p=PathBuilder([[xr,y]],0.);p.straight(x=xr+5)
        p.bend(signs[n]*math.pi,radius_for_displacement(2*abs(centers[signs[n]]-y)),'右端按端口排序折返')
        prefixes[n]=p
    old,_=inputs();ny=old['passive_rib']['ng_crystal_Y'];nz=old['passive_rib']['ng_crystal_Z'];fixed=fixed_delays(active,external)
    def compact(extension):
        p=PathBuilder(prefixes['loop2'].points,math.pi)
        p.straight(x=xl-extension)
        delta=ends['loop2']-p.points[-1][1]
        if abs(delta)>=2*Q:
            sign=1 if delta<0 else -1
            p.bend(sign*math.pi/2,80.,'左端短回路转向1')
            p.straight(y=ends['loop2']+(Q if delta<0 else -Q))
            p.bend(sign*math.pi/2,80.,'左端短回路转向2')
        elif abs(delta)>=D:
            p.bend((-1 if delta>0 else 1)*math.pi,radius_for_displacement(abs(delta)),'左端短回路180度')
        else:raise ValueError('短回路不能容纳R80折返')
        p.straight(x=xl)
        return p
    short=compact(5.)
    base=fixed[1]+route_delay(short.points,ny,nz)['directional_delay_ps']
    if base>254.5:raise ValueError('该端口分组的紧凑回路超过254.5ps')
    extension=5.+max(0,254.-base)*299.792458/ny/2
    short=compact(extension)
    routes={'loop2':short.points};bends={'loop2':prefixes['loop2'].bends+short.bends}
    configs={1:dict(down_col=-10845.,entry_y=650.,return_col=-10200.,u_count=1,fold_left=None,
                    final_rectangular=True,final_start_y=1050.,final_radius=83.,final_x=-10728.638053930058),
             3:dict(down_col=-10890.,entry_y=-400.,return_col=-10530.,u_count=1,fold_left=None,
                    final_rectangular=False,final_start_y=None,final_radius=80.,final_x=-10671.)}
    for i in (1,3):
        name=f'loop{i}';fold=-8800. if i==1 else -10100.
        for _ in range(4):
            pts,b=build_left(prefixes[name].points,i,fold,xl,ends[name],configs[i])
            fold+=(300-fixed[i-1]-route_delay(pts,ny,nz)['directional_delay_ps'])*299.792458/ny/2
        routes[name]=pts;bends[name]=prefixes[name].bends+b
    out=PathBuilder(prefixes['output'].points,math.pi);out.straight(x=-10877.+Q);out.bend(math.pi/2,80.,'输出向下')
    out.straight(y=100.+Q);out.bend(math.pi/2,80.,'输出向东');out.straight(x=-9500.)
    out.bend(-math.pi,radius_for_displacement(600.),'输出折返');out.straight(x=-10380.)
    routes['output']=out.points;bends['output']=prefixes['output'].bends+out.bends
    inc=PathBuilder([[-10380.,-1700.]],0.);inc.straight(x=-7000.);inc.bend(math.pi,80.,'输入折返')
    inc.straight(x=-7500.+Q);inc.bend(-math.pi/2,80.,'输入上行');inc.straight(y=1120.-Q)
    inc.bend(math.pi/2,80.,'输入向西');inc.straight(x=-10729.638053930058)
    inc.bend(-math.pi/2,80.,'输入接续1');inc.straight(y=ends['input']-Q);inc.bend(-math.pi/2,80.,'输入接续2');inc.straight(x=xl)
    routes['input']=inc.points;bends['input']=inc.bends
    result=finish_report(routes,bends,fixed,mux_mirror_flags=dict(zip(['left_upper','right_upper','left_lower','right_lower'],flips)),
        right_escape_signs=signs,right_port_order=order,left_port_y=ends,shift_um=[dx,dy],
        active_taper_um=active,external_taper_um=external,extension_um=extension)
    return result


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output-dir',type=Path,required=True);a=ap.parse_args()
    a.output_dir.mkdir(parents=True,exist_ok=False)
    rows=[];best=None
    for flips in itertools.product([False,True],repeat=4):
        for split in range(5):
            try:r=build(flips,split)
            except ValueError as ex:
                rows.append(dict(flips=flips,split=split,rejected=str(ex)));continue
            bad=[x for x in r['intersections'] if not x.get('can_reserve_crossing_L25',False)]
            score=(sum(not x for x in r['windows_inside_block'].values()),len(r['self_crossing']),len(bad),len(r['intersections']),sum(flips))
            rows.append(dict(flips=flips,split=split,score=score,delay_ps=r['intervals'][1]['total_ps']))
            if best is None or score<best[0]:
                best=(score,r);print(json.dumps({'best':score,'flips':flips,'split':split},ensure_ascii=False),flush=True)
    with (a.output_dir/'端口分组试排汇总.json').open('x',encoding='utf-8') as f:json.dump(rows,f,ensure_ascii=False,indent=2)
    if best:
        r=best[1]
        with (a.output_dir/'候选路线.json').open('x',encoding='utf-8') as f:json.dump(r,f,ensure_ascii=False,indent=2)
        plt.rcParams.update({'font.sans-serif':['Microsoft YaHei'],'axes.unicode_minus':False})
        fig,ax=plt.subplots(figsize=(12,9),layout='constrained')
        for n,pts in r['routes'].items():
            p=np.asarray(pts);ax.plot(p[:,0],p[:,1],label=n,lw=1)
        for i,c in enumerate(r['intersections']):
            if 'x_um' in c:ax.scatter(c['x_um'],c['y_um'],c='green' if c['can_reserve_crossing_L25'] else 'red',s=25);ax.text(c['x_um'],c['y_um'],str(i))
        ax.set(xlim=(-10920,-9000),ylim=(-700,1920),aspect='equal',title='固定器件端口分组候选：仍须完整掩膜检查',xlabel='x (µm)',ylabel='y (µm)');ax.legend();ax.grid(alpha=.2)
        fig.savefig(a.output_dir/'左侧端口分组.png',dpi=150);plt.close(fig)
        print(json.dumps({'best_score':best[0],'intervals':r['intervals'],'bad':[c for c in r['intersections'] if not c.get('can_reserve_crossing_L25')]},ensure_ascii=False,indent=2))


if __name__=='__main__':main()
