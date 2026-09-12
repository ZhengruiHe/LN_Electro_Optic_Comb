"""横向S形接入：冻结MUX/电极，以独立周期预算联合重排五条光路。"""
import argparse,json,math
from pathlib import Path
from shapely.geometry import LineString,Point,box
from shapely.ops import substring
from route_15mm_left_delays import PathBuilder
from route_five_ports_global import s_offset
from route_compact_254_trial import inputs,fixed_delays
from euler_routes import bend_displacement,radius_for_semicircle_displacement,partial_euler_local
from route_geometry import check_route_network
from search_fixed_block_placement import geometry
from delay_budget import route_delay

ROOT=Path(__file__).resolve().parents[2]


def build(periods=(3.,3.5,4.),access_radius1=165.,right_fanout_um=80.):
    if periods[0]%1 or periods[2]%1 or abs(periods[1]%1-.5)>1e-9:
        raise ValueError('各段独立取整数/半整数/整数周期')
    old,_=inputs();ny=old['passive_rib']['ng_crystal_Y'];nz=old['passive_rib']['ng_crystal_Z']
    fixed=fixed_delays(100.,50.);dx,dy=-11620.,1620.;xl,xr=-10520.,6280.
    q=bend_displacement(math.pi/2,80.)[0];d=bend_displacement(math.pi,80.)[1]
    first_lane=dy+29.85+d
    prefixes={}
    for name,y,lane in [('loop1',dy+29.85,first_lane),('loop2',dy+26.6,first_lane+20.),
                       ('output',dy-26.6,first_lane+40.),('loop3',dy-29.85,first_lane+60.)]:
        p=PathBuilder([[xr,y]],0.)
        if name in ('output','loop3') and right_fanout_um:
            s_offset(p,right_fanout_um,5.875 if name=='output' else -5.875)
        p.straight(x=p.points[-1][0]+5)
        p.bend(math.pi,radius_for_semicircle_displacement(lane-p.points[-1][1]),'20um间距返程折返')
        prefixes[name]=p
    local,_=partial_euler_local(math.pi,access_radius1)
    min_x=xl-10.-float(local[:,0].max())
    first_down=-math.ceil((-min_x+16.)/5.)*5.
    configs=[dict(down_col=first_down,entry_y=650.,return_col=-8900.,return_y=1050.,u_count=1,final_radius=access_radius1),
             dict(down_col=first_down-20.,entry_y=300.,return_col=-8000.,return_y=1180.,u_count=3,final_radius=90.),
             dict(down_col=first_down-60.,entry_y=-400.,return_col=-9200.,return_y=dy-33.4-d,u_count=1,final_radius=80.)]
    ends=[dy+33.4,dy-30.15,dy-33.4]
    routes={};delays={}
    def loop(i,fold):
        cfg=configs[i];p=PathBuilder(prefixes[f'loop{i+1}'].points,math.pi)
        p.straight(x=cfg['down_col']+q);p.bend(math.pi/2,80.,'左侧转下')
        p.straight(y=cfg['entry_y']+q);p.bend(math.pi/2,80.,'进入横向补偿段')
        p.straight(x=fold);p.bend(math.pi,80.,'补偿折返1')
        if cfg['u_count']==3:
            p.straight(x=-7000.);p.bend(-math.pi,80.,'补偿折返2')
            p.straight(x=fold);p.bend(math.pi,80.,'补偿折返3')
        p.straight(x=cfg['return_col']+q);p.bend(-math.pi/2,80.,'进入独立上行通道')
        p.straight(y=cfg['return_y']-q);p.bend(math.pi/2,80.,'转入横向接入段')
        if i<2:
            target_y=ends[i]-bend_displacement(math.pi,cfg['final_radius'])[1]
            p.straight(x=-9600.);s_offset(p,930.,-(target_y-p.points[-1][1]))
        else:p.straight(x=xl-10.)
        p.bend(-math.pi,cfg['final_radius'],'横向接入后的欧拉折返')
        p.straight(x=xl)
        return p
    for i in range(3):
        fold=[-7000.,-6000.,-2000.][i];history=[]
        for _ in range(4):
            p=loop(i,fold);r=route_delay(p.points,ny,nz);error=periods[i]*100-fixed[i]-r['directional_delay_ps']
            history.append({'fold_x_um':fold,'error_ps':error})
            fold+=error*299.792458/ny/(configs[i]['u_count']+1)
        p=loop(i,fold);name=f'loop{i+1}';routes[name]=p.points
        delays[name]={**route_delay(p.points,ny,nz),'bends':prefixes[name].bends+p.bends,
            'config':configs[i],'right_fold_x_um':fold,'target_total_ps':100*periods[i],
            'fixed_delay_ps':fixed[i],'target_passive_delay_ps':100*periods[i]-fixed[i],'history':history}
    inc=PathBuilder([[-10380.,-1700.]],0.);inc.straight(x=-10280.);inc.straight(x=-7650.)
    inc.bend(math.pi,80.,'输入下部折返');inc.straight(x=-8150.+q);inc.bend(-math.pi/2,80.,'输入独立上行')
    inc.straight(y=1120.-q);inc.bend(math.pi/2,80.,'输入横向接入')
    inc.straight(x=-9600.);s_offset(inc,930.,-((dy+30.15-bend_displacement(math.pi,130.)[1])-1120.))
    inc.bend(-math.pi,130.,'输入接入MUX前折返');inc.straight(x=xl)
    routes['input']=inc.points;delays['input']={**route_delay(inc.points,ny,nz),'bends':inc.bends}
    out=PathBuilder(prefixes['output'].points,math.pi);out.straight(x=first_down-40.+q)
    out.bend(math.pi/2,80.,'输出左侧下行');out.straight(y=100.+q)
    out.bend(math.pi/2,80.,'输出横向避让');out.straight(x=-9600.-q)
    out.bend(-math.pi/2,80.,'输出正交下行');out.straight(y=-650.+q)
    out.bend(-math.pi/2,80.,'输出接左端面')
    out.straight(x=-10280.);out.straight(x=-10380.)
    routes['output']=out.points;delays['output']={**route_delay(out.points,ny,nz),'bends':prefixes['output'].bends+out.bends}
    topology=check_route_network(routes);crossings=[]
    for c in topology['unexpected_intersections']:
        r=dict(c)
        if 'x_um' in r:
            point=Point(r['x_um'],r['y_um']);r['straight_245um_each_route']={}
            for name in r['routes']:
                line=LineString(routes[name]);s=line.project(point);seg=substring(line,max(0,s-122.5),min(line.length,s+122.5))
                r['straight_245um_each_route'][name]=s>=122.5 and s+122.5<=line.length and abs(seg.length-math.dist(seg.coords[0],seg.coords[-1]))<1e-5
        crossings.append(r)
    source=json.loads((ROOT/'results/layout/_历史归档/迭代版本_20260912/合并YSJ_版图确认_20260909_v1/YSJ源版图几何.json').read_text(encoding='utf-8'))
    oc,os,oe=(geometry(source,k) for k in [(20,0),(10,0),(21,2)])
    checks={}
    for name,pts in routes.items():
        line=LineString(pts);c=line.buffer(.35,cap_style=2);w=line.buffer(8.6,cap_style=2)
        checks[name]={'YSJ_core_overlap_um2':c.intersection(oc).area,'SiN_removal_over_YSJ_SiN_um2':w.intersection(os).area,
            'LN1_over_YSJ_LN2_etch_um2':c.intersection(oe).area,'window_inside_block':box(-10900,-1900,10900,1900).covers(w)}
    return dict(status='横向接入布线预算；未写GDS',electrode_length_um=15000,shift_um=[dx,dy],active_taper_um=100.,external_taper_um=50.,
        loop_periods=list(periods),layout_variant='horizontal_access',expected_crossing_count=len(crossings),
        access_radius1_um=access_radius1,right_fanout_um=right_fanout_um,
        edge_coupler_internal_ports_um=[[-10380.,-1700.],[-10380.,-650.]],
        routes=routes,delays=delays,topology=topology,crossing_details=crossings,shape_checks=checks,
        additional_cleanup_rois=[[6280.,1500.,6880.,1900.]],
        window_merge_caps=[
            {'layer':[20,1],'bounds_um':[6372.,1658.,6385.,1670.],'reason':'右端相邻窗口收尖区域改为宽端圆角合并'},
            {'layer':[20,1],'bounds_um':[6306.,1587.,6320.,1597.],'reason':'右端S扇出窗口收尖合并，不改变核心'},
            {'layer':[10,2],'bounds_um':[6418.,1725.,6429.,1737.],'reason':'右端SiN挖除窗口小凹口合并'},
            {'layer':[10,2],'bounds_um':[-10674.,1480.,-10660.,1494.],'reason':'左端嵌套弯SiN窗口收尖合并'},
            {'layer':[10,2],'bounds_um':[6464.,1596.,6478.,1607.],'reason':'右端SiN窗口收尖合并'},
            {'layer':[20,1],'bounds_um':[-10636.,1563.,-10630.,1569.],'reason':'固定7.2um窗口在左下相邻欧拉弯汇合处封口'},
            {'layer':[10,2],'bounds_um':[6484.,1604.,6491.,1609.],'reason':'固定17.2um窗口在右下S扇出汇合处封口'},
        ],input_and_final_output_included=True,normalize_mux_windows_to_pdk=True,full_layout_pass=False)


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output-dir',type=Path,required=True)
    ap.add_argument('--periods',nargs=3,type=float,default=[3.,3.5,4.]);a=ap.parse_args()
    r=build(a.periods);a.output_dir.mkdir(parents=True,exist_ok=False)
    with (a.output_dir/'左侧延时路线试算.json').open('x',encoding='utf-8') as f:json.dump(r,f,ensure_ascii=False,indent=2)
    print(json.dumps({'intervals':[{k:v[k] for k in ['target_total_ps','fixed_delay_ps','directional_delay_ps','right_fold_x_um']} for n,v in r['delays'].items() if n.startswith('loop')],
        'self_crossing':r['topology']['self_crossing_routes'],'crossings':r['crossing_details'],'shape_checks':r['shape_checks']},ensure_ascii=False,indent=2))


if __name__=='__main__':main()
