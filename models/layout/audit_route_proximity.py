"""只做路由几何间距检查；不把间距阈值换算为已验证串扰。

MUX内部不在路由中心线中。其固定端口相距很近，端口外150um过渡
单独报告，不隐去靠近；普通路由仍与对方整条外部路径比较。
"""
import argparse,itertools,json,math
from pathlib import Path
from shapely.geometry import LineString,box
from shapely.ops import substring,unary_union,nearest_points


def analyze(data,gap_um=10.,port_transition_um=150.,minimum_parallel_um=100.):
    routes={n:LineString(p) for n,p in data['routes'].items()}
    bbs=unary_union([box(c['x_um']-22.5,c['y_um']-22.5,c['x_um']+22.5,c['y_um']+22.5)
                    for c in data['crossing_details'] if 'x_um' in c])
    full={n:r.difference(bbs) for n,r in routes.items()}
    ordinary={}
    for n,r in routes.items():
        start=port_transition_um if n.startswith('loop') or n=='output' else 0.
        stop=r.length-port_transition_um if n.startswith('loop') or n=='input' else r.length
        ordinary[n]=substring(r,start,stop).difference(bbs)
    pairs=[]
    for a,b in itertools.combinations(routes,2):
        near_a=ordinary[a].intersection(full[b].buffer(gap_um+.7))
        near_b=ordinary[b].intersection(full[a].buffer(gap_um+.7))
        pa,pb=nearest_points(full[a],full[b])
        pairs.append({'routes':[a,b],'minimum_nominal_edge_gap_um':full[a].distance(full[b])-.7,
            'ordinary_min_edge_gap_um':min(ordinary[a].distance(full[b]),ordinary[b].distance(full[a]))-.7,
            'full_near_length_um':{a:full[a].intersection(full[b].buffer(gap_um+.7)).length,
                                   b:full[b].intersection(full[a].buffer(gap_um+.7)).length},
            'ordinary_near_length_um':{a:near_a.length,b:near_b.length},
            'nearest_points_um':[list(pa.coords)[0],list(pb.coords)[0]],
            'near_boxes_um':[list(g.bounds) for g in [near_a,near_b] if not g.is_empty]})
    segments=[]
    for n,points in data['routes'].items():
        for i,(a,b) in enumerate(zip(points,points[1:])):
            if abs(a[0]-b[0])<1e-7:
                axis='y';coord=a[0];lo,hi=sorted([a[1],b[1]])
            elif abs(a[1]-b[1])<1e-7:
                axis='x';coord=a[1];lo,hi=sorted([a[0],b[0]])
            else:continue
            if hi-lo>=minimum_parallel_um:segments.append((n,i,axis,coord,lo,hi))
    close=[];long_pairs=[]
    for a,b in itertools.combinations(segments,2):
        if a[2]!=b[2]:continue
        overlap=min(a[5],b[5])-max(a[4],b[4])
        if overlap<minimum_parallel_um:continue
        edge_gap=abs(a[3]-b[3])-.7
        row={'routes':[a[0],b[0]],'segments':[a[1],b[1]],'axis':a[2],
             'nominal_edge_gap_um':edge_gap,'parallel_overlap_um':overlap}
        long_pairs.append(row)
        if edge_gap<gap_um-1e-6:close.append(row)
    return {'scope':'0.7um路由中心线几何筛查；黑盒内部及MUX内部不属此几何；不是光学串扰验证',
        'ordinary_edge_gap_target_um':gap_um,'port_transition_report_zone_um':port_transition_um,
        'minimum_parallel_overlap_um':minimum_parallel_um,'pairs':pairs,'close_long_parallel':close,
        'minimum_long_parallel_nominal_edge_gap_um':min((x['nominal_edge_gap_um'] for x in long_pairs),default=None),
        'ordinary_pairwise_clear':all(r['ordinary_min_edge_gap_um']>=gap_um-1e-6 and max(r['ordinary_near_length_um'].values())<.001 for r in pairs),
        'no_close_long_parallel':not close,'optical_crosstalk_verified':False}


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--routes',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    r=analyze(json.loads(a.routes.read_text(encoding='utf-8')))
    with a.output.open('x',encoding='utf-8') as f:json.dump(r,f,ensure_ascii=False,indent=2)
    print(json.dumps({'ordinary_clear':r['ordinary_pairwise_clear'],'close_long_parallel':r['close_long_parallel'],
        'minimum_long_parallel_gap':r['minimum_long_parallel_nominal_edge_gap_um'],
        'near_pairs':[p for p in r['pairs'] if max(p['ordinary_near_length_um'].values())>=.001]},ensure_ascii=False,indent=2))


if __name__=='__main__':main()
