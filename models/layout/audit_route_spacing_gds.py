"""实际GDS路由芯层的10um间距检查；固定MUX端口过渡区单独记录。"""
import json,math,itertools
import pya
with open(payload_file,encoding='utf-8') as f:p=json.load(f)
with open(p['routes_file'],encoding='utf-8') as f:routes=json.load(f)['routes']
l=pya.Layout();l.read(p['merged_gds']);ours=l.cell(p['ours_top'])
actual=pya.Region(ours.begin_shapes_rec(l.layer(20,0))).merged()
with open(p['build_payload'],encoding='utf-8') as f:bp=json.load(f)
bb=pya.Region()
for x,y in bp['crossing_centers_um']:
    bb+=pya.Region(pya.DBox(x-22.5,y-22.5,x+22.5,y+22.5).to_itype(.001))

def length(points):return sum(math.dist(a,b) for a,b in zip(points,points[1:]))
def trim(points,start,end):
    out=[];s=0.
    for a,b in zip(points,points[1:]):
        ds=math.dist(a,b)
        if ds<1e-12:continue
        lo=max(start,s);hi=min(end,s+ds)
        if hi>=lo and hi>s and lo<s+ds:
            for t in [(lo-s)/ds,(hi-s)/ds]:
                q=[a[0]+t*(b[0]-a[0]),a[1]+t*(b[1]-a[1])]
                if not out or math.dist(out[-1],q)>1e-9:out.append(q)
        s+=ds
        if s>end:break
    return out
def corridor(points):
    return pya.Region(pya.DPath([pya.DPoint(*q) for q in points],2.).to_itype(.001))
def region(points):
    return ((actual&corridor(points))-bb).merged()
full={};ordinary={}
for name,points in routes.items():
    full[name]=region(points)
    start=150. if name.startswith('loop') or name=='output' else 0.
    end=length(points)-(150. if name.startswith('loop') or name=='input' else 0.)
    ordinary[name]=region(trim(points,start,end))
self_counts={name:reg.space_check(10000).size() for name,reg in ordinary.items()}
pairs=[]
for a,b in itertools.combinations(routes,2):
    ab=ordinary[a].separation_check(full[b],10000)
    ba=ordinary[b].separation_check(full[a],10000)
    pairs.append({'routes':[a,b],'ordinary_A_to_full_B':ab.size(),'ordinary_B_to_full_A':ba.size(),
        'boxes_um':[[e.bbox().left*.001,e.bbox().bottom*.001,e.bbox().right*.001,e.bbox().top*.001]
                    for coll in [ab,ba] for e in coll.each()]})
# 不仅查五条新路由彼此，也查它们与冻结MUX/有源芯层及原YSJ芯层的距离。
all_corridors=pya.Region()
for points in routes.values():all_corridors+=corridor(points)
frozen_core=(actual-all_corridors-bb).merged()
ysj_core=pya.Region(l.cell(p['ysj_top']).begin_shapes_rec(l.layer(20,0))).merged()
global_pairs=[]
for name in routes:
    for label,reg,other in [('ordinary_to_frozen_core',ordinary[name],frozen_core),
                           ('full_to_YSJ_core',full[name],ysj_core)]:
        markers=reg.separation_check(other,10000)
        global_pairs.append({'route':name,'check':label,'count':markers.size(),
            'boxes_um':[[e.bbox().left*.001,e.bbox().bottom*.001,e.bbox().right*.001,e.bbox().top*.001]
                        for e in markers.each()]})
result={'scope':'实际GDS芯层，含真实1.2um接口宽度；不是光学串扰求解',
        'edge_gap_target_um':10.,'port_transition_report_zone_um':150.,
        'port_zone_note':'每个固定MUX外端口沿光路150um单独报告，不代表该区域无耦合；黑盒内部按PDK占位排除',
        'self_space_markers':self_counts,'pairwise':pairs,'frozen_and_YSJ_core_checks':global_pairs,
        'ordinary_route_spacing_pass':not any(self_counts.values()) and all(x['ordinary_A_to_full_B']==0 and x['ordinary_B_to_full_A']==0 for x in pairs)
            and all(x['count']==0 for x in global_pairs)}
with open(globals().get('report_path',p['spacing_report']),'x',encoding='utf-8') as f:json.dump(result,f,ensure_ascii=False,indent=2)
print(json.dumps(result,ensure_ascii=False,indent=2))
if not result['ordinary_route_spacing_pass']:raise RuntimeError('普通路由实际芯层间距未达到10um')
