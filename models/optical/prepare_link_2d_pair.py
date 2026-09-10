"""从实际布线路径截取成对接入弯，保留相对位置、晶向和近邻距离。"""
import argparse,json,math
from pathlib import Path
from shapely.geometry import LineString,box
from shapely.ops import substring,unary_union
from shapely.affinity import translate


def polygons(g):
    pieces=list(g.geoms) if hasattr(g,'geoms') else [g]
    if any(p.interiors for p in pieces):raise ValueError('有孔几何需要拆分，不能直接填成实体')
    return [list(p.exterior.coords) for p in pieces]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--routes',type=Path,required=True)
    p.add_argument('--pair',choices=['left_upper','left_lower','right_upper','right_lower'],required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();r=json.loads(a.routes.read_text(encoding='utf-8'))
    names={'left_upper':['loop1','input'],'left_lower':['loop2','loop3'],
           'right_upper':['loop1','loop2'],'right_lower':['output','loop3']}[a.pair]
    left=a.pair.startswith('left');plane=-10526. if left else 6282.
    lines=[];near=[];far=[]
    for name in names:
        line=LineString(r['routes'][name]);bends=r['delays'][name]['bends']
        if left:
            last=bends[-1]
            if abs(last['angle_deg'])!=180:raise ValueError('末段不是预期欧拉折返')
            segment=substring(line,line.length-last['centerline_length_um']-28,line.length)
            ny=line.coords[-1][1]
        else:
            distance=5.
            for b in bends:
                distance+=b['centerline_length_um']
                if b['label']=='20um间距返程折返':break
            # 右下先有80um横向扇出，折返终点比MUX端口更靠右；
            # 必须保留足够返程直线回到统一端口平面，不能只取20um。
            segment=substring(line,0,distance+120)
            coords=list(segment.coords);coords.insert(0,(coords[0][0]-8,coords[0][1]))
            segment=LineString(coords);ny=line.coords[0][1]
        cut=segment.intersection(LineString([(plane,-20000),(plane,20000)]))
        points=list(cut.geoms) if hasattr(cut,'geoms') else [cut]
        ys=[p.y for p in points if p.geom_type=='Point']
        if len(ys)!=2:raise ValueError(f'{name}不能唯一确定同平面的两个端口：{cut}')
        iy=min(ys,key=lambda v:abs(v-ny));oy=max(ys,key=lambda v:abs(v-ny))
        near.append(iy);far.append(oy);lines.append(segment)
    origin_y=sum(near)/2
    moved=[translate(g,xoff=-plane,yoff=-origin_y) for g in lines]
    clip=box(-2000,-2000,6,2000) if left else box(-6,-2000,2000,2000)
    moved=[g.intersection(clip) for g in moved]
    if any(g.geom_type!='LineString' for g in moved):raise ValueError('PML外截断后路线不连续')
    core=unary_union([g.buffer(.35,cap_style=2,join_style=1) for g in moved])
    slab=unary_union([g.buffer(3.6,cap_style=2,join_style=1) for g in moved])
    bounds=list(slab.bounds);bounds=[bounds[0]-5,bounds[1]-5,bounds[2]+5,bounds[3]+5]
    bounds[2 if left else 0]=4 if left else -4
    direction='Backward' if left else 'Forward'
    ports={'near_pair':['x',0.,0.,direction],
           'far_a':['x',0.,far[0]-origin_y,direction],
           'far_b':['x',0.,far[1]-origin_y,direction]}
    data={'scope':'实际中心线的成对接入，有限7.2um残余平台为名义模型，非真实LN2掩膜',
          'pair':a.pair,'route_names':names,'origin_global_um':[plane,origin_y],
          'near_centers_um':[v-origin_y for v in near],'far_centers_um':[v-origin_y for v in far],
          'ports':ports,'port_spans_um':{'near_pair':14.,'far_a':10.,'far_b':10.},'bounds_um':bounds,
          'centerlines_um':[list(g.coords) for g in moved],
          'core_polygons_um':polygons(core),'slab_polygons_um':polygons(slab),'source_route_json':str(a.routes),
          'notes':['用宽端口的复模基底读出近邻两路，不能把宽端口mode1当成某一条单独波导',
                   'far_a/far_b分别激励；近邻局域导模幅度需从复模场重构/投影，不按单个S元素直接判断串扰']}
    a.output.parent.mkdir(parents=True,exist_ok=True)
    with a.output.open('x',encoding='utf-8') as f:json.dump(data,f,ensure_ascii=False,indent=2)
    print(json.dumps({k:v for k,v in data.items() if 'polygons' not in k and 'centerlines' not in k},ensure_ascii=False,indent=2))


if __name__=='__main__':main()
