"""在已确认的v11工作稿上替换无源路由为欧拉弯，并另存新GDS。

运行前提：时延预算文件已经同时通过方向群时延闭合和中心线拓扑筛查。
脚本不改写输入GDS，只替换顶层20/0、20/1、10/2的直接图形；已有MUX、
有源区、电极、交叉器、拉锥和端面耦合器实例保持不变。
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path

from klink import KLinkClient
from shapely.geometry import LineString, box


ROOT=Path(__file__).resolve().parents[2]
DEFAULT_SOURCE=(ROOT/'results/layout/_历史归档/迭代版本_20260912/MUX外端公共窗口加宽_工作稿_v11_01/'
                '四程10GHz_15mm_MUX外端公共窗口加宽_待验证工作稿.gds')
DEFAULT_BUDGET=(ROOT/'results/layout/_历史归档/迭代版本_20260912/欧拉回路时延预算_v3_无斜直线/'
                '欧拉回路_10GHz方向时延预算.json')
SOURCE_TOP='EO4P_10G_15MM_EXTERNAL_WINDOWS_WIDENED_DRAFT'
OUTPUT_TOP='EO4P_10G_15MM_EULER_DELAY_CLOSED_DRAFT'
KLAYOUT_EXE=Path(r'C:/Users/PC/AppData/Roaming/KLayout/klayout_app.exe')
KLAYOUT_WORKER=Path(__file__).with_name('klayout_replace_euler.py')
LAYERS={'TRENCH_SIN':'10/2','WGCORE_LN1':'20/0','WGCLAD_LN1':'20/1'}
WIDTHS={'TRENCH_SIN':17.2,'WGCORE_LN1':.7,'WGCLAD_LN1':7.2}


def sha256(path: Path) -> str:
    digest=hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda:stream.read(1024*1024),b''):
            digest.update(chunk)
    return digest.hexdigest()


def parse_args() -> argparse.Namespace:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,default=DEFAULT_SOURCE)
    parser.add_argument('--budget',type=Path,default=DEFAULT_BUDGET)
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--host',default='127.0.0.1')
    parser.add_argument('--port',type=int,default=8765)
    parser.add_argument('--klayout-exe',type=Path,default=KLAYOUT_EXE)
    return parser.parse_args()


def as_lines(geometry) -> list[LineString]:
    if geometry.geom_type=='LineString':
        return [geometry]
    if geometry.geom_type=='MultiLineString':
        return list(geometry.geoms)
    if geometry.geom_type=='GeometryCollection':
        return [item for item in geometry.geoms if item.geom_type=='LineString']
    raise RuntimeError(f'裁切后得到不支持的几何类型：{geometry.geom_type}')


def direct_route_pieces(routes: dict) -> dict[str,list[LineString]]:
    """剔除由层次化实例承载的交叉器和端口拉锥区段。"""
    pieces={
        'loop1':[LineString(routes['loop1'])],
        'loop3':[LineString(routes['loop3'])],
        'output':[LineString([[18150.,-26.6],[18600.,-26.6]])],
        'input_west':[LineString([[100.,250.],[177.5,250.]])],
    }
    crossing_window=box(177.5,127.5,422.5,372.5)
    pieces['loop2']=as_lines(LineString(routes['loop2']).difference(crossing_window))
    input_east=LineString(routes['input']).difference(box(-1e6,127.5,422.5,372.5))
    pieces['input_east']=as_lines(input_east)
    expected={'loop1':1,'loop2':2,'loop3':1,'output':1,'input_west':1,'input_east':1}
    actual={key:len(value) for key,value in pieces.items()}
    if actual!=expected:
        raise RuntimeError(f'器件窗口裁切结果异常：{actual}，预期{expected}')
    return pieces


def geometry_for_piece(centerline: LineString,width: float,layer: str):
    polygon=centerline.buffer(width/2,quad_segs=8,cap_style=2,join_style=1)
    polygon=polygon.simplify(.001,preserve_topology=True)
    if polygon.geom_type!='Polygon' or not polygon.is_valid or polygon.is_empty:
        raise RuntimeError(f'{layer}缓冲后不是有效单多边形')
    return polygon


def polygon_points(polygon) -> list[list[float]]:
    return [[float(x),float(y)] for x,y in list(polygon.exterior.coords)[:-1]]


def main() -> None:
    args=parse_args()
    args.source=args.source.resolve(); args.budget=args.budget.resolve()
    args.output_dir=args.output_dir.resolve()
    if not args.source.is_file() or not args.budget.is_file():
        raise FileNotFoundError('输入GDS或时延预算不存在')
    if args.output_dir.exists():
        raise FileExistsError(f'输出目录已存在，拒绝覆盖：{args.output_dir}')
    budget=json.loads(args.budget.read_text(encoding='utf-8'))
    if not budget.get('delay_closed') or not budget.get('layout_generation_allowed'):
        raise RuntimeError('时延或拓扑门槛未关闭，拒绝生成欧拉版图')
    topology=budget.get('topology_screen',{})
    if not topology.get('topology_screen_pass'):
        raise RuntimeError('中心线拓扑筛查未通过')
    routes=budget['route_points_um']
    pieces=direct_route_pieces(routes)
    long_diagonal=[]
    for route_name,points in routes.items():
        for first,second in zip(points,points[1:]):
            dx=second[0]-first[0]; dy=second[1]-first[1]
            if (dx*dx+dy*dy)**.5>2.0 and abs(dx)>1e-9 and abs(dy)>1e-9:
                long_diagonal.append({'route':route_name,'start':first,'end':second})
    if long_diagonal:
        raise RuntimeError(f'发现长斜直线，拒绝生成版图：{long_diagonal}')

    args.output_dir.mkdir(parents=True,exist_ok=False)
    output_gds=args.output_dir/'四程10GHz_15mm_欧拉弯_时延闭合_待验证工作稿.gds'
    full_png=args.output_dir/'欧拉弯四程版图_全图.png'
    left_png=args.output_dir/'欧拉输入与三级落回_局部.png'
    right_png=args.output_dir/'欧拉右端折返_局部.png'
    report_path=args.output_dir/'欧拉版图生成检查.json'
    manifest_path=args.output_dir/'欧拉路由多边形清单.json'

    vertices={}; polygons={layer:[] for layer in LAYERS.values()}
    inserted_per_layer={layer:0 for layer in LAYERS.values()}
    for logical,layer in LAYERS.items():
        geometries=[]
        for route_name,segments in pieces.items():
            for index,segment in enumerate(segments,1):
                geometries.append((f'{route_name}:{index}',
                                   geometry_for_piece(segment,WIDTHS[logical],layer)))
        components=[(name,geometry) for name,geometry in geometries]
        if logical=='WGCLAD_LN1':
            # 近邻端口从公共窗口分开时，以规则裕量充足的矩形公共窗口过渡。
            # 矩形终点处两条7.2um窗口的净间距均大于0.30um。
            common_boxes=[
                box(785.,26.5,851.,42.),box(785.,-42.,851.,-26.5),
                box(18149.,23.,18216.,38.1),box(18149.,-38.1,18216.,-23.),
            ]
            components.extend((f'端口公共窗口:{index}',geometry)
                              for index,geometry in enumerate(common_boxes,1))
        elif logical=='TRENCH_SIN':
            # 矩形终点处两条17.2um窗口的净间距均大于0.20um。
            common_boxes=[
                box(757.,21.5,851.,57.1),box(757.,-57.1,851.,-21.5),
                box(18149.,18.,18244.,53.4),box(18149.,-53.4,18244.,-18.),
            ]
            components.extend((f'端口公共窗口:{index}',geometry)
                              for index,geometry in enumerate(common_boxes,1))
        for name,geometry in components:
            geometry=geometry.simplify(.001,preserve_topology=True)
            if geometry.geom_type!='Polygon' or not geometry.is_valid or geometry.is_empty:
                raise RuntimeError(f'{layer}公共窗口处理后不是有效单多边形')
            points=polygon_points(geometry)
            polygons[layer].append(points)
            vertices[f'{layer}:{name}']=len(points)
            inserted_per_layer[layer]+=1
    manifest_path.write_text(json.dumps({
        'source_gds':str(args.source),'source_top':SOURCE_TOP,
        'output_top':OUTPUT_TOP,'polygon_simplification_um':.001,
        'polygons':polygons,
    },ensure_ascii=False),encoding='utf-8')

    client=KLinkClient(host=args.host,port=args.port)
    client.connect()
    try:
        client.layout_show_file(str(args.source),mode='replace',keep_position=False)
        info=client.layout_info('full')
        if info.get('file')!=str(args.source) or info.get('cell')!=SOURCE_TOP:
            raise RuntimeError(f'KLayout载入对象不符：{info}')
        instances_before=client.instance_query(SOURCE_TOP,limit=100)
        instance_signature=[(item['child'],item['bbox_dbu'],item.get('klink_id'))
                            for item in instances_before['instances']]
        dry=client.shape_delete(SOURCE_TOP,layers=list(LAYERS.values()),
                                limit=1000,dry_run=True)
        if dry.get('matched')!=13 or dry.get('per_layer')!={'10/2':3,'20/0':7,'20/1':3}:
            raise RuntimeError(f'顶层待替换图形与v11基线不符：{dry}')
    finally:
        client.close()

    args.klayout_exe=args.klayout_exe.resolve()
    if not args.klayout_exe.is_file() or not KLAYOUT_WORKER.is_file():
        raise FileNotFoundError('KLayout批处理程序或欧拉替换脚本不存在')
    command=[str(args.klayout_exe),'-b','-r',str(KLAYOUT_WORKER.resolve()),
             '-rd',f'source={args.source}','-rd',f'output={output_gds}',
             '-rd',f'manifest={manifest_path}','-rd',f'source_top={SOURCE_TOP}',
             '-rd',f'output_top={OUTPUT_TOP}']
    completed=subprocess.run(command,capture_output=True,text=True,encoding='utf-8',
                             errors='replace',timeout=120,check=False)
    if completed.returncode!=0 or not output_gds.is_file():
        raise RuntimeError('KLayout批处理写入失败：\n'+completed.stdout+'\n'+completed.stderr)

    client=KLinkClient(host=args.host,port=args.port)
    client.connect()
    try:
        client.layout_show_file(str(output_gds),mode='replace',keep_position=False)
        output_info=client.layout_info('full')
        if output_info.get('file')!=str(output_gds) or output_info.get('cell')!=OUTPUT_TOP:
            raise RuntimeError(f'生成GDS回读对象不符：{output_info}')
        instances_after=client.instance_query(OUTPUT_TOP,limit=100)
        after_signature=[(item['child'],item['bbox_dbu'],item.get('klink_id'))
                         for item in instances_after['instances']]
        if after_signature!=instance_signature:
            raise RuntimeError('层次化器件实例在路由替换过程中发生变化')
        direct_counts={}
        for layer in LAYERS.values():
            query=client.shape_query(OUTPUT_TOP,layers=[layer],limit=100)
            direct_counts[layer]=query['returned']
            if query['truncated'] or query['returned']!=inserted_per_layer[layer]:
                raise RuntimeError(f'{layer}顶层直写图形数量异常：{query["returned"]}')

        client.call('view.show_cell',{'cell':OUTPUT_TOP,'zoom_fit':True})
        lyp=ROOT/'ioptee_sin_tfln_v1.0_beta/ioptee_sin_tfln_v1.0_beta/techfiles/klayout/ioptee_sin_tfln.lyp'
        if lyp.is_file():
            client.call('layer.load_lyp',{'path':str(lyp.resolve())})
        client.call('view.hier_levels',{'min':0,'max':10})
        client.call('layer.set_visible',{'layers':['101/0'],'visible':False})
        client.call('view.zoom_fit',{})
        client.call('view.screenshot',{'mode':'path','path':str(full_png),
                    'width_px':2400,'height_px':1100})
        client.call('view.screenshot',{'mode':'path','path':str(left_png),
                    'bbox_um':[0.,-80.,1800.,620.],
                    'width_px':1800,'height_px':900})
        client.call('view.screenshot',{'mode':'path','path':str(right_png),
                    'bbox_um':[16700.,-650.,18900.,1100.],
                    'width_px':1500,'height_px':1200})
        file_info=client.call('layout.file_info',{'path':str(output_gds),'detail':'counts'})
    finally:
        client.close()

    report={
        'status':'延时先闭合后生成的欧拉弯名义工作稿；不是流片签核版',
        'source_gds':str(args.source),'source_sha256':sha256(args.source),
        'delay_budget':str(args.budget),'delay_budget_sha256':sha256(args.budget),
        'output_gds':str(output_gds),'output_sha256':sha256(output_gds),
        'delay_closed':budget['delay_closed'],
        'topology_screen_pass':topology['topology_screen_pass'],
        'long_diagonal_straight_segment_count':len(long_diagonal),
        'long_segment_threshold_um':2.0,
        'target_delays_ps':[row['target_total_delay_ps'] for row in budget['loops']],
        'actual_delays_ps':[row['total_delay_ps'] for row in budget['loops']],
        'route_lengths_um':[row['solved_euler_route_length_um'] for row in budget['loops']],
        'minimum_radius_um':min(row['minimum_radius_um'] for row in budget['loops']),
        'euler_fraction':budget['euler']['fraction'],
        'euler_bend_count':sum(row['euler_bend_count'] for row in budget['loops'])+2,
        'direct_route_piece_counts':{key:len(value) for key,value in pieces.items()},
        'deleted_v11_top_shapes_preflight':dry,'inserted_top_shapes':inserted_per_layer,
        'direct_shape_counts_after':direct_counts,
        'hierarchical_instances_preserved':len(instance_signature),
        'polygon_simplification_um':.001,'polygon_vertices':vertices,
        'local_common_window_boxes_um':{
            '20/1':[[785.,26.5,851.,42.],[785.,-42.,851.,-26.5],
                    [18149.,23.,18216.,38.1],[18149.,-38.1,18216.,-23.]],
            '10/2':[[757.,21.5,851.,57.1],[757.,-57.1,851.,-21.5],
                    [18149.,18.,18244.,53.4],[18149.,-53.4,18244.,-18.]],
        },
        'file_info':file_info,
        'klayout_batch_stdout':completed.stdout,'klayout_batch_stderr':completed.stderr,
        'screenshots':[str(full_png),str(left_png),str(right_png)],
        'not_validated':['欧拉弯传播损耗和模式串扰','交叉器真实群时延和损耗',
                         '完整级联光学传输','工艺容差','代工方最终DRC','射频焊盘与终端'],
    }
    report_path.write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps({'output_gds':str(output_gds),'report':str(report_path),
                      'actual_delays_ps':report['actual_delays_ps'],
                      'route_lengths_um':report['route_lengths_um'],
                      'topology_screen_pass':report['topology_screen_pass'],
                      'direct_shape_counts_after':direct_counts},
                     ensure_ascii=False,indent=2))


if __name__=='__main__':
    main()
