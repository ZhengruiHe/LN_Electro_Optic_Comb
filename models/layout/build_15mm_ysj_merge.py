"""构建100um外接续的15mm功能候选，原样保留YSJ单元记录合并进一个block。"""
import argparse
import hashlib
import json
import math
import struct
import subprocess
from pathlib import Path
from shapely.geometry import LineString, box
from shapely.ops import unary_union
from route_geometry import check_route_network
from preview_15mm_upper_candidate import polygons

ROOT=Path(__file__).resolve().parents[2]
KLAYOUT=Path('C:/Users/PC/AppData/Roaming/KLayout/klayout_app.exe')
SOURCE=ROOT/'results/layout/双端GSG150_窗口修正_v18_11/四程10GHz_15mm_GSG150um_分叉窗口修正_待验证工作稿.gds'


def gds_records(data):
    pos=0
    while pos<len(data):
        if pos+4>len(data):raise ValueError('GDS记录头不完整')
        length,kind,dtype=struct.unpack('>HBB',data[pos:pos+4])
        if length<4 or length%2 or pos+length>len(data):raise ValueError('GDS记录长度无效')
        yield kind,dtype,data[pos:pos+length]
        pos+=length


def structures(data):
    cells={};current=None;name=None;units=None;endlib=None
    for kind,dtype,record in gds_records(data):
        if kind==3:units=record[4:]
        if kind==5:current=[];name=None
        if current is not None:current.append(record)
        if kind==6 and current is not None:name=record[4:].rstrip(b'\0').decode('ascii')
        if kind==7:
            if not name or name in cells:raise ValueError('GDS单元名缺失或重复')
            cells[name]=b''.join(current);current=None
        if kind==4:endlib=record
    if current is not None or endlib is None or units is None:raise ValueError('GDS库结构不完整')
    return cells,units,endlib


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--routes',type=Path,required=True)
    p.add_argument('--ysj',type=Path,required=True)
    p.add_argument('--output-dir',type=Path,required=True)
    a=p.parse_args();a.output_dir=a.output_dir.resolve()
    a.output_dir.mkdir(parents=True,exist_ok=False)
    source_report=json.loads((SOURCE.parent/'分叉窗口与GSG检查.json').read_text(encoding='utf-8'))
    if hashlib.sha256(SOURCE.read_bytes()).hexdigest()!=source_report['output_sha256']:
        raise ValueError('15mm源GDS与已复核基线不一致，拒绝继续')
    data=json.loads(a.routes.read_text(encoding='utf-8'))
    if not data['input_and_final_output_included']:raise ValueError('缺少输入输出')
    allowed=[]
    for c in data['crossing_details']:
        if abs(c.get('angle_deg',0)-90)>1e-3 or not all(c['straight_245um_each_route'].values()):
            raise ValueError('交点未具备正交直段，不得放置交叉黑盒')
        allowed.append({'routes':c['routes'],'center_um':[round(c['x_um'],3),round(c['y_um'],3)],'angle_deg':90.})
    topology=check_route_network(data['routes'],allowed)
    if not topology['topology_screen_pass']:raise ValueError('存在未指定交叉或自交')
    centers=sorted({tuple(c['center_um']) for c in allowed})
    expected_count=int(data.get('expected_crossing_count',7))
    if len(centers)!=expected_count:raise ValueError('交叉数量与批准的路线记录不一致')
    groups=[]
    for x in sorted({x for x,y in centers}):
        ys=sorted(y for xx,y in centers if xx==x)
        vertical_names=set()
        for spec in allowed:
            if spec['center_um'][0]!=x:continue
            point_y=spec['center_um'][1]
            from shapely.geometry import Point
            for name in spec['routes']:
                line=LineString(data['routes'][name]);s=line.project(Point(x,point_y))
                first=line.interpolate(s-1);last=line.interpolate(s+1)
                if abs(first.x-last.x)<1e-6:vertical_names.add(name)
        if len(vertical_names)!=1:raise ValueError('同列交叉不能确认共用同一竖向波导')
        line=LineString(data['routes'][next(iter(vertical_names))])
        s0=line.project(Point(x,ys[0]));s1=line.project(Point(x,ys[-1]))
        if abs(abs(s1-s0)-(ys[-1]-ys[0]))>1e-3:raise ValueError('串接交叉之间不是同一直线')
        if any(y1-y0<45 for y0,y1 in zip(ys,ys[1:])):raise ValueError('交叉器黑盒互相重叠')
        groups.append((x,ys))
    col_x=max(x for x,ys in groups)
    active_um=float(data.get('active_taper_um',200.))
    external_um=float(data.get('external_taper_um',100.))
    revised=active_um!=200. or external_um!=100.
    if not 100<=active_um<=200 or not 50<=external_um<=100:
        raise ValueError('接续长度超出本轮已检查范围')
    clear_boxes=[box(x-22.5,y-22.5,x+22.5,y+22.5) for x,y in centers]
    tapers=[]
    for x,y in centers:
        clear_boxes.append(box(x-122.5,y-8.6,x+122.5,y+8.6))
        for label,angle,origin in [('W',0,[x-122.5,y]),('E',180,[x+122.5,y])]:
            tapers.append({'name':f'C_{x:g}_{y:g}_{label}','angle_deg':angle,'narrow_um':origin})
    for x,ys in groups:
        clear_boxes.append(box(x-8.6,ys[0]-122.5,x+8.6,ys[-1]+122.5))
        tapers += [{'name':f'bus_{x:g}_S','angle_deg':90,'narrow_um':[x,ys[0]-122.5]},
                   {'name':f'bus_{x:g}_N','angle_deg':270,'narrow_um':[x,ys[-1]+122.5]}]
    edge_ports=data.get('edge_coupler_internal_ports_um',[[-10380.,-1700.],[-10380.,-500.]])
    for label,(x,y) in zip(['input','output'],edge_ports):
        if x!=-10380:raise ValueError('此构建器仅接受左端面耦合器')
        clear_boxes.append(box(x,y-8.6,x+100,y+8.6))
        tapers.append({'name':'edge_'+label,'angle_deg':180,'narrow_um':[x+100,y]})
    clear=unary_union(clear_boxes)
    bb=unary_union([box(x-22.5,y-22.5,x+22.5,y+22.5) for x,y in centers])
    direct=[]
    # PDK LN模板的7.2/17.2是固定窗口总宽，不随自定义核心宽度缩放。
    route_widths=[((20,0),.7),((20,1),7.2),((10,2),17.2)]
    for key,width in route_widths:
        region=unary_union([LineString(points).buffer(width/2,cap_style=2,join_style=2) for points in data['routes'].values()]).difference(clear)
        bus_width={(20,0):1.2,(20,1):7.2,(10,2):17.2}[key]
        bus=unary_union([box(x-bus_width/2,ys[0]-22.5,x+bus_width/2,ys[-1]+22.5) for x,ys in groups]).difference(bb)
        region=unary_union([region,bus])
        direct.append({'layer':list(key),'polygons':polygons(region)})
    ours=a.output_dir/('四程15mm_周期匹配_布线候选.gds' if revised else '四程15mm_外接续100_功能版图检查稿.gds')
    top=f'EO4P_15MM_GSG150_A{active_um:g}_E{external_um:g}_ROUTED_DRAFT' if revised else 'EO4P_15MM_GSG150_EXT100_ROUTED_DRAFT'
    horizontal=data.get('layout_variant')=='horizontal_access'
    if horizontal:top='EO4P_15MM_GSG150_HORIZONTAL_DRAFT'
    ysjd=a.ysj.read_bytes();ysj_cells,ysj_units,endlib=structures(ysjd)
    ysjtop='MAIN_HETERO_V1_RTC_COMPACT_SPC_R7A'
    if ysjtop not in ysj_cells:raise ValueError('不是指定YSJ顶层')
    payload={'source_gds':str(SOURCE),'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
             'top_name':top,'shift_um':data['shift_um'],'route_geometry':direct,
             'active_taper_um':active_um,'external_taper_um':external_um,
             'normalize_mux_windows_to_pdk':bool(data.get('normalize_mux_windows_to_pdk',False)),
             'loop_periods':data.get('loop_periods',[3.,3.5,3.]),
             'left_cleanup_xmax_um':col_x+400 if revised else -10500.,
             'window_merge_caps':data.get('window_merge_caps',[]),
             'additional_cleanup_rois':data.get('additional_cleanup_rois',[]),
             'attachment_points_um':[data['routes'][n][-1] for n in ['loop1','loop2','loop3','input']]+
                                    [data['routes'][n][0] for n in ['loop1','loop2','loop3','output']],
             'crossing_centers_um':centers,'edge_coupler_internal_ports_um':edge_ports,
             'tapers':tapers,'bus_link_core_boxes_um':[[x-.6,y0+22.5,x+.6,y1-22.5] for x,ys in groups for y0,y1 in zip(ys,ys[1:])],
             'ours_gds':str(ours),'report_json':str(a.output_dir/'原生连接检查.json'),
             'preview_json':str(a.output_dir/'四程实际掩膜预览.json'),
             'parent_gds':str(a.output_dir/'父单元临时库.gds'),'ysj_top':ysjtop,
             'merged_top':('YSJ_EO4P_15MM_HORIZONTAL_DRAFT' if horizontal else 'YSJ_EO4P_15MM_CYCLE_MATCHED_DRAFT') if revised else 'YSJ_EO4P_15MM_BLOCK_21800X3800_DRAFT'}
    payload_path=a.output_dir/'构建参数.json'
    payload_path.write_text(json.dumps(payload,ensure_ascii=False),encoding='utf-8')
    r=subprocess.run([str(KLAYOUT),'-b','-r',str(Path(__file__).with_name('klayout_build_upper_15mm.py')),
                      '-rd','payload_file='+str(payload_path)],capture_output=True,text=True,encoding='gb18030',errors='replace')
    report=json.loads(Path(payload['report_json']).read_text(encoding='utf-8')) if Path(payload['report_json']).exists() else {}
    print(json.dumps({k:v for k,v in report.items() if k not in ('taper_checks','window_marker_boxes_um')},ensure_ascii=False,indent=2))
    if r.returncode:
        print(r.stdout[-5000:]);print(r.stderr[-4000:]);raise RuntimeError('连接检查未通过，未合并YSJ')
    own_cells,own_units,_=structures(ours.read_bytes())
    parent_cells,parent_units,_=structures(Path(payload['parent_gds']).read_bytes())
    if own_units!=ysj_units or parent_units!=ysj_units:raise ValueError('GDS单位记录不同，拒绝直接拼接')
    if set(ysj_cells)&set(own_cells):raise ValueError('单元重名，拒绝覆盖原YSJ')
    if not ysjd.endswith(endlib):raise ValueError('YSJ ENDLIB后存在额外内容')
    combined=ysjd[:-len(endlib)]+b''.join(own_cells.values())+parent_cells[payload['merged_top']]+endlib
    merged=a.output_dir/('YSJ与四程15mm_周期匹配_布线候选.gds' if revised else 'YSJ与四程15mm_同block_功能版图候选.gds')
    with merged.open('xb') as f:f.write(combined)
    merged_cells,_,_=structures(merged.read_bytes())
    preserved={name:merged_cells[name]==raw for name,raw in ysj_cells.items()}
    manifest={'status':'merged_geometry_candidate_not_tapeout_signoff','merged_gds':str(merged),
              'merged_top':payload['merged_top'],'YSJ_cell_records_byte_identical':preserved,
              'YSJ_source_unchanged':a.ysj.read_bytes()==ysjd,
              'source_15mm_unchanged':hashlib.sha256(SOURCE.read_bytes()).hexdigest()==payload['source_sha256'],
              'sha256':hashlib.sha256(combined).hexdigest(),'crossings':allowed,'topology':topology,
              'crossing_count':len(centers),'crossing_traversals':len(centers)*2,
              'active_taper_um':active_um,'external_taper_um':external_um,
              'loop_periods':data.get('loop_periods',[3.,3.5,3.]),
              'native_drc_counts':{k:v for k,v in report.get('rules',{}).items() if k in ('core_width_030','core_space_030','etch_width_030','etch_space_030','trench_width_020','trench_space_020','metal_width_2','metal_space_3')},
              'layout_variant':data.get('layout_variant','cycle_matched'),
              'edge_coupler_internal_ports_um':edge_ports,
              'interface_taper_count':len(tapers),
              'serial_bus_length_um':sum(y1-y0-45. for x,ys in groups for y0,y1 in zip(ys,ys[1:])),
              'limitations':['黑盒及1.2um连接总线实际群时延未知，当前预算沿用0.7um方向ng占位',
                             f'{active_um:g}/{external_um:g}um外部接续；750um耦合核心及原辅助尾端保持，无新MUX场求解',
                             'LN2暂缓，但原YSJ已有21层完整保留','完整光学性能和GSG S参数未签核']}
    (a.output_dir/'合并清单.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
    if not all(preserved.values()) or not manifest['YSJ_source_unchanged'] or not manifest['source_15mm_unchanged']:
        raise RuntimeError('源记录保留检查失败')
    print('merged='+str(merged))


if __name__=='__main__':main()
