"""从已有v16中心线和已保存群折射率直接核算全光路；不调用求解器、不改GDS。"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import json
import math
from pathlib import Path

from shapely.geometry import LineString, Point
from shapely.ops import substring

from delay_budget import route_delay

ROOT = Path(__file__).resolve().parents[1]
LAYOUT_DIR = ROOT / 'results/layout/欧拉弯_时延闭合_无斜直线_工作稿_v16_01'
BUDGET = ROOT / 'results/layout/欧拉回路时延预算_v3_无斜直线/欧拉回路_10GHz方向时延预算.json'
MUX = ROOT / 'results/optical/模式复用器_70度_EME_群时延摘要.json'


def calculate(budget, mux):
    """参考面从输入端面耦合器片内端口，到输出端面耦合器片内端口。"""
    ng_y = budget['passive_rib']['ng_crystal_Y']
    ng_z = budget['passive_rib']['ng_crystal_Z']
    routes = budget['route_points_um']
    modes = ['TE0', 'TE1', 'TE0', 'TE1']
    loop_rows = budget['loops']
    active = {'TE0': loop_rows[0]['fixed_active_delay_ps'],
              'TE1': loop_rows[1]['fixed_active_delay_ps']}
    active_taper = {'TE0': loop_rows[0]['fixed_two_active_taper_delay_ps']/2,
                    'TE1': loop_rows[1]['fixed_two_active_taper_delay_ps']/2}
    external_taper = loop_rows[0]['fixed_two_external_taper_delay_ps']/2
    channels = mux['channels_at_center']
    ledger = []

    def add(block, kind, name, delay, length=None, basis='已有模型工程估算'):
        ledger.append({'block': block, 'kind': kind, 'name': name,
                       'length_um': length, 'delay_ps': delay, 'basis': basis})

    def add_route(block, name, line, intervals=()):
        # intervals中的区段已经包含在整条路由中，拆分列出，禁止重复叠加。
        intervals = sorted(intervals)
        cursor = 0.0
        def regular(first, last):
            if last-first < 1e-9:
                return
            points = list(substring(line, first, last).coords)
            straight_delay = bend_delay = straight_length = bend_length = 0.0
            for p, q in zip(points, points[1:]):
                length = math.dist(p, q)
                dx, dy = q[0]-p[0], q[1]-p[1]
                item = route_delay([p, q], ng_y, ng_z)
                if abs(dx)<1e-9 or abs(dy)<1e-9:
                    straight_length += length
                    straight_delay += item['directional_delay_ps']
                else:
                    if length>2:
                        raise ValueError('中心线存在超过2um的斜直线')
                    bend_length += length
                    bend_delay += item['directional_delay_ps']
            if straight_length:
                add(block, 'passive_straight', name+'：直段', straight_delay,
                    straight_length, '保存的0.7um波导Y/Z群折射率乘实际中心线长度')
            if bend_length:
                add(block, 'euler_bend', name+'：欧拉弯', bend_delay,
                    bend_length, '沿保存中心线积分ngY*cos²θ+ngZ*sin²θ；工程插值')
        for start, stop, kind, label in intervals:
            if start<cursor-1e-8 or stop<=start or stop>line.length+1e-8:
                raise ValueError('路由细分区间重叠、反向或越界')
            regular(cursor, start)
            part = substring(line, start, stop)
            delay = route_delay(list(part.coords), ng_y, ng_z)['directional_delay_ps']
            add(block, kind, label, delay, part.length,
                '沿用0.7um同方向ng占位；该器件真实ng未单独提取')
            cursor = stop
        regular(cursor, line.length)

    def mux_in(block, mode):
        add(block, 'mux_external_taper', 'MUX输入侧200um外端拉锥', external_taper, 200)
        channel = 'TE0_反向_s14' if mode=='TE0' else '辅助转TE1_反向_s25'
        add(block, 'mux_body', 'MUX输入：'+channel,
            channels[channel]['group_delay_ps'], 750, '沿用已有750um本体EME群时延')
        add(block, 'mux_active_taper', 'MUX输入侧200um有源拉锥：'+mode,
            active_taper[mode], 200)

    def mux_out(block, mode):
        add(block, 'mux_active_taper', 'MUX输出侧200um有源拉锥：'+mode,
            active_taper[mode], 200)
        channel = 'TE0_正向_s41' if mode=='TE0' else 'TE1转辅助_正向_s52'
        add(block, 'mux_body', 'MUX输出：'+channel,
            channels[channel]['group_delay_ps'], 750, '沿用已有750um本体EME群时延')
        add(block, 'mux_external_taper', 'MUX输出侧200um外端拉锥', external_taper, 200)

    add_route('输入至第1程', '输入引线', LineString(routes['input']), [
        (0, 100, 'interface_taper', '输入端面耦合器接续拉锥'),
        (177.5, 277.5, 'interface_taper', '交叉器西侧拉锥'),
        (277.5, 322.5, 'crossing_blackbox', '交叉器横向通过'),
        (322.5, 422.5, 'interface_taper', '交叉器东侧拉锥'),
    ])
    mux_in('输入至第1程', 'TE0')
    for i, mode in enumerate(modes):
        block = f'第{i+1}程至第{i+2}程' if i<3 else '第4程至输出'
        add(block, 'active', f'第{i+1}程15mm有源区：{mode}', active[mode], 15000,
            '沿用已保存有源模式群折射率')
        mux_out(block, mode)
        if i<3:
            line = LineString(routes['loop'+str(i+1)])
            intervals = []
            if i==1:
                distances = [line.project(Point(300, y)) for y in (372.5,272.5,227.5,127.5)]
                for start, stop, kind, label in zip(distances,distances[1:],
                    ['interface_taper','crossing_blackbox','interface_taper'],
                    ['交叉器北侧拉锥','交叉器纵向通过','交叉器南侧拉锥']):
                    intervals.append((start,stop,kind,label))
            add_route(block, f'无源回路{i+1}', line, intervals)
            mux_in(block, modes[i+1])
        else:
            add_route(block, '输出引线', LineString(routes['output']), [
                (450,550,'interface_taper','输出端面耦合器接续拉锥')])

    block_totals = {}
    cumulative = 0.0
    for item in ledger:
        item['start_time_ps'] = cumulative
        cumulative += item['delay_ps']
        item['end_time_ps'] = cumulative
        block_totals[item['block']] = block_totals.get(item['block'],0)+item['delay_ps']
    times = [item['start_time_ps'] for item in ledger if item['kind']=='active']
    intervals = [times[i+1]-times[i] for i in range(3)]
    counts = Counter(item['kind'] for item in ledger)
    if any(counts[k]!=n for k,n in {'active':4,'mux_body':8,'mux_active_taper':8,
        'mux_external_taper':8,'interface_taper':6,'crossing_blackbox':2}.items()):
        raise ValueError('光程中器件通过次数不闭合')
    reference_plane_length = 4*15000+8*(750+200+200)+sum(
        LineString(points).length for points in routes.values())
    if not math.isclose(sum(item['length_um'] for item in ledger),reference_plane_length,abs_tol=1e-6):
        raise ValueError('中心线总长度不守恒，可能漏计或重复计算器件')
    return {'status':'当前版图全光路名义时延直接核算；未重跑MUX或欧拉弯求解器',
        'reference_planes':{'input_um':[0,250],'output_um':[18700,-26.6],
                            'definition':'两只端面耦合器的片内端口；不含黑盒端面耦合器本体'},
        'ledger':ledger, 'component_traversal_counts':dict(counts),
        'count_convention':'passive_straight和euler_bend统计分类合计行；其余统计实际器件通过次数',
        'individual_euler_bend_count':sum(len(meta['bends'])
            for meta in budget['route_metadata'].values() if 'bends' in meta),
        'block_delays_ps':block_totals,'active_entrance_times_ps':times,
        'interpass_delays_ps':intervals,
        'interpass_phase_residual_deg_at_10GHz':[
            (actual-target)*3.6 for actual,target in zip(intervals,[300,350,300])],
        'port_to_port_delay_ps':cumulative,'port_to_port_centerline_length_um':reference_plane_length,
        'facet_to_facet_delay_ps':None,
        'limitations':['MUX本体及其时延沿用既有结果，按用户要求不重新仿真',
                       '任意晶向采用cos²插值；不包含曲率对局部群折射率的独立修正',
                       '交叉黑盒通过两次，六个接口拉锥分别计入；仍以0.7um同晶向ng占位',
                       'MUX两侧拉锥沿用平均ng；外端主/辅助支路缺少独立拉锥群时延',
                       '各回路出口与下一程入口的有源拉锥分别按相应模式计入，较原两端同ng预算有亚飞秒修正',
                       '两只端面耦合器黑盒本体时延未知，未计入；不能将端口间时延称为端面到端面时延']}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--budget',type=Path,default=BUDGET)
    parser.add_argument('--mux',type=Path,default=MUX)
    parser.add_argument('--layout-report',type=Path,default=LAYOUT_DIR/'欧拉版图生成检查.json')
    parser.add_argument('--output-dir',type=Path,required=True)
    args = parser.parse_args()
    report = json.loads(args.layout_report.read_text(encoding='utf-8'))
    gds = Path(report['output_gds'])
    digest = hashlib.sha256(gds.read_bytes()).hexdigest()
    if digest!=report['output_sha256']:
        raise ValueError('当前GDS内容与生成报告不一致')
    if hashlib.sha256(args.budget.read_bytes()).hexdigest()!=report['delay_budget_sha256']:
        raise ValueError('预算文件与当前GDS不是同一版')
    result = calculate(json.loads(args.budget.read_text(encoding='utf-8')),
                       json.loads(args.mux.read_text(encoding='utf-8')))
    result['gds'] = str(gds)
    result['gds_sha256'] = digest
    result['budget'] = str(args.budget.resolve())
    result['mux_source'] = str(args.mux.resolve())
    args.output_dir.mkdir(parents=True,exist_ok=False)
    (args.output_dir/'全光路时延核算.json').write_text(
        json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    with (args.output_dir/'逐段时延.csv').open('w',encoding='utf-8-sig',newline='') as stream:
        writer = csv.DictWriter(stream,fieldnames=list(result['ledger'][0]))
        writer.writeheader(); writer.writerows(result['ledger'])
    print(json.dumps({k:v for k,v in result.items() if k not in ('ledger','limitations')},
                     ensure_ascii=False,indent=2))


if __name__=='__main__':
    main()
