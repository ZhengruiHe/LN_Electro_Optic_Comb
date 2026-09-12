"""对新合并候选的完整输入、四程与输出作时延账目；不启动求解器。"""
import argparse
from collections import Counter
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]


def calculate(route, old, mux):
    active_um=float(route.get('active_taper_um',200.))
    external_um=float(route.get('external_taper_um',100.))
    periods=route.get('loop_periods',[3.,3.5,3.])
    targets=[100.*x for x in periods]
    mode_sequence=['TE0','TE1','TE0','TE1']
    active={'TE0':old['loops'][0]['fixed_active_delay_ps'],'TE1':old['loops'][1]['fixed_active_delay_ps']}
    active_taper={'TE0':old['loops'][0]['fixed_two_active_taper_delay_ps']/2*active_um/200.,
                  'TE1':old['loops'][1]['fixed_two_active_taper_delay_ps']/2*active_um/200.}
    external=old['loops'][0]['fixed_two_external_taper_delay_ps']/4*external_um/100.
    channels=mux['channels_at_center']
    def mux_in(mode):
        key='TE0_反向_s14' if mode=='TE0' else '辅助转TE1_反向_s25'
        return external+channels[key]['group_delay_ps']+active_taper[mode]
    def mux_out(mode):
        key='TE0_正向_s41' if mode=='TE0' else 'TE1转辅助_正向_s52'
        return active_taper[mode]+channels[key]['group_delay_ps']+external
    before=route['delays']['input']['directional_delay_ps']+mux_in('TE0')
    intervals=[]
    for i in range(3):
        total=(active[mode_sequence[i]]+mux_out(mode_sequence[i])+
               route['delays'][f'loop{i+1}']['directional_delay_ps']+mux_in(mode_sequence[i+1]))
        intervals.append({'loop':i+1,'target_ps':targets[i],
                          'total_ps':total,'residual_ps':total-targets[i],
                          'phase_residual_deg_at_10GHz':(total-targets[i])*3.6,
                          'passive_route_ps':route['delays'][f'loop{i+1}']['directional_delay_ps'],
                          'passive_route_length_mm':route['delays'][f'loop{i+1}']['centerline_length_um']/1000})
    after=active['TE1']+mux_out('TE1')+route['delays']['output']['directional_delay_ps']
    crossings=Counter()
    for c in route['crossing_details']:
        crossings.update(c['routes'])
    columns={}
    for c in route['crossing_details']:
        columns.setdefault(round(c['x_um'],3),set()).add(round(c['y_um'],3))
    serial_um=sum(max(ys)-min(ys)-45*(len(ys)-1) for ys in columns.values())
    taper_count=2*len(route['crossing_details'])+2*len(columns)+2
    total=before+sum(r['total_ps'] for r in intervals)+after
    length=4*15000+8*750+8*active_um+8*external_um+sum(r['centerline_length_um'] for r in route['delays'].values())
    return {'status':'工程时延账目；非黑盒/总线真实群时延验证','frequency_GHz':10.,'period_ps':100.,
            'electrode_length_mm':15.,'MUX_active_taper_length_um':active_um,
            'MUX_external_taper_length_um':external_um,'intervals':intervals,
            'input_to_first_active_ps':before,'fourth_active_to_output_ps':after,
            'internal_port_to_port_ps':total,'optical_traversal_centerline_length_mm':length/1000,
            'reference_planes':'两只端面耦合器片内端口；不包括两只黑盒本体',
            'physical_crossings':len(route['crossing_details']),
            'serial_bus_length_um':serial_um,'interface_taper_count':taper_count,
            'crossing_traversals_by_path':dict(crossings),
            'crossing_traversals_total':sum(crossings.values()),
            'limitations':[f'交叉黑盒、{taper_count}段0.7/1.2接口拉锥及{serial_um:g}um串接1.2um直段目前沿用0.7um方向ng占位',
                           'MUX外部接续按旧平均ng缩放作工程估计，未重新场求解',
                           '各向异性弯曲仍用ngY*cos²θ+ngZ*sin²θ积分',
                           '没有计入未知黑盒散射相位，不能把微小算术残差解释为器件实际精度']}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--routes',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    data=calculate(json.loads(a.routes.read_text(encoding='utf-8')),
        json.loads((ROOT/'results/layout/_历史归档/迭代版本_20260912/欧拉回路时延预算_v3_无斜直线/欧拉回路_10GHz方向时延预算.json').read_text(encoding='utf-8')),
        json.loads((ROOT/'results/optical/模式复用器_70度_EME_群时延摘要.json').read_text(encoding='utf-8')))
    with a.output.open('x',encoding='utf-8') as f:json.dump(data,f,ensure_ascii=False,indent=2)
    print(json.dumps(data,ensure_ascii=False,indent=2))


if __name__=='__main__':main()
