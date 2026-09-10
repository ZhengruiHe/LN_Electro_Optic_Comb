"""串行运行整链路2D验证的基础组件批次；不等同于完整链路完成。

这里只补齐直参考、0.7/1.2接口和各半径孤立欧拉弯。近邻联合场、
实际MUX外接续和全链路复S级联属于后续阶段，状态里显式保留。
"""
import argparse
import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path
import psutil

ROOT=Path(__file__).resolve().parents[2]


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--model',type=Path,required=True)
    ap.add_argument('--routes',type=Path,required=True)
    ap.add_argument('--output-dir',type=Path,required=True)
    ap.add_argument('--prior-batch',type=Path,help='只复用该批次已达到能量停止的同模型结果，不覆盖原批次')
    ap.add_argument('--time-multiplier',type=float,default=1.)
    ap.add_argument('--allow-window-screen',action='store_true',help='status=1只保留为初筛并继续后续案例，绝不标为收敛通过')
    a=ap.parse_args();a.output_dir.mkdir(parents=True,exist_ok=False)
    routes=json.loads(a.routes.read_text(encoding='utf-8'))
    cases=[{'name':'直参考Z_100nm','kind':'reference','axis':'Z','time':1.},
           {'name':'拉锥Y_100nm','kind':'taper','axis':'Y','time':3.},
           {'name':'拉锥Z_100nm','kind':'taper','axis':'Z','time':3.},
           {'name':'R80_90度Z_100nm','kind':'bend90','axis':'Z','time':4.}]
    radii=sorted({round(b['minimum_radius_um'],9) for d in routes['delays'].values() for b in d['bends']
                  if abs(b['angle_deg'])==180.})
    cases += [{'name':f'R{r:.6f}_180度Y_100nm','kind':'bend180','axis':'Y','radius':r,'time':max(5.,r*.055)} for r in radii]
    cases.append({'name':'右侧S扇出_80um_5p875um','kind':'s_bend','axis':'Y','length':80.,'offset':5.875,'time':3.})
    # 930um横向S的三个实际偏移量，按路由中已记录的两半欧拉曲线相加。
    offsets=[]
    for name in ('loop1','loop2','input'):
        pair=[b for b in routes['delays'][name]['bends'] if b['label'].startswith('路由S渐移')]
        if len(pair)!=2:raise ValueError('无法确认横向S的两段实际几何')
        offsets.append((name,abs(sum(b['displacement_y_um'] for b in pair))))
    cases += [{'name':n+'_横向S_930um','kind':'s_bend','axis':'Y','length':930.,'offset':v,'time':12.} for n,v in offsets]
    state={'status':'starting','controller_pid':psutil.Process().pid,'child_pid':None,'cases':cases,
           'completed':[],'full_chain_pass':False,'blackbox_crossings':'用户确认理想直通，未验证黑盒自身',
           'pending_after_library':['四组MUX近邻接入联合仿真（不能用孤立弯替代）','当前GDS MUX外接续段',
                                    '15mm有源段及既有MUX复S接入','完整相干链路模式功率与场图']}
    state['allow_window_screen']=a.allow_window_screen
    reused={}
    if a.prior_batch:
        previous=json.loads((a.prior_batch/'批次状态.json').read_text(encoding='utf-8'))
        digest=hashlib.sha256(a.model.read_bytes()).hexdigest()
        for entry in previous['completed']:
            if entry['fdtd_status']!=2:continue
            p=Path(entry['result']);s=json.loads((p.parent/'状态.json').read_text(encoding='utf-8'))
            if s['model_sha256']!=digest:raise ValueError('旧组件等效材料不同，不能复用')
            reused[entry['name']]={**entry,'reused_from_prior':True}
    statefile=a.output_dir/'批次状态.json'
    def save():
        state['updated_at']=time.strftime('%Y-%m-%d %H:%M:%S')
        statefile.write_text(json.dumps(state,ensure_ascii=False,indent=2),encoding='utf-8')
    save()
    for i,case in enumerate(cases):
        if case['name'] in reused:
            state['completed'].append(reused[case['name']]);save();continue
        d=a.output_dir/case['name']
        cmd=[sys.executable,'-X','utf8',str(ROOT/'models/optical/link_2d_component.py'),
             '--model',str(a.model.resolve()),'--output-dir',str(d.resolve()),'--kind',case['kind'],
             '--axis',case['axis'],'--mesh-nm','100','--time-ps',str(case['time']*a.time_multiplier)]
        for key,flag in [('radius','--radius-um'),('length','--length-um'),('offset','--offset-um')]:
            if key in case:cmd.extend([flag,str(case[key])])
        state.update(status='running',current_case=case['name'],current_index=i+1,total=len(cases));save()
        with (a.output_dir/(case['name']+'_控制器输出.log')).open('x',encoding='utf-8') as log:
            child=subprocess.Popen(cmd,stdout=log,stderr=subprocess.STDOUT,cwd=ROOT)
            state['child_pid']=child.pid;save()
            code=child.wait()
        state['child_pid']=None
        if code:
            state.update(status='failed',failed_case=case['name'],returncode=code);save();raise SystemExit(code)
        r=json.loads((d/'结果.json').read_text(encoding='utf-8'))
        state['completed'].append({'name':case['name'],'fdtd_status':r['fdtd_status'],'result':str(d/'结果.json')})
        if r['fdtd_status']!=2 and not a.allow_window_screen:
            state.update(status='needs_time_window_review',reason='未触发能量停止；不将时间窗结束当通过，暂停后续组件');save();return
        print(f'completed {i+1}/{len(cases)} {case["name"]}',flush=True);save()
    state.update(status='library_screening_completed' if any(x['fdtd_status']!=2 for x in state['completed']) else 'library_completed',current_case=None);save()
    print('基础组件批次完成；近邻联合场和整链路级联尚未完成',flush=True)


if __name__=='__main__':main()
