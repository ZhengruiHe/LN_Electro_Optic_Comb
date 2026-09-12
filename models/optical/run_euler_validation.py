"""串行完成欧拉弯名义数值验证算例；保留进度、原始结果，绝不自动宣称全部通过。"""
from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time
import traceback

import numpy as np
import psutil

ROOT=Path(__file__).resolve().parents[2]
WORKER=Path(__file__).with_name('euler_bend_eme.py')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seed-dir',type=Path,required=True)
    parser.add_argument('--seed-pid',type=int,required=True)
    parser.add_argument('--output-dir',type=Path,required=True)
    args=parser.parse_args()
    args.output_dir=args.output_dir.resolve()
    args.output_dir.mkdir(parents=True,exist_ok=False)
    state_path=args.output_dir/'批次状态.json'
    state={'status':'waiting_seed','controller_pid':os.getpid(),'started_at':time.strftime('%Y-%m-%d %H:%M:%S'),
           'seed_dir':str(args.seed_dir.resolve()),'completed':[],'full_validation_pass':False}
    def save(**kw):
        state.update(kw);state['updated_at']=time.strftime('%Y-%m-%d %H:%M:%S')
        state_path.write_text(json.dumps(state,ensure_ascii=False,indent=2),encoding='utf-8')
    save()
    try:
        while True:
            seed=json.loads((args.seed_dir/'状态.json').read_text(encoding='utf-8'))
            if seed['status']=='single_run_completed':break
            if seed['status']=='failed':raise RuntimeError('首轮接口核查失败，请查看首轮状态')
            if not psutil.pid_exists(args.seed_pid):raise RuntimeError('首轮控制进程异常消失，停止后续算例')
            time.sleep(20)
        power=np.asarray(seed['output_power_by_mode'],dtype=float).ravel()
        reflected=np.asarray(seed['reflection_power_by_mode'],dtype=float).ravel()
        if not np.all(np.isfinite(np.r_[power,reflected])):
            raise RuntimeError('首轮端口功率存在非有限数')
        if power.sum()+reflected.sum()>1.01:
            raise RuntimeError('首轮功率超出输入超过1%，需先核对模型/模式基底，未继续批次')
        if not seed.get('port_mode_identity_verified') or power.sum()<1e-6:
            raise RuntimeError('首轮模式身份未验证或透射异常低，先核对接口，不排入长批次')
        # 前6项用于最小半径90度的直参考、模式数/分段/网格/边界检查。
        plans=[
            ('直参考_Y',90.,80.,21,12,320,160,16.,True,0.),
            ('直参考_Z',90.,80.,21,12,320,160,16.,True,90.),
            ('90度_R80_分段21_模式12',90.,80.,21,12,320,160,16.,False,0.),
            ('90度_R80_分段41_模式12',90.,80.,41,12,320,160,16.,False,0.),
            ('90度_R80_分段41_模式20',90.,80.,41,20,320,160,16.,False,0.),
            ('90度_R80_细网格',90.,80.,41,20,400,200,16.,False,0.),
            ('90度_R80_扩大边界',90.,80.,41,20,480,200,19.2,False,0.),
        ]
        budget_path=ROOT/'results/layout/_历史归档/迭代版本_20260912/欧拉回路时延预算_v3_无斜直线/欧拉回路_10GHz方向时延预算.json'
        budget=json.loads(budget_path.read_text(encoding='utf-8'))
        radius_special=budget['route_metadata']['loop1']['start_radius_um']
        for label,radius in [('R80',80.),('R81p67',radius_special),('R270',270.)]:
            plans.extend([(f'180度_{label}_分段41',180.,radius,41,20,400,200,16.,False,0.),
                          (f'180度_{label}_分段81',180.,radius,81,20,400,200,16.,False,0.)])
        save(status='running',planned_cases=len(plans),plans=[list(p) for p in plans])
        for name,angle,radius,segments,modes,mesh_y,mesh_z,span,reference,heading in plans:
            folder=args.output_dir/name
            command=[sys.executable,'-X','utf8',str(WORKER),'--angle-deg',str(angle),
                     '--radius-um',str(radius),'--segments',str(segments),'--modes',str(modes),
                     '--mesh-y',str(mesh_y),'--mesh-z',str(mesh_z),'--y-span-um',str(span),
                     '--start-heading-deg',str(heading),'--processes','1','--threads','6',
                     '--output-dir',str(folder)]
            if reference:command.append('--reference')
            save(current_case=name)
            with (args.output_dir/(name+'.log')).open('w',encoding='utf-8') as stream:
                process=subprocess.Popen(command,cwd=ROOT,stdout=stream,stderr=subprocess.STDOUT)
                save(child_pid=process.pid)
                code=process.wait()
            if code!=0:raise RuntimeError(f'{name}未成功结束，保留错误日志，停止后续计算')
            result=json.loads((folder/'状态.json').read_text(encoding='utf-8'))
            if result['status']!='single_run_completed':raise RuntimeError(f'{name}没有有效计算结果')
            power=np.asarray(result['output_power_by_mode'],dtype=float).ravel()
            back=np.asarray(result['reflection_power_by_mode'],dtype=float).ravel()
            if not np.all(np.isfinite(np.r_[power,back])):raise RuntimeError(f'{name}功率非有限')
            target=result.get('output_target_mode_number',1)-1
            if not result.get('port_mode_identity_verified'):
                raise RuntimeError(f'{name}模式身份未验证')
            if reference and abs(-10*np.log10(max(power[target],1e-30)))>.03:
                raise RuntimeError(f'{name}直参考损耗超过0.03dB，禁止继续弯曲扫描')
            if not reference and power.sum()<1e-6:
                raise RuntimeError(f'{name}透射异常低，先排查模型')
            record={'name':name,'folder':str(folder),'TE0_target_mode':target+1,
                    'TE0_candidate_mode1_power':float(power[target]),
                    'TE0_candidate_insertion_loss_db':float(-10*np.log10(max(power[target],1e-30))),
                    'other_output_modes_power':float(power.sum()-power[target]),
                    'total_reflection_power':float(back.sum()),'sum_power':float(power.sum()+back.sum()),
                    'elapsed_s':result['elapsed_s'],'mode_identity_verified':True}
            state['completed'].append(record)
            save(child_pid=None)
            if record['sum_power']>1.01:
                raise RuntimeError(f'{name}非物理功率超限，停止后续算例，待分析')
        comparisons=[]
        for first,second in [('90度_R80_分段21_模式12','90度_R80_分段41_模式12'),
                             ('90度_R80_分段41_模式12','90度_R80_分段41_模式20'),
                             ('90度_R80_分段41_模式20','90度_R80_细网格'),
                             ('90度_R80_细网格','90度_R80_扩大边界')]+[
                             (f'180度_{r}_分段41',f'180度_{r}_分段81') for r in ['R80','R81p67','R270']]:
            a=next(row for row in state['completed'] if row['name']==first)
            b=next(row for row in state['completed'] if row['name']==second)
            delta=abs(a['TE0_candidate_insertion_loss_db']-b['TE0_candidate_insertion_loss_db'])
            comparisons.append({'first':first,'second':second,'delta_IL_db':delta,
                                'candidate_delta_below_0p01db':delta<.01})
        save(status='numerical_batch_completed_pending_mode_review',comparisons=comparisons,
             current_case=None,full_validation_pass=False)
    except Exception:
        save(status='failed',error=traceback.format_exc(),full_validation_pass=False)
        raise


if __name__=='__main__':
    main()
