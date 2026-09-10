"""先验证双晶向直参考，再串行计算真实90度欧拉弯及网格复核。"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import traceback
import numpy as np
import psutil
from bend_mode_analysis import analyze_folder
from euler_status import inspect

ROOT=Path(__file__).resolve().parents[2]


def validate(folder,reference=False):
    s=json.loads((folder/'状态.json').read_text(encoding='utf-8'))
    if s['status']!='single_run_completed':raise RuntimeError('算例未有效完成：'+str(folder))
    p=np.asarray(s['output_power_by_mode'],dtype=float).ravel()
    r=np.asarray(s['reflection_power_by_mode'],dtype=float).ravel()
    stop=int(np.asarray(s['fdtd_status']).ravel()[0])
    if stop!=2:raise RuntimeError('时间窗耗尽或停止状态异常，不继续后续算例：'+str(folder))
    if not np.all(np.isfinite(np.r_[p,r])) or p.sum()+r.sum()>1.01:
        raise RuntimeError('端口功率有效性检查失败')
    analysis=analyze_folder(folder)
    target=analysis['output_TE0_number']-1
    if reference and (abs(-10*np.log10(max(p[target],1e-30)))>.03 or r.sum()>1e-3):
        raise RuntimeError('直参考插损/反射检查未通过，停止弯曲计算')
    return {'folder':str(folder),'fdtd_status':stop,'TE0_candidate_power':float(p[target]),
            'TE0_candidate_IL_db':float(-10*np.log10(max(p[target],1e-30))),
            'other_mode_power':float(p.sum()-p[target]),'reflection_power':float(r.sum()),
            'elapsed_s':s['elapsed_s'],'mode_identity_verified':True}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--seed-dir',type=Path,required=True)
    p.add_argument('--bend-seed-dir',type=Path,help='等待已有130nm完整弯曲，再补110nm直参考和细网格弯曲')
    p.add_argument('--base-mesh-nm',type=float,default=130.)
    p.add_argument('--fine-mesh-nm',type=float,default=110.)
    p.add_argument('--standalone',action='store_true',help='完整弯曲先建模关闭CAD，再调用独立引擎')
    p.add_argument('--output-dir',type=Path,required=True)
    a=p.parse_args();a.output_dir=a.output_dir.resolve();a.output_dir.mkdir(parents=True,exist_ok=False)
    path=a.output_dir/'批次状态.json'
    state={'status':'waiting_Y_reference','controller_pid':os.getpid(),'completed':[],
           'validation_pass':False,'planned_cases':['Z直参考_130nm','90度欧拉弯_130nm','90度欧拉弯_110nm']}
    def save(**kw):
        state.update(kw);state['updated_at']=time.strftime('%Y-%m-%d %H:%M:%S')
        path.write_text(json.dumps(state,ensure_ascii=False,indent=2),encoding='utf-8')
    save()
    try:
        while True:
            seed=json.loads((a.seed_dir/'状态.json').read_text(encoding='utf-8'))
            if seed['status']=='single_run_completed':break
            if seed['status']=='failed':raise RuntimeError('Y直参考失败')
            if not psutil.pid_exists(seed['controller_pid']):raise RuntimeError('Y直参考控制器异常消失')
            time.sleep(15)
        state['completed'].append(validate(a.seed_dir,True));save(status='running')
        if a.bend_seed_dir:
            save(status='waiting_existing_bend',current_case=str(a.bend_seed_dir),
                 planned_cases=['Y直参考_110nm','Z直参考_110nm','90度欧拉弯_110nm'])
            while True:
                seed=json.loads((a.bend_seed_dir/'状态.json').read_text(encoding='utf-8'))
                if seed['status']=='single_run_completed':break
                if seed['status']=='failed' or not inspect(a.bend_seed_dir)['controller_identity_alive']:
                    raise RuntimeError('130nm弯曲未正常结束，未启动后续网格')
                time.sleep(30)
            state['completed'].append(validate(a.bend_seed_dir))
            plan=[('Y直参考_110nm',110,'Y'),('Z直参考_110nm',110,'Z'),('90度欧拉弯_110nm',110,None)]
        else:
            base=a.base_mesh_nm;fine=a.fine_mesh_nm
            plan=[(f'Z直参考_{base:g}nm',base,'Z'),(f'90度欧拉弯_{base:g}nm',base,None),
                  (f'Y直参考_{fine:g}nm',fine,'Y'),(f'Z直参考_{fine:g}nm',fine,'Z'),
                  (f'90度欧拉弯_{fine:g}nm',fine,None)]
            save(planned_cases=[name for name,_,_ in plan])
        for name,mesh,reference in plan:
            folder=a.output_dir/name
            command=[sys.executable,'-X','utf8',str(Path(__file__).with_name('euler_bend_fdtd.py')),
                     '--mesh-nm',str(mesh),'--time-ps','1' if reference else '3','--output-dir',str(folder)]
            if reference:command+=['--reference-axis',reference]
            use_standalone=not reference and (a.bend_seed_dir or a.standalone)
            if use_standalone:command+=['--build-only']
            save(status='running',current_case=name)
            with (a.output_dir/(name+'.log')).open('w',encoding='utf-8') as log:
                proc=subprocess.Popen(command,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
                save(child_pid=proc.pid);code=proc.wait()
            if code!=0:raise RuntimeError(name+'未成功结束；保留状态/日志，不启动后续仿真')
            if use_standalone:
                solved_folder=a.output_dir/(name+'_独立引擎')
                command=[sys.executable,'-X','utf8',str(Path(__file__).with_name('run_saved_euler_fdtd.py')),
                         '--source-dir',str(folder),'--output-dir',str(solved_folder)]
                with (a.output_dir/(name+'_独立引擎.log')).open('w',encoding='utf-8') as log:
                    proc=subprocess.Popen(command,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
                    save(child_pid=proc.pid);code=proc.wait()
                if code!=0:raise RuntimeError('细网格独立引擎未完成，查看其内存/错误状态')
                folder=solved_folder
            state['completed'].append(validate(folder,bool(reference)));save(child_pid=None)
        if a.bend_seed_dir:
            coarse=validate(a.bend_seed_dir)
            fine=state['completed'][-1]
            save(delta_TE0_IL_db=abs(coarse['TE0_candidate_IL_db']-fine['TE0_candidate_IL_db']),
                 delta_other_mode_power=abs(coarse['other_mode_power']-fine['other_mode_power']))
        else:
            bends=[row for row in state['completed'] if '90度欧拉弯' in Path(row['folder']).name]
            if len(bends)==2:
                save(delta_TE0_IL_db=abs(bends[0]['TE0_candidate_IL_db']-bends[1]['TE0_candidate_IL_db']),
                     delta_other_mode_power=abs(bends[0]['other_mode_power']-bends[1]['other_mode_power']))
        save(status='completed_pending_modal_and_mesh_review',current_case=None)
    except Exception:
        save(status='stopped_for_review',error=traceback.format_exc());raise


if __name__=='__main__':main()
