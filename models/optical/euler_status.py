"""只读核对欧拉弯状态、真实进程、求解日志和输出；不回显进程参数。"""
import argparse
import json
from pathlib import Path
import re
import time
import psutil


def inspect(folder):
    folder=Path(folder)
    data=json.loads((folder/'状态.json').read_text(encoding='utf-8'))
    pid=data.get('controller_pid')
    live=False;children=[]
    try:
        process=psutil.Process(pid)
        # 防止PID复用误判；只读取参数确认身份，不打印参数或许可证信息。
        live=any(any(script in item for script in ('euler_bend_fdtd.py','run_saved_euler_fdtd.py'))
                 for item in process.cmdline())
        if data.get('controller_created'):
            live=live and abs(process.create_time()-data['controller_created'])<1
        if live:
            children=[{'pid':p.pid,'name':p.name(),'cpu_s':round(sum(p.cpu_times()[:2]),2)}
                       for p in process.children(recursive=True)
                       if 'fdtd' in p.name().lower()]
    except (psutil.NoSuchProcess,psutil.AccessDenied,TypeError):pass
    state=data['status']
    if state in ('opening','building','checking_memory','built','running') and not live:
        state='interrupted_or_stale'
    logs=list(folder.glob('*_p0.log'))
    progress=None;checkpoint_failures=0
    if logs:
        log=max(logs,key=lambda p:p.stat().st_mtime)
        # 只提取进度字段，不输出日志头部的许可证上下文。
        text=log.read_text(encoding='utf-8',errors='replace')
        checkpoint_failures=text.count('Checkpointing failed')
        matches=list(re.finditer(r'([\d.]+)% complete\. Elapsed simulation time: ([\deE.+-]+) secs\..*?Auto Shutoff: ([\deE.+-]+)',text))
        if matches:
            match=matches[-1]
            progress={'percent':float(match[1]),'physical_time_ps':float(match[2])*1e12,
                      'energy_ratio':float(match[3]),'log_age_seconds':round(time.time()-log.stat().st_mtime,1)}
    return {'recorded_state':data['status'],'effective_state':state,'controller_pid':pid,
            'controller_identity_alive':live,'solver_processes':children,'progress':progress,
            'results_exported':(folder/'结果.json').is_file(),
            'modal_analysis_exported':(folder/'模式功率分析.json').is_file(),
            'checkpoint_files':[p.name for p in folder.iterdir() if p.name.endswith('.ckpt')],
            'partial_checkpoint_files':[p.name for p in folder.iterdir() if p.name.endswith('.ckpt.part')],
            'checkpoint_failure_count':checkpoint_failures,
            'validation_pass':data.get('validation_pass',False)}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('folder',type=Path)
    a=p.parse_args();print(json.dumps(inspect(a.folder),ensure_ascii=False,indent=2))
