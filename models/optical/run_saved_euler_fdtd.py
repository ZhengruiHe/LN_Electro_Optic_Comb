"""独立MPI引擎运行已预检的欧拉FSP；关闭CAD后启动，保留检查点与状态。"""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import time
import traceback

import numpy as np
import psutil
from mode_pdk_sweep import DEFAULT_CONFIG,load_config,load_lumapi
from crossing_fdtd import assert_wavelength,serial
from bend_mode_analysis import analyze_folder

MPI=Path('C:/Program Files/Microsoft MPI/Bin/mpiexec.exe')
ENGINE=Path('C:/Program Files/Lumerical/v232/bin/fdtd-engine-msmpi.exe')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source-dir',type=Path,required=True)
    p.add_argument('--output-dir',type=Path,required=True)
    a=p.parse_args();a.source_dir=a.source_dir.resolve();a.output_dir=a.output_dir.resolve()
    source_state=json.loads((a.source_dir/'状态.json').read_text(encoding='utf-8'))
    if not source_state.get('checkpoint_enabled') or not source_state.get('port_mode_identity_verified'):
        raise ValueError('源工程必须完成端口身份和检查点设置')
    source=a.source_dir/'真实欧拉弯.fsp'
    memory=json.loads((a.source_dir/'内存预检.json').read_text(encoding='utf-8'))
    required=float(memory.get('Memory_Recommended',memory['Approximate_Memory_Requirements']['Running_Simulation_Bytes']))
    if not MPI.is_file() or not ENGINE.is_file():raise FileNotFoundError('已安装MPI/引擎未找到')
    if a.output_dir.exists():raise FileExistsError('输出目录已存在，拒绝覆盖/重复求解')
    a.output_dir.mkdir(parents=True,exist_ok=False)
    state={'status':'preflight','controller_pid':os.getpid(),'controller_created':psutil.Process().create_time(),
           'source_project':str(source),'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
           'checkpoint_enabled':True,'port_mode_identity_verified':True,'validation_pass':False,
           'recommended_memory_bytes':required,'available_memory_bytes':psutil.virtual_memory().available,
           'resources':{'processes':1,'threads':6,'MPI':'Microsoft MPI','engine_version':'v232'}}
    def save(**kw):
        state.update(kw);state['updated_at']=time.strftime('%Y-%m-%d %H:%M:%S')
        (a.output_dir/'状态.json').write_text(json.dumps(serial(state),ensure_ascii=False,indent=2),encoding='utf-8')
    save()
    try:
        if psutil.virtual_memory().available<required+max(1024**3,.2*required):
            raise RuntimeError('独立求解可用内存仍不足，未启动')
        for name in ['真实欧拉弯.fsp','参数与几何.json','内存预检.json','端口模式预检查.json']:
            shutil.copy2(a.source_dir/name,a.output_dir/name)
        project=a.output_dir/'真实欧拉弯.fsp'
        source_check=hashlib.sha256(source.read_bytes()).hexdigest()
        if hashlib.sha256(project.read_bytes()).hexdigest()!=source_check:raise ValueError('工程复制不一致')
        info=json.loads((a.output_dir/'参数与几何.json').read_text(encoding='utf-8'))
        start=time.time()
        with (a.output_dir/'独立引擎输出.log').open('w',encoding='utf-8') as log:
            proc=subprocess.Popen([str(MPI),'-n','1',str(ENGINE),'-t','6',str(project)],
                                   stdout=log,stderr=subprocess.STDOUT,cwd=a.output_dir,
                                   creationflags=subprocess.CREATE_NO_WINDOW)
            save(status='running',project=str(project),mpi_pid=proc.pid)
            while proc.poll() is None:
                try:
                    children=psutil.Process(proc.pid).children(recursive=True)
                    save(elapsed_s=time.time()-start,
                         engine_pids=[c.pid for c in children if 'fdtd-engine' in c.name()])
                except psutil.NoSuchProcess:pass
                time.sleep(15)
        if proc.returncode!=0:raise RuntimeError(f'引擎退出码{proc.returncode}，保留日志和检查点')
        save(status='collecting',engine_exit_code=proc.returncode)
        api=load_lumapi(load_config(DEFAULT_CONFIG))
        with api.FDTD(filename=str(project),hide=True) as f:
            status=f.getresult('FDTD','status')
            if int(np.asarray(status).ravel()[0]) not in [1,2]:
                raise RuntimeError('求解器没有有效完成结果')
            results={}
            for name in info['ports']:
                obj='FDTD::ports::'+name
                s=f.getresult(obj,'S');assert_wavelength(s,1550.,name);results[name]=s
                for kind in ['mode profiles','neff','expansion for port monitor']:
                    data=f.getresult(obj,kind);assert_wavelength(data,1550.,name+' '+kind)
                    np.savez_compressed(a.output_dir/(name+'_'+kind.replace(' ','_')+'.npz'),
                                        **{k:v for k,v in data.items() if isinstance(v,np.ndarray)})
            data=f.getresult('field_xy','E');assert_wavelength(data,1550.,'field_xy')
            np.savez_compressed(a.output_dir/'场分布.npz',**{k:v for k,v in data.items() if isinstance(v,np.ndarray)})
        (a.output_dir/'结果.json').write_text(json.dumps(serial({'ports':results,
            'fdtd_status':status,'elapsed_s':time.time()-start,'validation_pass':False}),
            ensure_ascii=False,indent=2),encoding='utf-8')
        save(status='single_run_completed',fdtd_status=status,elapsed_s=time.time()-start,
             output_power_by_mode=abs(np.asarray(results['output']['S']).ravel())**2,
             reflection_power_by_mode=abs(np.asarray(results['input']['S']).ravel())**2)
        analysis=analyze_folder(a.output_dir)
        (a.output_dir/'模式功率分析.json').write_text(json.dumps(serial(analysis),ensure_ascii=False,indent=2),encoding='utf-8')
        try:
            from plot_euler_fdtd import plot
            save(field_plot=plot(a.output_dir,a.output_dir/'传播场与输出模式功率.png'))
        except Exception:
            save(plot_warning=traceback.format_exc())
        if hashlib.sha256(source.read_bytes()).hexdigest()!=source_check:raise RuntimeError('源工程变化')
        print(json.dumps({k:v for k,v in analysis.items() if k!='port_modes'},ensure_ascii=False))
    except Exception:
        save(status='failed',error=traceback.format_exc());raise


if __name__=='__main__':main()
