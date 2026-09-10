"""固定全局LN晶轴、真实平面90°欧拉几何的3D FDTD；先检查内存再求解。"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
import re
from pathlib import Path
import sys
import time
import traceback
import numpy as np
import psutil
from shapely.geometry import LineString,box
from mode_pdk_sweep import DEFAULT_CONFIG,load_config,load_lumapi,zelmon_ln_indices
from crossing_fdtd import add_layer,configure_exact_frequency,assert_wavelength,serial
from bend_mode_analysis import mode_metrics,analyze_folder

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'models/layout'))
from euler_routes import partial_euler_local


def build(f,args,cfg):
    points,_=partial_euler_local(math.pi/2,args.radius_um,.4)
    ex,ey=points[-1];lead=12.
    if args.reference_axis=='Y':
        points=np.array([[-lead,0.],[lead,0.]])
        bounds=(-10.,10.,-8.,8.)
        ports={'input':('x',-6.,0.,'Forward'),'output':('x',6.,0.,'Backward')}
    elif args.reference_axis=='Z':
        points=np.array([[0.,-lead],[0.,lead]])
        bounds=(-8.,8.,-10.,10.)
        ports={'input':('y',-6.,0.,'Forward'),'output':('y',6.,0.,'Backward')}
    else:
        points=np.vstack(([[-lead,0]],points,[[ex,ey+lead]]))
        bounds=(-10.,ex+8.,-8.,ey+10.)
        ports={'input':('x',-6.,0.,'Forward'),'output':('y',ey+6.,ex,'Backward')}
    core=LineString(points).buffer(.35,cap_style=2,join_style=1)
    slab=LineString(points).buffer(3.6,cap_style=2,join_style=1)
    f.switchtolayout();f.deleteall();f.addfdtd()
    xmin,xmax,ymin,ymax=bounds
    for key,value in {'dimension':'3D','x min':xmin*1e-6,'x max':xmax*1e-6,
        'y min':ymin*1e-6,'y max':ymax*1e-6,'z min':-1.5e-6,'z max':3.2e-6,
        'simulation time':args.time_ps*1e-12,'auto shutoff min':1e-6,'mesh accuracy':2}.items():f.set(key,value)
    for axis in 'xyz':
        for side in ['min','max']:f.set(f'{axis} {side} bc','PML')
    f.addrect()
    for key,value in {'name':'SiO2_fill_and_top_cladding','x min':(xmin-3)*1e-6,
        'x max':(xmax+3)*1e-6,'y min':(ymin-3)*1e-6,'y max':(ymax+3)*1e-6,
        'z min':-4e-6,'z max':2.1e-6,'material':cfg['material_models']['sio2'],
        'override mesh order from material database':True,'mesh order':3}.items():f.set(key,value)
    no,ne=[v[0] for v in zelmon_ln_indices(np.array([1.55]))]
    index=f'{no:.12g};{ne:.12g};{no:.12g}'
    add_layer(f,slab,.7,.2,args.slices,70,index,'LN_residual')
    add_layer(f,core,.9,.2,args.slices,70,index,'LN_rib')
    f.addmesh()
    for key,value in {'name':'LN_mesh','x min':xmin*1e-6,'x max':xmax*1e-6,
        'y min':ymin*1e-6,'y max':ymax*1e-6,'z min':.5e-6,'z max':1.3e-6,
        'dx':args.mesh_nm*1e-9,'dy':args.mesh_nm*1e-9,'dz':args.mesh_z_nm*1e-9}.items():f.set(key,value)
    for name,(axis,position,other_position,direction) in ports.items():
        f.addport();other='y' if axis=='x' else 'x'
        for key,value in {'name':name,'injection axis':axis+'-axis','direction':direction,
            axis:position*1e-6,other:other_position*1e-6,other+' span':10e-6,
            'z min':-1e-6,'z max':2.9e-6,'mode selection':'user select'}.items():f.set(key,value)
        f.seteigensolver('number of trial modes',max(12,args.modes))
    f.select('FDTD::ports');f.set('source port','input');f.set('source mode','1')
    frequency=configure_exact_frequency(f,1550.)
    for name in ports:
        f.select('FDTD::ports::'+name)
        f.updateportmodes(np.arange(1,args.modes+1,dtype=float))
    f.addpower()
    for key,value in {'name':'field_xy','monitor type':'2D Z-normal','z':1e-6,
        'x min':(xmin+1)*1e-6,'x max':(xmax-1)*1e-6,
        'y min':(ymin+1)*1e-6,'y max':(ymax-1)*1e-6,
        'down sample x':4,'down sample y':4}.items():f.set(key,value)
    return {'centerline_um':points.tolist(),'core_polygon_um':list(core.exterior.coords),
            'platform_polygon_um':list(slab.exterior.coords),'ports':ports,'frequency':frequency,
            'axes':{'x':'crystal_Y','y':'crystal_Z','z':'crystal_X'},
            'tensor_index_xyz':[no,ne,no],'reference_axis':args.reference_axis,
            'scope':'独立90度欧拉弯或12um端口间直参考；有限7.2um残余平台名义截面'}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--radius-um',type=float,default=80.)
    p.add_argument('--mesh-nm',type=float,default=130.)
    p.add_argument('--mesh-z-nm',type=float,default=50.)
    p.add_argument('--slices',type=int,default=8)
    p.add_argument('--modes',type=int,default=8)
    p.add_argument('--time-ps',type=float,default=3.)
    p.add_argument('--reference-axis',choices=['Y','Z'])
    p.add_argument('--build-only',action='store_true')
    p.add_argument('--source-project',type=Path,help='从已建未完成工程另存后求解；保留源工程')
    p.add_argument('--source-geometry',type=Path,help='源工程对应参数与几何JSON')
    p.add_argument('--memory-headroom-gib',type=float,default=1.)
    a=p.parse_args();a.output_dir.mkdir(parents=True,exist_ok=False)
    state_path=a.output_dir/'状态.json'
    state={'status':'opening','controller_pid':psutil.Process().pid,'validation_pass':False,
           'args':{k:str(v) if isinstance(v,Path) else v for k,v in vars(a).items()}}
    def save(**kw):
        state.update(kw);state['updated_at']=time.strftime('%Y-%m-%d %H:%M:%S')
        state_path.write_text(json.dumps(serial(state),ensure_ascii=False,indent=2),encoding='utf-8')
        print(state['status'],flush=True)
    save()
    cfg=load_config(DEFAULT_CONFIG);api=load_lumapi(cfg)
    try:
        with api.FDTD(hide=True) as f:
            save(status='building',version=f.version())
            if a.source_project:
                if not a.source_geometry or not a.source_project.is_file():
                    raise ValueError('从保留工程运行需要对应几何文件')
                source_digest=hashlib.sha256(a.source_project.read_bytes()).hexdigest()
                f.load(str(a.source_project.resolve()))
                if not f.layoutmode():
                    raise RuntimeError('源工程已有分析结果，应先提取，拒绝重新求解覆盖')
                info=json.loads(a.source_geometry.read_text(encoding='utf-8'))
                save(source_project=str(a.source_project.resolve()),source_sha256=source_digest,
                     numerical_settings_source='从源FSP直接读取；命令行网格/时间参数不改写源模型')
            else:
                info=build(f,a,cfg)
            checkpoint_before={key:f.getnamed('FDTD',key) for key in
                ['checkpoint during simulation','checkpoint at shutoff','checkpoint period']}
            f.setnamed('FDTD','checkpoint during simulation',True)
            f.setnamed('FDTD','checkpoint at shutoff',True)
            port_checks={}
            for name,spec in info['ports'].items():
                obj='FDTD::ports::'+name
                port_checks[name]=mode_metrics(f.getresult(obj,'mode profiles'),f.getresult(obj,'neff'),spec[0])
            source_selection=str(f.getnamed('FDTD::ports','source mode'))
            match=re.fullmatch(r'(?:mode\s+)?(\d+)(?:\.0+)?',source_selection.strip())
            if not match:raise ValueError('不能识别源模式选择：'+source_selection)
            source_mode=int(match.group(1))
            if source_mode!=port_checks['input']['target_mode_number']:
                raise ValueError('实际源模式与TE0身份不一致，拒绝求解')
            info['source_mode_number']=source_mode
            (a.output_dir/'端口模式预检查.json').write_text(
                json.dumps(serial(port_checks),ensure_ascii=False,indent=2),encoding='utf-8')
            wavelength=float(f.getglobalmonitor('frequency center'))
            if not np.isclose(299792458./wavelength*1e9,1550.,rtol=0,atol=1e-6):
                raise ValueError('源工程监视器频率不匹配1550nm')
            save(checkpoint_enabled=True,checkpoint_period_default=checkpoint_before['checkpoint period'],
                 port_mode_identity_verified=True)
            project=(a.output_dir/'真实欧拉弯.fsp').resolve()
            f.save(str(project))
            (a.output_dir/'参数与几何.json').write_text(json.dumps(serial(info),ensure_ascii=False,indent=2),encoding='utf-8')
            save(status='checking_memory',project=str(project))
            check=f.runsystemcheck()
            (a.output_dir/'内存预检.json').write_text(json.dumps(serial(check),ensure_ascii=False,indent=2),encoding='utf-8')
            if 'Memory_Recommended_Bytes' in check:
                required=float(check['Memory_Recommended_Bytes'])
            elif 'Memory_Recommended' in check:
                required=float(check['Memory_Recommended'])
            else:
                # 2023R2没有新版顶层推荐内存键；按运行/初始化较大者增加20%余量。
                memory=check['Approximate_Memory_Requirements']
                required=1.2*max(float(memory['Running_Simulation_Bytes']),
                                 float(memory['Initialization_and_Mesh_Bytes']))
            available=psutil.virtual_memory().available
            save(status='built',recommended_memory_bytes=required,available_memory_bytes=available)
            if a.build_only:return
            headroom=max(a.memory_headroom_gib*1024**3,.2*required)
            if required>available-headroom:
                raise RuntimeError('三维FDTD内存预估超过当前安全余量，工程已保存，未启动求解')
            old={}
            try:
                for key,value in [('processes','1'),('threads','6')]:
                    old[key]=f.getresource('FDTD',1,key);f.setresource('FDTD',1,key,value)
                start=time.time();save(status='running')
                f.run()
            finally:
                for key,value in old.items():f.setresource('FDTD',1,key,value)
            results={}
            for name in info['ports']:
                obj='FDTD::ports::'+name
                s=f.getresult(obj,'S');assert_wavelength(s,1550.,name)
                results[name]=s
                for kind in ['mode profiles','neff','expansion for port monitor']:
                    data=f.getresult(obj,kind);assert_wavelength(data,1550.,name+' '+kind)
                    np.savez_compressed(a.output_dir/(name+'_'+kind.replace(' ','_')+'.npz'),
                                        **{k:v for k,v in data.items() if isinstance(v,np.ndarray)})
            field=f.getresult('field_xy','E');assert_wavelength(field,1550.,'field_xy')
            np.savez_compressed(a.output_dir/'场分布.npz',**{k:v for k,v in field.items() if isinstance(v,np.ndarray)})
            status=f.getresult('FDTD','status');f.save(str(project))
            (a.output_dir/'结果.json').write_text(json.dumps(serial({'ports':results,'fdtd_status':status,
                'elapsed_s':time.time()-start,'validation_pass':False}),ensure_ascii=False,indent=2),encoding='utf-8')
            save(status='single_run_completed',fdtd_status=status,elapsed_s=time.time()-start,
                 output_power_by_mode=abs(np.asarray(results['output']['S']).ravel())**2,
                 reflection_power_by_mode=abs(np.asarray(results['input']['S']).ravel())**2)
            analysis=analyze_folder(a.output_dir)
            (a.output_dir/'模式功率分析.json').write_text(json.dumps(serial(analysis),ensure_ascii=False,indent=2),encoding='utf-8')
            if a.source_project and hashlib.sha256(a.source_project.read_bytes()).hexdigest()!=source_digest:
                raise RuntimeError('源工程哈希发生变化')
    except Exception:
        save(status='failed',error=traceback.format_exc());raise


if __name__=='__main__':main()
