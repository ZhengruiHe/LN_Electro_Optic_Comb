"""无源链路组件的固定全局各向异性2D FDTD；保存复S、端口模与传播场。

这是经截面校准的等效模型，厚度已折合；不重新拟合真实LN材料，
不旋转晶轴，不覆盖源版图/旧仿真，不代表完整矢量或流片验证。
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys
import time
import traceback
import numpy as np
import psutil
from shapely.geometry import LineString, Polygon
from shapely.ops import unary_union
from mode_pdk_sweep import DEFAULT_CONFIG, load_config, load_lumapi
from crossing_fdtd import configure_exact_frequency, assert_wavelength, serial
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'models/layout'))
from euler_routes import partial_euler_local
from route_15mm_left_delays import PathBuilder
from route_five_ports_global import s_offset


def geometry(kind,radius,axis,offset_um=5.875,length_um=80.):
    if kind=='reference':
        pts=np.array([[-12.,0.],[12.,0.]])
        ports={'input':('x',-6.,0.,'Forward'),'output':('x',6.,0.,'Backward')}
        core=LineString(pts).buffer(.35,cap_style=2)
    elif kind=='taper':
        x=np.linspace(0,100,401);u=x/100
        w=.7+.5*(3*u*u-2*u*u*u)
        x=np.r_[-12,x,112];w=np.r_[.7,w,1.2]
        core=Polygon(np.vstack((np.column_stack((x,-w/2)),np.column_stack((x[::-1],w[::-1]/2)))))
        pts=np.array([[-12.,0.],[112.,0.]])
        ports={'input':('x',-6.,0.,'Forward'),'output':('x',106.,0.,'Backward')}
    elif kind=='s_bend':
        path=PathBuilder([[-12.,0.],[0.,0.]],0.)
        s_offset(path,length_um,offset_um)
        path.straight(x=length_um+12)
        pts=np.array(path.points)
        ports={'input':('x',-6.,0.,'Forward'),'output':('x',length_um+6,offset_um,'Backward')}
        core=LineString(pts).buffer(.35,cap_style=2,join_style=1)
    else:
        angle=math.pi/2 if kind=='bend90' else math.pi
        p,_=partial_euler_local(angle,radius,.4)
        end=p[-1]
        extension=end+12*np.array([math.cos(angle),math.sin(angle)])
        pts=np.vstack(([[-12.,0.]],p,[extension]))
        ports={'input':('x',-6.,0.,'Forward')}
        if kind=='bend90':ports['output']=('y',end[1]+6,end[0],'Backward')
        else:ports['output']=('x',end[0]-6,end[1],'Forward')
        core=LineString(pts).buffer(.35,cap_style=2,join_style=1)
    slab=LineString(pts).buffer(3.6,cap_style=2,join_style=1)
    if axis=='Z':
        # 只旋转几何90度，保持材料在全局x=晶体Y、y=晶体Z。
        from shapely.affinity import rotate
        core=rotate(core,90,origin=(0,0));slab=rotate(slab,90,origin=(0,0))
        pts=np.column_stack((-pts[:,1],pts[:,0]))
        updated={}
        for n,(a,v,o,d) in ports.items():
            updated[n]=('y',v,-o,d) if a=='x' else ('x',-v,o,'Forward' if d=='Backward' else 'Backward')
        ports=updated
    # 输入/输出在器件延伸段内，PML比端口再外移4um。
    b=slab.bounds
    bounds=[b[0]-5,b[1]-5,b[2]+5,b[3]+5]
    for normal,position,_,direction in ports.values():
        # 几何延伸到端口外6um，边界只外移4um，确保波导穿过PML。
        offset=0 if normal=='x' else 1
        if direction=='Forward':bounds[offset]=position-4
        else:bounds[offset+2]=position+4
    return core,slab,pts,ports,bounds


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--model',type=Path,required=True)
    ap.add_argument('--output-dir',type=Path,required=True)
    ap.add_argument('--kind',choices=['reference','taper','bend90','bend180','s_bend','near_pair'],required=True)
    ap.add_argument('--radius-um',type=float,default=80.)
    ap.add_argument('--axis',choices=['Y','Z'],default='Y')
    ap.add_argument('--mesh-nm',type=float,default=100.)
    ap.add_argument('--time-ps',type=float,default=5.)
    ap.add_argument('--offset-um',type=float,default=5.875)
    ap.add_argument('--length-um',type=float,default=80.)
    ap.add_argument('--pair-config',type=Path)
    ap.add_argument('--source-port',default='input')
    ap.add_argument('--source-mode',type=int,default=1)
    ap.add_argument('--build-only',action='store_true')
    a=ap.parse_args();a.output_dir.mkdir(parents=True,exist_ok=False)
    model=json.loads(a.model.read_text(encoding='utf-8'))
    if model.get('calibration_target')!='0.7um两晶向的相位折射率和RMS场宽':
        raise ValueError('仅拟合neff的探索模型未通过场形检查，不允许用于路由组件')
    state={'status':'opening','controller_pid':psutil.Process().pid,'full_chain_pass':False,
           'args':{k:str(v) if isinstance(v,Path) else v for k,v in vars(a).items()},
           'model_sha256':hashlib.sha256(a.model.read_bytes()).hexdigest()}
    def save(**kw):
        state.update(kw);state['updated_at']=time.strftime('%Y-%m-%d %H:%M:%S')
        (a.output_dir/'状态.json').write_text(json.dumps(serial(state),ensure_ascii=False,indent=2),encoding='utf-8')
        print(state['status'],flush=True)
    pair=None
    if a.kind=='near_pair':
        if not a.pair_config:raise ValueError('near_pair需要实际路线截取文件')
        pair=json.loads(a.pair_config.read_text(encoding='utf-8'))
        if pair.get('geometry_verified_against_gds'):
            current_hash=hashlib.sha256(Path(pair['source_gds']).read_bytes()).hexdigest()
            if current_hash!=pair['source_gds_sha256']:
                raise ValueError('当前GDS与提取的实际芯层哈希不同，不启动旧几何')
        core=unary_union([Polygon(v) for v in pair['core_polygons_um']])
        slab=unary_union([Polygon(v) for v in pair['slab_polygons_um']])
        pts=None;ports=pair['ports'];bounds=pair['bounds_um']
    else:core,slab,pts,ports,bounds=geometry(a.kind,a.radius_um,a.axis,a.offset_um,a.length_um)
    if not core.is_valid or not slab.is_valid:
        raise ValueError('二维几何无效，不启动求解')
    def pieces(poly):return list(poly.geoms) if hasattr(poly,'geoms') else [poly]
    cfg=load_config(DEFAULT_CONFIG);api=load_lumapi(cfg);save()
    try:
        with api.FDTD(hide=True) as f:
            save(status='building',version=f.version())
            f.addfdtd();xmin,ymin,xmax,ymax=bounds
            for k,v in {'dimension':'2D','x min':xmin*1e-6,'x max':xmax*1e-6,'y min':ymin*1e-6,'y max':ymax*1e-6,
                        'z':0.,'index':model['cladding_index'],'simulation time':a.time_ps*1e-12,
                        'auto shutoff min':1e-6,'mesh accuracy':2}.items():f.set(k,v)
            for name,poly,idx,order in [('Effective_slab',slab,model['slab_effective_index_xyz'],2),
                                       ('Effective_rib',core,model['core_effective_index_xyz'],1)]:
                for i,part in enumerate(pieces(poly)):
                    if part.interiors:raise ValueError('不支持有孔多边形直接导入，拒绝错误填充')
                    f.addpoly()
                    for k,v in {'name':name+'_'+str(i),'vertices':np.array(part.exterior.coords[:-1])*1e-6,'z span':2e-6,
                                'material':'<Object defined dielectric>','index':';'.join(str(n) for n in idx),
                                'override mesh order from material database':True,'mesh order':order}.items():f.set(k,v)
            f.addmesh()
            for k,v in {'name':'plane_mesh','x min':xmin*1e-6,'x max':xmax*1e-6,'y min':ymin*1e-6,'y max':ymax*1e-6,
                        'dx':a.mesh_nm*1e-9,'dy':a.mesh_nm*1e-9,'override z mesh':False}.items():f.set(k,v)
            for name,(axis,pos,other,direction) in ports.items():
                f.addport();cross='y' if axis=='x' else 'x'
                for k,v in {'name':name,'injection axis':axis+'-axis','direction':direction,
                            axis:pos*1e-6,cross:other*1e-6,cross+' span':(pair['port_spans_um'][name] if pair else 10.)*1e-6,
                            'z':0.,'mode selection':'user select'}.items():f.set(k,v)
                f.seteigensolver('number of trial modes',12)
            if a.source_port not in ports:raise ValueError('激励端口不属于该组件')
            f.select('FDTD::ports');f.set('source port',a.source_port);f.set('source mode',str(a.source_mode))
            configure_exact_frequency(f,1550.)
            checks={}
            for name,spec in ports.items():
                f.select('FDTD::ports::'+name);f.updateportmodes(np.arange(1,5,dtype=float))
                p=f.getresult('FDTD::ports::'+name,'mode profiles');n=f.getresult('FDTD::ports::'+name,'neff')
                assert_wavelength(p,1550.,name+' modes');assert_wavelength(n,1550.,name+' neff')
                mode_number=a.source_mode if name==a.source_port else 1
                component=np.sum(abs(p[f'E{mode_number}'].reshape(-1,3))**2,axis=0)
                e_inplane=float((component[0]+component[1])/sum(component))
                checks[name]={'neff':serial(n),'first_mode_inplane_E_fraction':e_inplane}
                np.savez_compressed(a.output_dir/(name+'_mode_profiles.npz'),**{k:v for k,v in p.items() if isinstance(v,np.ndarray)})
                if e_inplane<.9:raise ValueError('源/端口第一模不是选定面内Hz分支，不启动求解')
            f.addpower()
            for k,v in {'name':'field_xy','monitor type':'2D Z-normal','z':0.,'x min':(xmin+2)*1e-6,
                        'x max':(xmax-2)*1e-6,'y min':(ymin+2)*1e-6,'y max':(ymax-2)*1e-6,
                        'down sample x':3,'down sample y':3}.items():f.set(k,v)
            project=(a.output_dir/'二维等效组件.fsp').resolve();f.save(str(project))
            (a.output_dir/'参数与几何.json').write_text(json.dumps(serial({'centerline_um':pts,'pair':pair,'ports':ports,'bounds_um':bounds,
                'core_polygons_um':[list(g.exterior.coords) for g in pieces(core)],
                'slab_polygons_um':[list(g.exterior.coords) for g in pieces(slab)],'port_checks':checks,
                'physical_layout_unchanged':True,'scope':'名义有限平台的等效2D组件'}),ensure_ascii=False,indent=2),encoding='utf-8')
            save(status='checking_memory',project=str(project));check=f.runsystemcheck()
            (a.output_dir/'内存预检.json').write_text(json.dumps(serial(check),ensure_ascii=False,indent=2),encoding='utf-8')
            mem=check['Approximate_Memory_Requirements'];required=1.2*max(float(mem['Running_Simulation_Bytes']),float(mem['Initialization_and_Mesh_Bytes']))
            available=psutil.virtual_memory().available
            save(status='built',recommended_memory_bytes=required,available_memory_bytes=available)
            if a.build_only:return
            if required>available-max(1024**3,.2*required):raise RuntimeError('内存余量不足，保留工程但不启动')
            previous={}
            try:
                for k,v in [('processes','1'),('threads','4')]:previous[k]=f.getresource('FDTD',1,k);f.setresource('FDTD',1,k,v)
                start=time.time();save(status='running');f.run()
            finally:
                for k,v in previous.items():f.setresource('FDTD',1,k,v)
            results={}
            for name in ports:
                obj='FDTD::ports::'+name;s=f.getresult(obj,'S');assert_wavelength(s,1550.,name+' S');results[name]=s
                for label in ['expansion for port monitor','neff']:
                    ds=f.getresult(obj,label);assert_wavelength(ds,1550.,name+' '+label)
                    np.savez_compressed(a.output_dir/(name+'_'+label.replace(' ','_')+'.npz'),**{k:v for k,v in ds.items() if isinstance(v,np.ndarray)})
            field=f.getresult('field_xy','E');assert_wavelength(field,1550.,'field_xy')
            np.savez_compressed(a.output_dir/'场分布.npz',**{k:v for k,v in field.items() if isinstance(v,np.ndarray)})
            result={'ports':results,'fdtd_status':f.getresult('FDTD','status'),'elapsed_s':time.time()-start,
                    'validation_pass':False,'full_chain_pass':False,'scope':'未做直参考归一化的等效2D单组件'}
            (a.output_dir/'结果.json').write_text(json.dumps(serial(result),ensure_ascii=False,indent=2),encoding='utf-8')
            f.save(str(project));save(status='component_completed',fdtd_status=result['fdtd_status'],elapsed_s=result['elapsed_s'])
    except Exception as e:
        save(status='failed',error=str(e),traceback=traceback.format_exc());raise


if __name__=='__main__':main()
