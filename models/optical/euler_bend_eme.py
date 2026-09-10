"""1550nm独立欧拉弯的多模EME；弧长展开并逐段变换固定LN晶轴。

参考Ansys曲率单元EME示例。输出单次计算结果，收敛判断由多网格比较完成。
不修改MUX、v16版图或现有时延预算，不运行制造容差。
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import time
import traceback

import numpy as np

from mode_pdk_sweep import DEFAULT_CONFIG, load_config, load_lumapi
from mode_mux_eme import add_layered_trapezoid, ln_index_string_eme, extract_s_matrix, write_s_matrix
from crossing_fdtd import serial


def identify_port(mode,port,modes):
    """在耗时传播求解前核对真实TE0与E/H一致性，禁止默认为mode1。"""
    obj='EME::Ports::'+port
    data=mode.getresult(obj,'mode profiles')
    neff=np.asarray(mode.getresult(obj,'neff')['neff']).ravel()
    y=np.asarray(data['y']).ravel();z=np.asarray(data['z']).ravel()
    weight=np.gradient(y)[:,None]*np.gradient(z)[None,:]
    rows=[]
    for i in range(1,modes+1):
        E=np.asarray(data[f'E{i}']).squeeze();H=np.asarray(data[f'H{i}']).squeeze()
        intensity=abs(E)**2;norm=np.sum(intensity*weight[:,:,None])
        te=np.sum(intensity[:,:,1]*weight)/np.sum((intensity[:,:,1]+intensity[:,:,2])*weight)
        central=np.sum(intensity[abs(y)<1.5e-6]*weight[abs(y)<1.5e-6,:,None])/norm
        ratio=float(np.max(abs(H))/np.max(abs(E)))
        rows.append({'number':i,'neff':neff[i-1],'te_fraction':float(te),
                     'central_fraction':float(central),'max_H_over_max_E':ratio})
    eligible=[r for r in rows if r['te_fraction']>.8 and r['central_fraction']>.9]
    if not eligible:
        raise RuntimeError(port+'没有唯一强局域TE0候选，拒绝按mode1继续')
    target=max(eligible,key=lambda r:r['neff'].real)
    if not 1e-5<target['max_H_over_max_E']<.1:
        raise RuntimeError(port+'的TE0电磁场比例异常：'+str(target['max_H_over_max_E']))
    return {'target_mode':target['number'],'modes':rows,'electromagnetic_field_sanity_pass':True}


def profile(radius, angle_deg, fraction, segments, lead_um, reference=False, start_heading_deg=0.):
    angle=math.radians(abs(angle_deg)); sign=1 if angle_deg>=0 else -1
    ramp=radius*fraction*angle
    circle=radius*(1-fraction)*angle
    total=2*ramp+circle
    counts=[max(2,round(segments*fraction/(1+fraction))),
            max(2,round(segments*(1-fraction)/(1+fraction)))]
    spans=np.r_[lead_um,np.full(counts[0],ramp/counts[0]),
                np.full(counts[1],circle/counts[1]),np.full(counts[0],ramp/counts[0]),lead_um]
    def state(s):
        if s<=0:return 0.,0.
        if s>=total:return sign*angle,0.
        if s<ramp:return sign*s*s/(2*radius*ramp),sign*s/(radius*ramp)
        if s<ramp+circle:return sign*(ramp/(2*radius)+(s-ramp)/radius),sign/radius
        u=s-ramp-circle
        return sign*(ramp/(2*radius)+circle/radius+u/radius-u*u/(2*radius*ramp)),sign*(1-u/ramp)/radius
    positions=np.cumsum(spans)-spans/2
    states=[state(s-lead_um) for s in positions]
    headings=np.array([s[0] for s in states]); curvature=np.array([s[1] for s in states])
    if reference:
        spans=np.array([lead_um,total,lead_um])
        positions=np.cumsum(spans)-spans/2
        headings=np.zeros(3);curvature=np.zeros(3)
    headings+=math.radians(start_heading_deg)
    if not reference and not np.isclose(np.sum(curvature*spans),sign*angle,rtol=0,atol=1e-12):
        raise ValueError('曲率积分与目标转角不闭合')
    return {'spans_um':spans,'midpoints_um':positions,'heading_rad':headings,
            'curvature_per_um':curvature,'arc_length_um':total,'total_length_um':float(sum(spans))}


def build(mode,cfg,args):
    pr=profile(args.radius_um,args.angle_deg,.4,args.segments,5.,args.reference,args.start_heading_deg)
    spans=pr['spans_um']; count=len(spans); length=sum(spans)
    mode.switchtolayout(); mode.deleteall()
    # EME局部坐标：x切向，y平面内法向，z竖直。
    # 物理晶轴固定；只将其在局部坐标系中的表示变换phi=-heading。
    index=ln_index_string_eme(1.55)
    if args.tensor_representation=='legacy_profile':
        mode.putv('arc_coordinate_m',pr['midpoints_um']*1e-6)
        mode.putv('phi_profile_rad',-pr['heading_rad'])
        mode.eval('PR=rectilineardataset("LN local orientation",arc_coordinate_m,0,0);'
                  'PR.addattribute("phi",phi_profile_rad);'
                  'PR.addattribute("theta",0*phi_profile_rad);'
                  'PR.addattribute("psi",0*phi_profile_rad);'
                  'addgridattribute("permittivity rotation",PR);'
                  'set("name","LN_fixed_crystal_in_local_frame");')
        sections=[(-1.,length+1.,'LN_fixed_crystal_in_local_frame')]
    else:
        sections=[];edges=np.r_[0,np.cumsum(spans)]
        # 使用每段常值、严格正交的矩阵，避免空间旋转数据的接口异常。
        # 全局晶体没有旋转；U把局部坐标向量变回原主轴坐标。
        for i,theta in enumerate(pr['heading_rad']):
            c,s=math.cos(theta),math.sin(theta)
            U=np.array([[c,-s,0],[s,c,0],[0,0,1.]])
            U[abs(U)<1e-14]=0
            if not np.allclose(U.T@U,np.eye(3),atol=1e-12):
                raise ValueError('晶向变换矩阵不正交')
            attribute=f'LN_crystal_U_{i+1}'
            mode.addgridattribute('matrix transform');mode.set('name',attribute);mode.set('U',U)
            sections.append((float(edges[i])-(1 if i==0 else 0),
                             float(edges[i+1])+(1 if i==count-1 else 0),attribute))
    for section,(left,right,attribute) in enumerate(sections,1):
        for family,width,z in [('LN2_residual',7.2,.7),('LN1_rib',.7,.9)]:
            name=f'{family}_group{section}'
            add_layered_trapezoid(mode,name,np.array([left,right]),np.zeros(2),np.full(2,width),
                                  z,.2,70.,args.slices,index)
            for i in range(args.slices):
                mode.setnamed(f'{name}_slice_{i+1}','grid attribute name',attribute)
    mode.addrect()
    for key,value in {'name':'Oxide_fill_and_top_cladding','x min':-1e-6,
        'x max':(length+1)*1e-6,'y':0.,'y span':(args.y_span_um+4)*1e-6,
        'z min':-4e-6,'z max':2.1e-6,'material':cfg['material_models']['sio2'],
        'override mesh order from material database':True,'mesh order':3}.items():
        mode.set(key,value)
    mode.addeme()
    for key,value in {'solver type':'3D','x min':0.,'y':0.,'y span':args.y_span_um*1e-6,
        'z min':-2.5e-6,'z max':3.2e-6,'wavelength':1.55e-6,
        'background index':1.,'number of cell groups':count,'group spans':spans*1e-6,
        'cells':np.ones(count),'subcell method':np.zeros(count),
        'number of modes for all cell groups':args.modes,
        'mesh cells y':args.mesh_y,'mesh cells z':args.mesh_z,
        'convergence tolerance':1e-9,'energy conservation':'none',
        'allow custom eigensolver settings':True}.items():
        mode.set(key,value)
    for side in ['y min','y max','z min','z max']:
        mode.set(side+' bc','PML')
    for i,k in enumerate(pr['curvature_per_um'],1):
        mode.select(f'EME::Cells::cell_{i}')
        mode.seteigensolver('bent waveguide',bool(abs(k)>1e-14))
        if abs(k)>1e-14:
            mode.seteigensolver('bend radius',1/abs(k)*1e-6)
            # 正曲率圆心在+y(局部法向)，Ansys对应180度。
            mode.seteigensolver('bend orientation',180. if k>0 else 0.)
    for port in ['port_1','port_2']:
        mode.select('EME::Ports::'+port)
        mode.set('use full simulation span',True)
        mode.set('mode selection','user select')
        if int(mode.updateportmodes(np.arange(1,args.modes+1,dtype=float)))!=1:
            raise RuntimeError('EME端口模式更新失败：'+port)
    return pr


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--radius-um',type=float,default=80.)
    p.add_argument('--angle-deg',type=float,default=180.)
    p.add_argument('--start-heading-deg',type=float,default=0.)
    p.add_argument('--tensor-representation',choices=['cell_matrix','legacy_profile'],default='cell_matrix')
    p.add_argument('--segments',type=int,default=21)
    p.add_argument('--modes',type=int,default=12)
    p.add_argument('--mesh-y',type=int,default=320)
    p.add_argument('--mesh-z',type=int,default=160)
    p.add_argument('--y-span-um',type=float,default=16.)
    p.add_argument('--slices',type=int,default=8)
    p.add_argument('--processes',type=int,default=1)
    p.add_argument('--threads',type=int,default=6)
    p.add_argument('--reference',action='store_true')
    p.add_argument('--output-dir',type=Path,required=True)
    args=p.parse_args()
    if args.radius_um<80 or args.modes<4 or args.segments<5:
        raise ValueError('半径或求解离散参数无效')
    args.output_dir.mkdir(parents=True,exist_ok=False)
    state_path=args.output_dir/'状态.json'
    state={'status':'opening','args':vars(args),'scope':'单次欧拉弯EME；尚未完成模式数/分段/网格收敛'}
    def update(status,**kw):
        state.update(status=status,**kw)
        state_path.write_text(json.dumps(serial({**state,'args':{k:str(v) if isinstance(v,Path) else v
            for k,v in vars(args).items()}}),ensure_ascii=False,indent=2),encoding='utf-8')
        print(status,flush=True)
    update('opening'); start=time.time()
    cfg=load_config(DEFAULT_CONFIG);api=load_lumapi(cfg)
    try:
        with api.MODE(hide=True) as mode:
            update('building',version=mode.version())
            pr=build(mode,cfg,args)
            port_checks={port:identify_port(mode,port,args.modes) for port in ['port_1','port_2']}
            (args.output_dir/'端口身份与电磁场检查.json').write_text(
                json.dumps(serial(port_checks),ensure_ascii=False,indent=2),encoding='utf-8')
            (args.output_dir/'曲率与晶向.json').write_text(json.dumps(serial(pr),indent=2),encoding='utf-8')
            project=(args.output_dir/'欧拉弯多模.lms').resolve()
            mode.save(str(project))
            update('solving_modes',project=str(project),groups=len(pr['spans_um']))
            old={}
            try:
                for key,value in [('processes',args.processes),('threads',args.threads)]:
                    old[key]=mode.getresource('EME',1,key)
                    mode.setresource('EME',1,key,str(value))
                mode.run()
            finally:
                for key,value in old.items():
                    mode.setresource('EME',1,key,value)
            update('propagating')
            mode.emepropagate()
            raw=extract_s_matrix(mode.getresult('EME','user s matrix'))
            normalized=extract_s_matrix(mode.getresult('EME','power normalized user s matrix'))
            if normalized.shape!=(2*args.modes,2*args.modes):
                raise RuntimeError(f'端口S矩阵维数异常：{normalized.shape}')
            write_s_matrix(args.output_dir/'原始S矩阵.csv',raw)
            write_s_matrix(args.output_dir/'功率归一化S矩阵.csv',normalized)
            np.savez_compressed(args.output_dir/'S矩阵.npz',raw=raw,power_normalized=normalized)
            # 端口模式先保存完整信息，不能把总透射直接当成目标TE0透射。
            result_names={}
            for obj in ['EME::Ports::port_1','EME::Ports::port_2','EME::Cells::cell_1',
                        'EME::Cells::cell_'+str(len(pr['spans_um']))]:
                try:result_names[obj]=str(mode.getresult(obj))
                except Exception as exc:result_names[obj]=str(exc)
                if 'mode fields' in result_names[obj]:
                    field=mode.getresult(obj,'mode fields')
                    arrays={k:v for k,v in field.items() if isinstance(v,np.ndarray)}
                    np.savez_compressed(args.output_dir/(obj.split('::')[-1]+'_mode_fields.npz'),**arrays)
            (args.output_dir/'可用模式结果.json').write_text(json.dumps(result_names,ensure_ascii=False,indent=2),encoding='utf-8')
            mode.save(str(project))
            input_index=port_checks['port_1']['target_mode']-1
            output_index=port_checks['port_2']['target_mode']-1
            power=abs(normalized[:,input_index])**2
            update('single_run_completed',elapsed_s=time.time()-start,
                   matrix_shape=list(normalized.shape),input_mode_number=input_index+1,
                   output_target_mode_number=output_index+1,port_mode_identity_verified=True,
                   target_TE0_power=float(power[args.modes+output_index]),
                   output_power_by_mode=power[args.modes:],reflection_power_by_mode=power[:args.modes],
                   sum_output_and_reflected_power=float(sum(power)),
                   raw_sum_output_and_reflected_power=float(sum(abs(raw[:,input_index])**2)),
                   full_bend_validation_pass=False)
    except Exception:
        update('failed',error=traceback.format_exc())
        raise


if __name__=='__main__':
    main()
