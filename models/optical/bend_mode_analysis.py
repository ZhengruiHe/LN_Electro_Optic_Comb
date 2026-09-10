"""欧拉弯端口模式身份与输出功率分类；FDTD固定全局晶轴，1550nm。"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path
import numpy as np
from crossing_fdtd import assert_wavelength,serial


def mode_metrics(profile,neff_dataset,axis):
    assert_wavelength(profile,1550.,'port mode profiles')
    assert_wavelength(neff_dataset,1550.,'port neff')
    normal='y' if axis=='x' else 'x'
    horizontal=np.asarray(profile[normal]).ravel()
    vertical=np.asarray(profile['z']).ravel()
    midpoint=(horizontal[0]+horizontal[-1])/2
    weight=np.abs(np.gradient(horizontal))[:,None]*np.abs(np.gradient(vertical))[None,:]
    neff=np.asarray(neff_dataset['neff']).ravel()
    e_trans=1 if axis=='x' else 0
    records=[]
    for i,n in enumerate(neff,1):
        E=np.asarray(profile[f'E{i}']).squeeze()
        H=np.asarray(profile[f'H{i}']).squeeze()
        if E.shape!=(len(horizontal),len(vertical),3) or H.shape!=E.shape:
            raise ValueError('端口场维度与坐标不一致')
        intensity=abs(E)**2
        full=float(np.sum(intensity*weight[:,:,None]))
        if full<=0 or not np.isfinite(full):raise ValueError('模式场强度无效')
        te=float(np.sum(intensity[:,:,e_trans]*weight)/
                 np.sum((intensity[:,:,e_trans]+intensity[:,:,2])*weight))
        mask=abs(horizontal-midpoint)<=1.5e-6
        center=float(np.sum(intensity[mask]*weight[mask,:,None])/full)
        lateral=np.sum(intensity*weight[:,:,None],axis=(1,2))
        rms=float(np.sqrt(np.sum((horizontal-midpoint)**2*lateral)/full)*1e6)
        # 这是 |E|²空间局域指标，不是严格色散电磁能量占比。
        records.append({'number':i,'neff':n,'TE_transverse_fraction':te,
            'central_E2_fraction_within_1p5um':center,'rms_lateral_um':rms,
            'max_H_over_max_E':float(np.max(abs(H))/np.max(abs(E))),
            'label':'TM_or_hybrid' if te<.55 and center>.7 else 'platform_or_other'})
    candidates=[r for r in records if r['TE_transverse_fraction']>.8
                and r['central_E2_fraction_within_1p5um']>.9]
    if not candidates:raise ValueError('未找到强局域TE0，不能按模式序号直接判断')
    target=max(candidates,key=lambda r:r['neff'].real)
    target['label']='TE0'
    if not 1e-5<target['max_H_over_max_E']<.1:
        raise ValueError('TE0电磁场比例异常，结果不可用')
    return {'target_mode_number':target['number'],'modes':records,'mode_identity_verified':True}


def analyze_folder(folder):
    folder=Path(folder)
    state=json.loads((folder/'状态.json').read_text(encoding='utf-8'))
    if state['status']!='single_run_completed':raise ValueError('工程尚未正常导出完整结果')
    info=json.loads((folder/'参数与几何.json').read_text(encoding='utf-8'))
    ports={}
    for name,spec in info['ports'].items():
        with np.load(folder/(name+'_mode_profiles.npz')) as p,np.load(folder/(name+'_neff.npz')) as n:
            ports[name]=mode_metrics(p,n,spec[0])
    incoming=ports['input']['target_mode_number']
    # 历史脚本统一激励模式1；若识别到非1，不能后处理伪造为TE0输入。
    excited=info.get('source_mode_number',1)
    if incoming!=excited:raise ValueError('实际输入模式并非识别出的TE0')
    powers=np.asarray(state['output_power_by_mode']).ravel()
    reflect=np.asarray(state['reflection_power_by_mode']).ravel()
    outgoing=ports['output']['target_mode_number']-1
    if len(powers)!=len(ports['output']['modes']):raise ValueError('输出模式数与功率数组不匹配')
    target=float(powers[outgoing]);others=float(powers.sum()-target)
    tm=sum(float(powers[i]) for i,r in enumerate(ports['output']['modes']) if r['label']=='TM_or_hybrid')
    stop=int(np.asarray(state['fdtd_status']).ravel()[0])
    return {'folder':str(folder.resolve()),'wavelength_nm':1550.,'port_modes':ports,
        'input_TE0_number':incoming,'output_TE0_number':outgoing+1,'mode_identity_verified':True,
        'TE0_power':target,'raw_TE0_insertion_loss_db':-10*math.log10(max(target,1e-300)),
        'sum_other_output_modes_power':others,'TM_or_hybrid_output_power':tm,
        'sum_reflection_power':float(reflect.sum()),
        'other_mode_crosstalk_db_relative_TE0':(10*math.log10(max(others,1e-300)/target)
                                               if target>0 else None),
        'reflection_db_relative_input':10*math.log10(max(float(reflect.sum()),1e-300)),
        'sum_output_and_reflection':float(powers.sum()+reflect.sum()),
        'fdtd_status':stop,'energy_autoshutoff_reached':stop==2,'validation_pass':False,
        'scope':'单次网格端口结果；模式身份已核对，尚不代表弯曲网格/PML/时间收敛通过'}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('folder',type=Path);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError('分析文件已存在，请指定新路径')
    result=analyze_folder(a.folder)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(serial(result),ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k!='port_modes'},ensure_ascii=False,indent=2))


if __name__=='__main__':main()
