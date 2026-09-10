"""只读审计MUX外部接续结果：实际GDS、LMS、端口复场和CSV交叉核对。

只写新的核对报告/图，不修改仿真工程，不调用run/emepropagate。
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import numpy as np
from shapely.geometry import Polygon,LineString
from shapely.ops import unary_union
from mode_pdk_sweep import DEFAULT_CONFIG,load_config,load_lumapi
from crossing_fdtd import serial

ROOT=Path(__file__).resolve().parents[2]

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,r):
    with p.open('x',encoding='utf-8') as f:json.dump(serial(r),f,ensure_ascii=False,indent=2)
def geometry_family(objects,family):
    obj=max((o for o in objects if o['family']==family),key=lambda o:o['z_max_um'])
    # 分层多边形是薄层中面宽度；向内补偿得到名义顶边，不能拿不同高度比较。
    offset=(obj['z_max_um']-obj['z_min_um'])/2/np.tan(np.deg2rad(70))
    return Polygon(obj['vertices_um']),offset
def section(g,x,offset=0):
    s=g.intersection(LineString([(x,-20),(x,20)]))
    if s.is_empty:return None
    lo,hi=s.bounds[1]+offset,s.bounds[3]-offset
    return {'center_um':(lo+hi)/2,'width_um':hi-lo,'lower_um':lo,'upper_um':hi}

def snapshot(api,project,output):
    before=digest(project);objects=[]
    print('只读LMS：'+project.name,flush=True)
    with api.MODE(filename=str(project.resolve()),hide=True) as m:
        info={'source_project':str(project),'sha256':before,'version':m.version(),
              'layoutmode':bool(m.layoutmode()),'objects':objects,'EME':{}}
        for family in ('LN2_Residual_Platform','LN1_Main_Rib','LN1_Auxiliary_Rib'):
            for i in range(1,33):
                name=f'{family}_slice_{i}'
                if not m.getnamednumber(name):break
                v=np.array(m.getnamed(name,'vertices'))*1e6
                v+=np.array([m.getnamed(name,'x'),m.getnamed(name,'y')],dtype=float).ravel()*1e6
                objects.append({'family':family,'name':name,'vertices_um':v,
                                'z_min_um':float(m.getnamed(name,'z min'))*1e6,
                                'z_max_um':float(m.getnamed(name,'z max'))*1e6,
                                'index':m.getnamed(name,'index')})
        for k in ('wavelength','x min','y span','z min','z max','cells','group spans','energy conservation',
                  'number of modes for all cell groups'):
            info['EME'][k]=m.getnamed('EME',k)
        try:info['available_results']=m.getresult('EME')
        except Exception as e:info['available_results_error']=str(e)
    info['unchanged_after_read']=digest(project)==before
    if not info['unchanged_after_read']:raise RuntimeError('LMS源文件发生变化')
    save(output,info);return info

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--directory',type=Path,required=True)
    ap.add_argument('--gds',type=Path,required=True)
    ap.add_argument('--output-dir',type=Path,required=True)
    ap.add_argument('--reuse-snapshots',type=Path)
    a=ap.parse_args();a.output_dir.mkdir(parents=True,exist_ok=False)
    native=a.output_dir/'实际MUX绘制层.json'
    cmd=['C:/Users/PC/AppData/Roaming/KLayout/klayout_app.exe','-b','-r',str(ROOT/'models/layout/klayout_extract_mux_validation.py'),
         '-rd','gds_path='+str(a.gds.resolve()),'-rd','mux_name=MUX_TE01_A70_A100_E50_BODY750_DRAFT',
         '-rd','report_path='+str(native.resolve())]
    sub=subprocess.run(cmd,capture_output=True)
    if sub.returncode:raise RuntimeError(sub.stderr.decode('utf-8',errors='replace'))
    gds=json.loads(native.read_text(encoding='utf-8'))
    geometries={tuple(l['layer']):unary_union([Polygon(p['exterior'],p['holes']) for p in l['polygons']]) for l in gds['layers']}
    baseline=ROOT/'results/optical/模式复用器_70度_EME_精细检查.lms'
    candidate=a.directory/'MUX外部接续_100um_50um_EME.lms'
    if a.reuse_snapshots:
        old=json.loads((a.reuse_snapshots/'原MUX_LMS快照.json').read_text(encoding='utf-8'))
        new=json.loads((a.reuse_snapshots/'外接续_LMS快照.json').read_text(encoding='utf-8'))
        if digest(baseline)!=old['sha256'] or digest(candidate)!=new['sha256']:raise ValueError('快照哈希不匹配')
    else:
        api=load_lumapi(load_config(DEFAULT_CONFIG))
        old=snapshot(api,baseline,a.output_dir/'原MUX_LMS快照.json')
        new=snapshot(api,candidate,a.output_dir/'外接续_LMS快照.json')
    comparisons=[]
    for x in (-99.99,-50.,.01,187.5,375.,562.5,749.99,775.,799.99):
        parts=geometries[(20,0)].intersection(LineString([(x,-10),(x,10)]))
        pieces=list(parts.geoms) if hasattr(parts,'geoms') else [parts]
        physical=sorted((section(p,x) for p in pieces if not p.is_empty),key=lambda r:r['center_um'])
        entry={'x_um':x,'actual_GDS_ribs':physical,'external_simulation':{},'baseline_body':{}}
        for family in ('LN1_Main_Rib','LN1_Auxiliary_Rib'):
            poly,off=geometry_family(new['objects'],family);entry['external_simulation'][family]=section(poly,x,off)
            if 0<=x<=750:
                poly,off=geometry_family(old['objects'],family);entry['baseline_body'][family]=section(poly,x,off)
        comparisons.append(entry)
    save(a.output_dir/'截面坐标逐点比较.json',comparisons)
    body=[]
    for family in ('LN1_Main_Rib','LN1_Auxiliary_Rib'):
        pa,oa=geometry_family(old['objects'],family);pb,ob=geometry_family(new['objects'],family)
        samples=[]
        for x in np.linspace(.01,749.99,1001):
            sa,sb=section(pa,x,oa),section(pb,x,ob)
            samples.append([x,sb['width_um']-sa['width_um'],sb['center_um']-sa['center_um']])
        v=np.array(samples)
        body.append({'family':family,'max_width_difference_um':float(max(abs(v[:,1]))),
                     'max_center_difference_um':float(max(abs(v[:,2])))})
    profiles=[]
    for name in ('port_1','port_2'):
        with np.load(a.directory/(name+'_mode_profiles.npz')) as p:
            y=p['y'].ravel()*1e6;z=p['z'].ravel()*1e6;w=np.gradient(y)[:,None]*np.gradient(z)[None,:]
            if not np.allclose(p['lambda']*1e9,1550,rtol=0,atol=1e-6):raise ValueError('源模频点不是1550nm')
            for i in range(1,5):
                e=p[f'E{i}'].squeeze();power=abs(e)**2;comp=np.sum(power*w[:,:,None],axis=(0,1));comp/=sum(comp)
                lat=np.sum(power*w[:,:,None],axis=(1,2))
                profiles.append({'port':name,'full_mode':i,'electric_component_fraction_xyz':comp,
                                 'centroid_y_um':float(sum(lat*y)/sum(lat)),
                                 'peak_y_um':float(y[lat.argmax()]),
                                 'label':'vertical_polarization_dominant' if comp[2]>.5 else 'inplane_or_hybrid'})
    with (a.directory/'MUX外部接续_S矩阵.csv').open(encoding='utf-8-sig') as f:rows=list(csv.DictReader(f))
    s=np.zeros((8,8),complex)
    for r in rows:s[int(r['输出索引'])-1,int(r['输入索引'])-1]=complex(float(r['实部']),float(r['虚部']))
    power=abs(s)**2
    summary={'scope':'旧MUX外接续仿真结果审计；不启动求解、不改GDS/LMS',
        'current_layout_validated':False,'body_geometry_matches_baseline':all(v['max_width_difference_um']<.002 and v['max_center_difference_um']<.002 for v in body),
        'body_geometry_difference':body,'profile_wavelength_nm':1550.,'port_profiles':profiles,
        'first_input_raw_power':{'output_mode1':power[4,0],'output_mode2':power[5,0],
            'output_first_two_sum':sum(power[4:6,0]),'all_four_output_sum':sum(power[4:,0]),
            'all_four_reflection_sum':sum(power[:4,0])},
        'EME_settings':new['EME'],'source_files_unchanged':True,
        'findings':['外接续脚本改变了原750um核心的宽度变化规律，不能只归因于端口漏选。',
                    '实际GDS有源侧主路中心为-1.95um，外侧为-1.8/1.45um；旧仿真添加了额外横向漂移。',
                    '输出模式3以竖直电场为主，不能按旧s52编号称TE1转辅助TE0。',
                    '旧仿真残余LN平台固定居中7.2um；不能把新20/1绘制窗口直接当成21层残余LN边界。',
                    'state中的fdtd_status=2是脚本硬编码；EME没有FDTD能量自动停止含义。',
                    '37个EME单元、每单元8模式、每端4模式没有附独立模式数/分段收敛证据。'],
        'recommended_use':'这些结果仅用于排查旧建模流程，不用于判定v13_22外接续通过或失败。'}
    save(a.output_dir/'审计结论.json',summary)
    print(json.dumps(serial({k:v for k,v in summary.items() if k not in ('port_profiles','EME_settings')}),ensure_ascii=False,indent=2))


if __name__=='__main__':main()
