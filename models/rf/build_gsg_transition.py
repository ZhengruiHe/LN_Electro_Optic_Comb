"""从v17实际GDS导出的M1轮廓建立单侧GSG接入HFSS工程；默认仅建模。"""
from __future__ import annotations
import argparse
import copy
import hashlib
import json
from pathlib import Path
import time
import traceback

import psutil
from shapely.geometry import Polygon,box

from build_hfss_cpw import load_config,DEFAULT_CONFIG,installed_aedt_roots,create_stack_and_cpw,add_solution

ROOT=Path(__file__).resolve().parents[2]


def extract(directory):
    report=json.loads((directory/'GSG连接与版图检查.json').read_text(encoding='utf-8'))
    if not report['gsg_geometry_pass']:raise ValueError('输入GSG几何未通过检查')
    gds=Path(report['output_gds'])
    if hashlib.sha256(gds.read_bytes()).hexdigest()!=report['output_sha256']:
        raise ValueError('GDS哈希已改变')
    layers=json.loads((directory/'版图预览几何.json').read_text(encoding='utf-8'))
    metal=next(layer for layer in layers if layer['layer']==[42,0])
    x0=2000-report['connection_parameters']['total_length_um']
    x1=2500. # 保留原T形电极前10周期，末端取周期边界。
    window=box(x0,-300,x1,300)
    polygons={}
    for row in metal['polygons']:
        poly=Polygon(row['exterior'],row['holes']).intersection(window)
        if poly.is_empty:continue
        if poly.geom_type!='Polygon' or poly.interiors:raise ValueError('截取M1不是无孔单多边形')
        if poly.bounds[1]<0<poly.bounds[3]:name='signal'
        elif poly.bounds[1]>0:name='ground_right'
        else:name='ground_left'
        if name in polygons:raise ValueError('M1导体分裂')
        # HFSS x=横向晶体Z，y=传播晶体Y；仅变换坐标，不旋转晶体取向约定。
        polygons[name]=[[float(y),float(x-x0)] for x,y in list(poly.exterior.coords)[:-1]]
    if set(polygons)!={'signal','ground_left','ground_right'}:raise ValueError('三导体提取不完整')
    return report,polygons,x1-x0


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--layout-dir',type=Path,required=True)
    p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--solve',action='store_true')
    a=p.parse_args();a.layout_dir=a.layout_dir.resolve();a.output_dir=a.output_dir.resolve()
    if not str(a.output_dir).isascii():
        raise ValueError('AEDT2023 gRPC保存路径需使用ASCII目录；中文说明仍由Python正常写入')
    report,polygons,length=extract(a.layout_dir)
    if a.output_dir.exists():raise FileExistsError('输出目录已存在，拒绝覆盖')
    a.output_dir.mkdir(parents=True,exist_ok=False)
    state={'status':'preflight','controller_pid':psutil.Process().pid,'geometry_source_gds':report['output_gds'],
           'gds_sha256':report['output_sha256'],'pitch_um':report['pitch_um'],
           'total_length_um':length,'active_T_stub_um':500.,'frequency_ghz':10.,'rf_performance_pass':False}
    path=a.output_dir/'GSG射频工程状态.json'
    def save(**kw):
        state.update(kw);state['updated_at']=time.strftime('%Y-%m-%d %H:%M:%S')
        path.write_text(json.dumps(state,ensure_ascii=False,indent=2),encoding='utf-8')
        print(state['status'],flush=True)
    save()
    cfg=load_config(DEFAULT_CONFIG)
    cfg['aedt']['frequency_start_ghz']=8.;cfg['aedt']['frequency_stop_ghz']=12.
    profile=copy.deepcopy(cfg['simulation_profiles']['screening'])
    profile.update(adaptive_frequencies_ghz=[10.],frequency_stop_ghz=12.,sweep_points=41)
    model={'metal_polygons_HFSS_xy_um':polygons,'config':cfg,'simulation_profile':profile,
           'scope':'实际M1接入635um+有源T电极500um；两个平面波端口，不含实际探针针头',
           'port_planes_GDS_x_um':[2000-report['connection_parameters']['total_length_um'],2500.],
           'limitations':['介质层沿用已有射频材料模型；系数未代工标定',
                          'LN部分为两条1.33um脊+7.2um有限平台的参考截面，未纳入MUX局部变化或远处回路',
                          'SiN移除区沿用原射频模型的全宽SiO2填充近似',
                          '不含真实探针针头、接触电阻、表面压痕和实际背面载台',
                          '本轮只准备几何及端口；建模验证不等于S参数通过']}
    (a.output_dir/'GSG射频模型参数.json').write_text(json.dumps(model,ensure_ascii=False,indent=2),encoding='utf-8')
    if not installed_aedt_roots():raise RuntimeError('未找到已安装HFSS')
    if psutil.virtual_memory().available<2*1024**3:raise RuntimeError('可用内存不足以安全打开AEDT建模')
    from ansys.aedt.core import Hfss
    project=a.output_dir/f'GSG{report["pitch_um"]:g}_Transition_10GHz.aedt'
    hfss=None
    try:
        save(status='opening_AEDT')
        hfss=Hfss(project=str(project),design='GSG_Launch_Actual_M1',solution_type='Terminal',
                  version='2023.1',non_graphical=True,new_desktop=True,close_on_exit=True)
        save(status='building')
        # regular g17只用于复用层叠/参考LN位置；端口前会以实际GDS金属整体替换。
        solids=create_stack_and_cpw(hfss,cfg,length,43.,17.,100.,'regular',profile,
                                   metal_polygons_xy=polygons)
        setup=add_solution(hfss,cfg,profile)
        validation=hfss.validate_simple(str(a.output_dir/'hfss_validation.log'))
        saved=hfss.save_project(str(project))
        if not saved or not project.is_file():raise RuntimeError('AEDT工程未实际保存')
        save(status='model_built' if validation==1 else 'model_validation_failed',
             project=str(project),geometry_validation_code=validation,
             solid_names=[obj.name for obj in solids.values()],limitations=model['limitations'])
        if validation!=1:raise RuntimeError('HFSS几何/端口检查未通过')
        if a.solve:
            if any('fdtd-engine' in proc.info['name'].lower() for proc in psutil.process_iter(['name'])):
                raise RuntimeError('FDTD正在运行，本脚本不并发启动HFSS求解')
            save(status='solving')
            if not hfss.analyze(setup=setup,cores=4):raise RuntimeError('HFSS求解失败')
            hfss.export_touchstone(setup=setup,sweep='Sweep_RF',output_file=str(a.output_dir/'GSG过渡.s4p'),
                                   renormalization=False,gamma_impedance_comments=True)
            hfss.save_project(str(project));save(status='solved_pending_review')
    except Exception:
        save(status='failed',error=traceback.format_exc());raise
    finally:
        if hfss:hfss.release_desktop(close_projects=True,close_desktop=True)


if __name__=='__main__':main()
