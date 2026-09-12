"""HFSS screening for the straight O-style GSG150 input launch.
This is a local launch model, not a full 15 mm tapeout signoff.
"""
from __future__ import annotations
import argparse, copy, json, hashlib
from pathlib import Path
from ansys.aedt.core import Hfss
from build_hfss_cpw import (
    DEFAULT_CONFIG, load_config, create_stack_and_cpw, add_solution,
    installed_aedt_roots,
)

ROOT=Path(__file__).resolve().parents[2]
p=argparse.ArgumentParser()
p.add_argument('--output-dir',type=Path,required=True)
p.add_argument('--non-graphical',action='store_true')
p.add_argument('--cores',type=int,default=4)
a=p.parse_args(); a.output_dir=a.output_dir.resolve(); a.output_dir.mkdir(parents=True,exist_ok=False)
if not installed_aedt_roots(): raise RuntimeError('未找到Ansys Electronics Desktop')
cfg=load_config(DEFAULT_CONFIG)
cfg['aedt']['frequency_start_ghz']=8.; cfg['aedt']['frequency_stop_ghz']=12.
profile=copy.deepcopy(cfg['simulation_profiles']['screening'])
profile.update({'adaptive_frequencies_ghz':[10.0],'frequency_stop_ghz':12.0,'sweep_points':81})
# GDS O型输入的等效局部坐标：321um输入过渡 + 500um均匀主干。
L=821.0
x=[0.,120.,220.,321.,821.]
sw=[25.,25.,21.5,21.5,21.5]
si=[38.5]*5
so=[190.,190.,138.5,138.5,138.5]
def poly(left,right):
    # HFSS x=横向，y=传播；left/right 已分别给出左右边界坐标。
    return [[left[i],x[i]] for i in range(len(x))]+[[right[i],x[i]] for i in range(len(x)-1,-1,-1)]
polygons={
    'signal':poly([-sw[i] for i in range(5)],sw),
    'ground_left':poly([-so[i] for i in range(5)],[-si[i] for i in range(5)]),
    'ground_right':poly(si,so),
}
model={
 'status':'HFSS_O_STYLE_GSG150_INPUT_SCREENING',
 'line_length_um':L,
 'input_transition_um':321.,
 'signal_pad_width_um':50.,
 'signal_trunk_width_um':43.,
 'gsg_pitch_um':150.,
 'electrode_style':'regular_local_launch',
 'geometry_source':'v6c O型输入候选的参数化等效几何；未复制黑盒内部',
 'materials':cfg['materials'],
 'limitations':['仅局部输入过渡和821um均匀CPW','未含15mm周期T帽','未含HTR终端与接触寄生','LN RF张量/M1电导率仍为占位值','未含实际探针和背面载台']
}
(a.output_dir/'HFSS模型参数.json').write_text(json.dumps(model,ensure_ascii=False,indent=2),encoding='utf-8')
project=a.output_dir/'O_input_GSG150_HFSS_screening.aedt'
state={'status':'opening','project':str(project),'output_dir':str(a.output_dir),'frequency_ghz':[8.,12.],
       'geometry_sha256':hashlib.sha256(json.dumps(polygons,sort_keys=True).encode()).hexdigest(),
       'materials':{'LN':'Project_LN_RF_Anisotropic','M1':'Project_M1_RF','HTR':'not_in_local_launch_model'}}
(a.output_dir/'HFSS状态.json').write_text(json.dumps(state,ensure_ascii=False,indent=2),encoding='utf-8')
hfss=None
try:
    hfss=Hfss(project=str(project),design='O_INPUT_GSG150',solution_type='Terminal',
              version='2023.1',non_graphical=a.non_graphical,new_desktop=True,close_on_exit=True)
    solids=create_stack_and_cpw(hfss,cfg,L,43.,5.,100.,'regular',profile,metal_polygons_xy=polygons)
    setup=add_solution(hfss,cfg,profile)
    log=a.output_dir/'HFSS几何验证.log'
    code=hfss.validate_simple(str(log))
    # PyAEDT 1.4 在部分 AEDT 版本上验证成功时返回 None；只有明确的非1/None才阻断。
    if code not in (None,1): raise RuntimeError(f'HFSS几何验证失败: {code}')
    hfss.save_project(str(project))
    state.update(status='solving',geometry_validation_code=code)
    (a.output_dir/'HFSS状态.json').write_text(json.dumps(state,ensure_ascii=False,indent=2),encoding='utf-8')
    if not hfss.analyze(setup=setup,cores=a.cores): raise RuntimeError('HFSS求解失败')
    touch=a.output_dir/'O_input_GSG150_8_12GHz.s4p'
    hfss.export_touchstone(setup=setup,sweep='Sweep_RF',output_file=str(touch),renormalization=False,gamma_impedance_comments=True)
    hfss.save_project(str(project))
    state.update(status='solved_pending_review',touchstone=str(touch),geometry_validation_code=code)
    (a.output_dir/'HFSS状态.json').write_text(json.dumps(state,ensure_ascii=False,indent=2),encoding='utf-8')
finally:
    if hfss: hfss.release_desktop(close_projects=True,close_desktop=True)
print(json.dumps(state,ensure_ascii=False,indent=2))
