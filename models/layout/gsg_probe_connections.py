"""按150/120um相邻针距给v16加入两端GSG接口，保留光学层和有源T电极。"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import subprocess

import numpy as np

ROOT=Path(__file__).resolve().parents[2]
SOURCE=ROOT/'results/layout/欧拉弯_时延闭合_无斜直线_工作稿_v16_01/四程10GHz_15mm_欧拉弯_时延闭合_待验证工作稿.gds'
KLAYOUT=Path('C:/Users/PC/AppData/Roaming/KLayout/klayout_app.exe')


def connector(pitch_um=150.,signal_pad_width_um=50.,pad_length_um=120.,transition_um=515.):
    if pitch_um not in (120.,150.):raise ValueError('本轮仅支持用户指定的120/150um针距')
    if not 40<=signal_pad_width_um<=50:raise ValueError('当前光路限制信号焊盘宽40–50um')
    x=np.r_[0.,np.linspace(pad_length_um,pad_length_um+transition_um,129)]
    u=np.clip((x-pad_length_um)/transition_um,0,1)
    s=u*u*(3-2*u)
    width=signal_pad_width_um+(43.-signal_pad_width_um)*s
    ground_inner=np.full_like(x,38.5)
    outer=(pitch_um+40.)+(138.5-pitch_um-40.)*s
    def polygon(lower,upper):
        return np.r_[np.c_[x,lower],np.c_[x[::-1],upper[::-1]]].tolist()
    return {'pitch_um':pitch_um,'pad_length_um':pad_length_um,'transition_length_um':transition_um,
            'total_length_um':float(x[-1]),'signal_pad_width_um':signal_pad_width_um,
            'signal_trunk_width_um':43.,'ground_inner_edge_um':38.5,
            'ground_outer_edge_at_pad_um':pitch_um+40.,'ground_outer_edge_at_trunk_um':138.5,
            'polygons':{'S':polygon(-width/2,width/2),
                        'Gupper':polygon(ground_inner,outer),'Glower':polygon(-outer,-ground_inner)},
            'contact_centers_um':{'Gupper':[25.,pitch_um],'S':[25.,0.],'Glower':[25.,-pitch_um]},
            'contact_landing_size_um':[50.,50.],
            'approach':'输入探针从左向右，输出探针从右向左；G/S/G接触中心共线',
            'impedance_status':'几何接续候选，50ohm/S参数尚未电磁验证'}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,default=SOURCE)
    p.add_argument('--pitch-um',type=float,default=150.,choices=[120.,150.])
    p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--klayout-exe',type=Path,default=KLAYOUT)
    a=p.parse_args();a.source=a.source.resolve();a.output_dir=a.output_dir.resolve()
    if not a.source.is_file():raise FileNotFoundError(a.source)
    if a.output_dir.exists():raise FileExistsError('输出已存在，拒绝覆盖')
    a.output_dir.mkdir(parents=True,exist_ok=False)
    output=a.output_dir/f'四程10GHz_15mm_双端GSG{a.pitch_um:g}um_待验证工作稿.gds'
    report=a.output_dir/'GSG连接与版图检查.json'
    config=connector(a.pitch_um)
    payload={'source':str(a.source),'output':str(output),'report':str(report),
             'top':f'EO4P_10G_15MM_GSG{a.pitch_um:g}_DRAFT',
             'cell':f'GSG{a.pitch_um:g}_LAUNCH_50UM_DRAFT','connector':config,
             'preview_output':str(a.output_dir/'版图预览几何.json')}
    source_digest=hashlib.sha256(a.source.read_bytes()).hexdigest()
    payload['source_sha256']=source_digest
    manifest=a.output_dir/'GSG连接参数.json'
    manifest.write_text(json.dumps(payload,ensure_ascii=False,indent=2),encoding='utf-8')
    run=subprocess.run([str(a.klayout_exe),'-b','-r',str(Path(__file__).with_name('klayout_gsg_connections.py')),
                        '-rd','payload_file='+str(manifest)],capture_output=True,text=True,encoding='utf-8',errors='replace')
    if hashlib.sha256(a.source.read_bytes()).hexdigest()!=source_digest:raise RuntimeError('源GDS变化')
    if run.returncode or not report.exists():raise RuntimeError(run.stdout+'\n'+run.stderr)
    result=json.loads(report.read_text(encoding='utf-8'))
    print(json.dumps({k:v for k,v in result.items() if k not in ('preview','connectivity','inherited_cell_signatures')},ensure_ascii=False,indent=2))
    if not result['gsg_geometry_pass']:raise SystemExit(2)


if __name__=='__main__':main()
