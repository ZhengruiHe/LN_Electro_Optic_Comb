"""只读检查候选MUX的端口模式场形和波长；不运行/传播/保存LMS。"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
from mode_pdk_sweep import DEFAULT_CONFIG, load_config, load_lumapi


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--project',type=Path,required=True)
    a=p.parse_args();cfg=load_config(DEFAULT_CONFIG);api=load_lumapi(cfg)
    out={'project':str(a.project),'ports':[],'read_only':True}
    with api.MODE(filename=str(a.project.resolve()),hide=True) as m:
        for port in ('port_1','port_2'):
            row={'port':port}
            try:
                data=m.getresult('EME::Ports::'+port,'mode profiles')
                row['result_keys']=sorted(k for k in data if isinstance(data[k],np.ndarray))
                row['wavelength_nm']=float(np.asarray(data['lambda']).ravel()[0]*1e9)
                y=np.asarray(data['y']).ravel()*1e6
                for i in range(1,9):
                    if f'E{i}' not in data:continue
                    e=np.asarray(data[f'E{i}']).squeeze()
                    if e.ndim<2:continue
                    intensity=np.sum(abs(e)**2,axis=-1) if e.shape[-1]==3 else abs(e)**2
                    intensity=np.asarray(intensity).squeeze()
                    if intensity.ndim==2 and intensity.shape[0]==y.size:
                        intensity=np.sum(intensity,axis=1)
                    if intensity.ndim!=1 or intensity.size!=y.size:
                        row.setdefault('mode_profiles',[]).append({'solver_mode':i,'shape':list(e.shape)})
                        continue
                    total=float(np.sum(intensity)); centroid=float(np.sum(y*intensity)/total)
                    row.setdefault('mode_profiles',[]).append({'solver_mode':i,'shape':list(e.shape),
                        'centroid_y_um':centroid,'peak_y_um':float(y[int(np.argmax(intensity))]),
                        'component_fraction':(np.sum(abs(e)**2,axis=tuple(range(e.ndim-1)))/np.sum(abs(e)**2)).tolist()})
            except Exception as exc:row['error']=str(exc)
            out['ports'].append(row)
    print(json.dumps(out,ensure_ascii=False,indent=2))


if __name__=='__main__':main()
