"""核对MODE的LN晶向旋转、端口多模和局部曲率支持；不改已有工程。"""
import argparse
import json
import math
from pathlib import Path
import numpy as np

from passive_directional_ng import build
from mode_pdk_sweep import DEFAULT_CONFIG, load_config, load_lumapi, scalar, mode_localization_metrics
from crossing_fdtd import serial


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,required=True)
    args=parser.parse_args()
    args.output_dir.mkdir(parents=True,exist_ok=False)
    cfg=load_config(DEFAULT_CONFIG)
    api=load_lumapi(cfg)
    result={'status':'running','cases':[],'scope':'局部FDE模式检查；不等于完整欧拉弯传输'}
    path=args.output_dir/'局部模式检查.json'
    with api.MODE(hide=True) as mode:
        result['version']=mode.version()
        for heading,radius in [(0,None),(90,None),(45,None),(0,80),(45,80),(90,80)]:
            print(f'FDE heading={heading} radius={radius}',flush=True)
            build(mode,cfg,1.55,'Y')
            theta=math.radians(heading)
            U=np.array([[math.cos(theta),0,math.sin(theta)],
                        [0,1,0],[-math.sin(theta),0,math.cos(theta)]])
            mode.addgridattribute('matrix transform')
            mode.set('name','LN_crystal_transform')
            mode.set('U',U)
            for name in ['LN_Residual_Slab','LN_Rib']:
                mode.setnamed(name,'grid attribute name','LN_crystal_transform')
            mode.select('FDE')
            mode.set('number of trial modes',12)
            if radius:
                mode.set('bent waveguide',True)
                mode.set('bend radius',radius*1e-6)
                mode.set('bend orientation',0)
            row={'heading_deg':heading,'minimum_radius_um':radius,
                 'matrix_U':U,'modes':[]}
            n=int(mode.findmodes())
            for i in range(1,n+1):
                obj=f'FDE::data::mode{i}'
                neff=scalar(mode.getdata(obj,'neff'))
                te=scalar(mode.getdata(obj,'TE polarization fraction')).real
                loss=scalar(mode.getdata(obj,'loss')).real
                loc=mode_localization_metrics(mode,obj,1.5)
                row['modes'].append({'number':i,'neff':neff,'te_fraction':te,
                                     'loss_db_m':loss,'loss_db_cm':loss/100,**loc})
                if loc['central_energy_fraction']>.5:
                    field={key:mode.getdata(obj,key) for key in ('x','y','Ex','Ey','Ez','Hx','Hy','Hz')}
                    np.savez_compressed(args.output_dir/f'heading{heading}_R{radius}_mode{i}.npz',**field)
            mode.save(str((args.output_dir/f'局部模式_heading{heading}_R{radius}.lms').resolve()))
            result['cases'].append(row)
            path.write_text(json.dumps(serial(result),ensure_ascii=False,indent=2),encoding='utf-8')
            print(json.dumps(serial(row['modes']),ensure_ascii=False),flush=True)
        result['status']='local_mode_probe_completed'
        path.write_text(json.dumps(serial(result),ensure_ascii=False,indent=2),encoding='utf-8')


if __name__=='__main__':
    main()
