"""2D等效模型的独立宽度验证：真实厚度/侧壁的全矢量截面FDE。"""
import argparse
import json
import time
import traceback
from pathlib import Path
import numpy as np
import psutil
from mode_pdk_sweep import (DEFAULT_CONFIG, load_config, load_lumapi, build_case,
                            scalar, mode_localization_metrics)
from passive_directional_ng import index_string


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--width-um',type=float,default=1.2)
    a=p.parse_args();a.output_dir.mkdir(parents=True,exist_ok=False)
    cfg=load_config(DEFAULT_CONFIG);cfg['fde']['trial_modes']=12
    state={'status':'opening','controller_pid':psutil.Process().pid,'full_chain_pass':False,'results':{}}
    def save():
        state['updated_at']=time.strftime('%Y-%m-%d %H:%M:%S')
        (a.output_dir/'状态.json').write_text(json.dumps(state,ensure_ascii=False,indent=2),encoding='utf-8')
        print(state['status'],flush=True)
    save()
    try:
        with load_lumapi(cfg).MODE(hide=True) as m:
            for axis in ('Y','Z'):
                build_case(m,cfg,a.width_um,.2,False,1.55)
                for name in ('LN_Rib','LN_Residual_Slab'):
                    m.setnamed(name,'index',index_string(1.55,axis))
                project=(a.output_dir/f'{axis}_w{a.width_um:g}_独立截面.lms').resolve()
                m.save(str(project));state.update(status='solving_'+axis);save()
                count=int(m.findmodes());candidates=[]
                for i in range(1,count+1):
                    path=f'FDE::data::mode{i}'
                    te=float(scalar(m.getdata(path,'TE polarization fraction')).real)
                    if te<.55:continue
                    loc=mode_localization_metrics(m,path,1.5)
                    row={'solver_mode':i,'neff':scalar(m.getdata(path,'neff')).real,'te_fraction':te,**loc}
                    candidates.append(row)
                    if loc['central_energy_fraction']>.7:
                        field={k:m.getdata(path,k) for k in ['x','y','Ex','Ey','Ez','neff']}
                        np.savez_compressed(a.output_dir/f'{axis}_mode{i}_全矢量场.npz',**field)
                m.save(str(project))
                state['results'][axis]={'wavelength_nm':float(m.getnamed('FDE','wavelength'))*1e9,
                                        'width_um':a.width_um,'modes':candidates}
                state['status']='completed_'+axis;save()
        state['status']='reference_completed';save()
    except Exception as e:
        state.update(status='failed',error=str(e),traceback=traceback.format_exc());save();raise


if __name__=='__main__':main()
