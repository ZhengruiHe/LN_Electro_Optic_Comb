"""用相同X-normal截面对照LN张量表示，检查磁场和功率，不修改既有工程。"""
import argparse
import json
import math
from pathlib import Path
import numpy as np
from mode_pdk_sweep import DEFAULT_CONFIG,load_config,load_lumapi,zelmon_ln_indices,scalar
from mode_mux_eme import add_layered_trapezoid
from crossing_fdtd import serial


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output-dir',type=Path,required=True)
    a=p.parse_args();a.output_dir.mkdir(parents=True,exist_ok=False)
    cfg=load_config(DEFAULT_CONFIG);api=load_lumapi(cfg)
    no,ne=[v[0] for v in zelmon_ln_indices(np.array([1.55]))]
    result=[]
    with api.MODE(hide=True) as mode:
        for name,angle,representation in [('Y_diagonal',0,'diagonal'),('Z_diagonal',90,'diagonal'),
            ('Z_matrix',90,'matrix'),('45_matrix',45,'matrix')]:
            print(name,flush=True)
            mode.switchtolayout();mode.deleteall()
            if representation=='diagonal':
                index=f'{no};{ne};{no}' if angle==0 else f'{ne};{no};{no}'
            else:
                index=f'{no};{ne};{no}'
                t=math.radians(angle);c=math.cos(t);s=math.sin(t)
                U=np.array([[c,-s,0],[s,c,0],[0,0,1.]])
                U[abs(U)<1e-14]=0
                mode.addgridattribute('matrix transform');mode.set('name','LN_rotation');mode.set('U',U)
            for family,width,z in [('slab',7.2,.7),('rib',.7,.9)]:
                add_layered_trapezoid(mode,family,np.array([-1.,1.]),np.zeros(2),np.full(2,width),
                                      z,.2,70.,8,index)
                if representation!='diagonal':
                    for i in range(8):mode.setnamed(f'{family}_slice_{i+1}','grid attribute name','LN_rotation')
            mode.addrect()
            for k,v in {'name':'oxide','x span':2e-6,'y span':20e-6,'z min':-4e-6,'z max':2.1e-6,
                        'material':cfg['material_models']['sio2'],'override mesh order from material database':True,
                        'mesh order':3}.items():mode.set(k,v)
            mode.addfde()
            for k,v in {'solver type':'2D X normal','x':0.,'y':0.,'y span':16e-6,
                        'z min':-2.5e-6,'z max':3.2e-6,'mesh cells y':320,'mesh cells z':160,
                        'wavelength':1.55e-6,'number of trial modes':4}.items():mode.set(k,v)
            for side in ['y min','y max','z min','z max']:mode.set(side+' bc','PML')
            count=int(mode.findmodes());rows=[]
            for i in range(1,count+1):
                obj=f'FDE::data::mode{i}';yy=np.asarray(mode.getdata(obj,'y')).ravel()
                zz=np.asarray(mode.getdata(obj,'z')).ravel();w=np.gradient(yy)[:,None]*np.gradient(zz)[None,:]
                E=np.stack([mode.getdata(obj,k).squeeze() for k in ['Ex','Ey','Ez']],axis=-1)
                H=np.stack([mode.getdata(obj,k).squeeze() for k in ['Hx','Hy','Hz']],axis=-1)
                I=abs(E)**2
                row={'mode':i,'neff':scalar(mode.getdata(obj,'neff')),
                     'maxE':np.max(abs(E),axis=(0,1)),'maxH':np.max(abs(H),axis=(0,1)),
                     'TEfraction':float(np.sum(I[:,:,1]*w)/np.sum((I[:,:,1]+I[:,:,2])*w)),
                     'central_fraction':float(np.sum(I[abs(yy)<1.5e-6]*w[abs(yy)<1.5e-6,:,None])/np.sum(I*w[:,:,None]))}
                rows.append(row)
            result.append({'case':name,'modes':rows})
            mode.save(str((a.output_dir/(name+'.lms')).resolve()))
            (a.output_dir/'张量对照.json').write_text(json.dumps(serial(result),ensure_ascii=False,indent=2),encoding='utf-8')
            print(json.dumps(serial(rows),ensure_ascii=False),flush=True)


if __name__=='__main__':main()
