"""读取已求解的欧拉弯EME工程，导出输入/输出/中间单元模式，不重新求解。"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from mode_pdk_sweep import DEFAULT_CONFIG,load_config,load_lumapi
from crossing_fdtd import serial


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project',type=Path,required=True)
    parser.add_argument('--groups',type=int,required=True)
    parser.add_argument('--output-dir',type=Path,required=True)
    args=parser.parse_args()
    args.output_dir.mkdir(parents=True,exist_ok=False)
    digest=hashlib.sha256(args.project.read_bytes()).hexdigest()
    api=load_lumapi(load_config(DEFAULT_CONFIG))
    result={'project':str(args.project.resolve()),'sha256':digest,'objects':{}}
    with api.MODE(filename=str(args.project.resolve()),hide=True) as mode:
        result['wavelength_nm']=float(mode.getnamed('EME','wavelength'))*1e9
        if not np.isclose(result['wavelength_nm'],1550,rtol=0,atol=1e-6):
            raise ValueError('工程波长不是1550nm')
        for obj in ['EME::Ports::port_1','EME::Ports::port_2',
                    'EME::Cells::cell_1',f'EME::Cells::cell_{args.groups//2}',
                    f'EME::Cells::cell_{args.groups}']:
            names=str(mode.getresult(obj))
            row={'available_results':names}
            for name in ['mode fields','mode profiles','neff']:
                if name not in names:
                    continue
                dataset=mode.getresult(obj,name)
                if isinstance(dataset,dict):
                    arrays={key:value for key,value in dataset.items() if isinstance(value,np.ndarray)}
                    row[name]={key:{'shape':list(value.shape),'complex':bool(np.iscomplexobj(value))}
                               for key,value in arrays.items()}
                    filename=obj.split('::')[-1]+'_'+name.replace(' ','_')+'.npz'
                    np.savez_compressed(args.output_dir/filename,**arrays)
                    for key,value in arrays.items():
                        if value.size<=50:
                            row[name][key]['values']=serial(value)
            result['objects'][obj]=row
    if hashlib.sha256(args.project.read_bytes()).hexdigest()!=digest:
        raise RuntimeError('源工程发生变化')
    (args.output_dir/'模式数据索引.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))


if __name__=='__main__':
    main()
