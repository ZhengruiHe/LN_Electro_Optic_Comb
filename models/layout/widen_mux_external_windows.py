"""在v10上加宽MUX外端公共窗口，覆盖两条0.7um路由；LN1脊及回路不动。"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "results/layout/_历史归档/迭代版本_20260912/直脊平行公共窗口_工作稿_v10/四程10GHz_15mm_直脊平行公共窗口_待验证工作稿.gds"
KLAYOUT = Path(r"C:\Users\PC\AppData\Roaming\KLayout\klayout_app.exe")

JOB = r'''
import hashlib
import json
import pya

def signature(layout,cell,skip_text=False):
    rows=[]
    for li in layout.layer_indices():
        info=layout.get_info(li)
        if skip_text and (info.layer,info.datatype)==(101,0):
            continue
        values=sorted(str(s) for s in cell.shapes(li).each())
        if values:
            rows.append((str(info),values))
    return hashlib.sha256(json.dumps(rows).encode()).hexdigest()

layout=pya.Layout()
layout.read(payload['input'])
top=layout.top_cell()
muxes=[c for c in layout.each_cell() if c.name.startswith('MUX_TE01_A70_TAIL030_PARALLEL_WINDOWS')]
if top is None or len(muxes)!=1 or layout.dbu!=.001 or len(layout.top_cells())!=1:
    raise RuntimeError('v10结构、格点或顶层不匹配')
mux=muxes[0]
instances=[i for i in top.each_inst() if i.cell==mux]
if len(instances)!=4:
    raise RuntimeError('MUX实例数不是四个')
other_cells={c.name:signature(layout,c) for c in layout.each_cell() if c not in (top,mux)}
before_core=pya.Region(top.begin_shapes_rec(layout.layer(20,0))).merged()
before_metal=pya.Region(top.begin_shapes_rec(layout.layer(42,0))).merged()
before_instances=sorted(str(i.dcplx_trans) for i in top.each_inst())

specs={
 (20,1):{'old':pya.Box(-200000,-6000,950000,2100),'patch':pya.Box(750000,2100,950000,5050),
         'base_bounds':[-6.0,2.1],'external_bounds':[-6.0,5.05],'branch_width':7.2},
 (10,2):{'old':pya.Box(-200000,-10550,950000,6650),'patch':pya.Box(750000,6650,950000,10050),
         'base_bounds':[-10.55,6.65],'external_bounds':[-10.55,10.05],'branch_width':17.2},
}
before_local={}
for layer_key,spec in specs.items():
    li=layout.layer(*layer_key)
    shapes=list(mux.shapes(li).each())
    if len(shapes)!=1 or not shapes[0].is_box() or shapes[0].box!=spec['old']:
        raise RuntimeError('v10 MUX窗口不匹配：'+str(layer_key))
    before_local[layer_key]=pya.Region(mux.shapes(li)).merged()
    region=before_local[layer_key]+pya.Region(spec['patch'])
    mux.shapes(li).clear()
    mux.shapes(li).insert(region)

mux.name='MUX_TE01_A70_TAIL030_EXTERNAL_WINDOWS_WIDENED_UNVERIFIED'
top.name='EO4P_10G_15MM_EXTERNAL_WINDOWS_WIDENED_DRAFT'
mux.shapes(layout.layer(101,0)).insert(pya.Text('EXTERNAL_COMMON_WINDOWS_WIDENED_X750_950',pya.Trans(760000,9000)))

core=pya.Region(top.begin_shapes_rec(layout.layer(20,0))).merged()
metal=pya.Region(top.begin_shapes_rec(layout.layer(42,0))).merged()
clad=pya.Region(top.begin_shapes_rec(layout.layer(20,1))).merged()
trench=pya.Region(top.begin_shapes_rec(layout.layer(10,2))).merged()
checks={
 'four_mux_instances':len(instances)==4,
 'other_cells_unchanged':all(signature(layout,layout.cell(n))==h for n,h in other_cells.items()),
 'all_LN1_core_unchanged':(core^before_core).is_empty(),
 'all_M1_unchanged':(metal^before_metal).is_empty(),
 'all_instance_transforms_unchanged':before_instances==sorted(str(i.dcplx_trans) for i in top.each_inst()),
 'no_metal_core_overlap':(metal&core).is_empty(),
}
local_rows={}
for layer_key,spec in specs.items():
    li=layout.layer(*layer_key)
    region=pya.Region(mux.shapes(li)).merged()
    changed=region^before_local[layer_key]
    allowed=pya.Region(spec['patch'])
    row={'changed_only_in_external_200um':(changed-allowed).is_empty(),
         'old_region_preserved':(before_local[layer_key]-region).is_empty(),
         'single_connected_region':region.size()==1,
         'bbox_um':[region.bbox().left*.001,region.bbox().bottom*.001,
                    region.bbox().right*.001,region.bbox().top*.001]}
    for branch,center in [('main',-1.8),('aux',1.45)]:
        width=spec['branch_width']
        # 这里只检查MUX单元内部最后2nm；跨单元接续在下方用顶层递归区域另查。
        test=pya.Region(pya.DBox(949.998,center-width/2,950.000,center+width/2).to_itype(layout.dbu))
        row['external_'+branch+'_full_width_covered']=(test-region).is_empty()
    for key,value in row.items():
        if isinstance(value,bool):
            checks[str(layer_key)+'_'+key]=value
    local_rows[str(layer_key)]=row

joint_rows=[]
for index,inst in enumerate(instances,1):
    row={'index':index,'transform':str(inst.dcplx_trans)}
    for layer_key,spec in specs.items():
        full=pya.Region(top.begin_shapes_rec(layout.layer(*layer_key))).merged()
        for branch,center in [('main',-1.8),('aux',1.45)]:
            width=spec['branch_width']
            test=pya.Region(pya.DPolygon(pya.DBox(949.999,center-width/2,950.001,center+width/2)).transformed(inst.dcplx_trans).to_itype(layout.dbu))
            ok=(test-full).is_empty()
            row[str(layer_key)+'_'+branch]=ok
            checks['join_'+str(index)+'_'+str(layer_key)+'_'+branch]=ok
    joint_rows.append(row)
if not all(checks.values()):
    raise RuntimeError('MUX外端加宽检查失败：'+json.dumps(checks))
etch1=clad-core
markers={'core_width_030':core.width_check(300).size(),'core_space_030':core.space_check(300).size(),
         'etch1_width_030':etch1.width_check(300).size(),'etch1_space_030':etch1.space_check(300).size(),
         'sin_trench_width_020':trench.width_check(200).size(),'sin_trench_space_020':trench.space_check(200).size(),
         'metal_width_2':metal.width_check(2000).size(),'metal_space_3':metal.space_check(3000).size()}
layout.write(payload['output'])
result={'status':'MUX外端公共窗口局部加宽工作稿；非光学或流片签核','top_cell':top.name,'mux_cell':mux.name,
        'instances_updated':4,'widening_x_um':[750,950],'local_windows':local_rows,'joint_checks':joint_rows,
        'checks':checks,'manual_rule_markers':markers,'LN1_core_changed':False,'M1_changed':False,
        'optical_validation_pass':False,'foundry_signoff':False}
with open(payload['report'],'w',encoding='utf-8') as stream:
    json.dump(result,stream,ensure_ascii=False,indent=2)
'''


def main() -> None:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,default=BASE)
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--klayout-exe',type=Path,default=KLAYOUT)
    args=parser.parse_args()
    if not args.input.is_file() or not args.klayout_exe.is_file():
        raise FileNotFoundError('找不到v10输入GDS或KLayout')
    args.output_dir.mkdir(parents=True,exist_ok=False)
    output=(args.output_dir/'四程10GHz_15mm_MUX外端公共窗口加宽_待验证工作稿.gds').resolve()
    report_path=(args.output_dir/'MUX外端窗口加宽检查.json').resolve()
    job_path=(args.output_dir/'KLayout批处理_外端窗口加宽.py').resolve()
    source_hash=hashlib.sha256(args.input.read_bytes()).hexdigest()
    payload={'input':str(args.input.resolve()),'output':str(output),'report':str(report_path)}
    job_path.write_text('payload = '+repr(payload)+'\n'+JOB,encoding='utf-8')
    subprocess.run([str(args.klayout_exe.resolve()),'-b','-r',str(job_path)],check=True,
                   creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0))
    if not output.is_file() or not report_path.is_file():
        raise RuntimeError('KLayout没有产生输出GDS或报告')
    if hashlib.sha256(args.input.read_bytes()).hexdigest()!=source_hash:
        raise RuntimeError('v10源文件发生变化')
    report=json.loads(report_path.read_text(encoding='utf-8'))
    report.update({'working_gds':str(output),'source_gds':str(args.input.resolve()),
                   'source_sha256':source_hash,'source_unchanged':True,
                   'limitations':['局部加宽改变20/1与10/2窗口，需按最终截面核对光学模式和射频加载',
                                  '200um矩形加宽区的边界散射尚未仿真',
                                  'MUX本体辅助脊仍有0.3um法向宽度/格点标记',
                                  '21层、圆弯、10GHz时延、RF焊盘与终端仍暂缓']})
    report_path.write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k not in ('checks','joint_checks')},ensure_ascii=False,indent=2))


if __name__=='__main__':
    main()
