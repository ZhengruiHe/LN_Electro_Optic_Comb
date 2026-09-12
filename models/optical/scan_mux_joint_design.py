"""小规模MUX联合参数筛选：间隙、渐变轨迹、耦合长度和辅助脊末端宽度。

每个候选单独运行EME，结果只写入新目录。该批次是名义筛选，不是容差扫描；
S矩阵的物理通道仍需结合端口场形确认，不能只按模式编号下结论。
"""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]


def read_s(path: Path):
    with path.open(encoding='utf-8-sig', newline='') as stream:
        return {(int(r['输出索引']), int(r['输入索引'])): float(r['功率'])
                for r in csv.DictReader(stream)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--preset', choices=('initial', 'deep', 'fixed_length', 'refined', 'main_width'), default='initial')
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=False)
    candidates = [
        {'id':'C01', 'minimum_gap_um':.65, 'aux_start_um':.35, 'aux_stop_um':.55,
         'progress_exponent':2., 'progress_ratio':2.3333333333, 'length_um':550.},
        {'id':'C02', 'minimum_gap_um':.80, 'aux_start_um':.35, 'aux_stop_um':.55,
         'progress_exponent':2., 'progress_ratio':2.3333333333, 'length_um':650.},
        {'id':'C03', 'minimum_gap_um':.65, 'aux_start_um':.35, 'aux_stop_um':.60,
         'progress_exponent':2., 'progress_ratio':2.3333333333, 'length_um':650.},
        {'id':'C04', 'minimum_gap_um':.65, 'aux_start_um':.35, 'aux_stop_um':.50,
         'progress_exponent':1.5, 'progress_ratio':2., 'length_um':750.},
        {'id':'C05', 'minimum_gap_um':.80, 'aux_start_um':.35, 'aux_stop_um':.50,
         'progress_exponent':1.5, 'progress_ratio':2., 'length_um':750.},
        {'id':'C06', 'minimum_gap_um':.65, 'aux_start_um':.40, 'aux_stop_um':.60,
         'progress_exponent':2., 'progress_ratio':1.5, 'length_um':650.},
    ] if args.preset == 'initial' else [
        {'id':'D01', 'minimum_gap_um':.65, 'aux_start_um':.35, 'aux_stop_um':.40,
         'progress_exponent':2., 'progress_ratio':5., 'length_um':1000.},
        {'id':'D02', 'minimum_gap_um':.65, 'aux_start_um':.35, 'aux_stop_um':.45,
         'progress_exponent':2., 'progress_ratio':5., 'length_um':1200.},
        {'id':'D03', 'minimum_gap_um':.80, 'aux_start_um':.35, 'aux_stop_um':.40,
         'progress_exponent':2., 'progress_ratio':5., 'length_um':1500.},
        {'id':'D04', 'minimum_gap_um':.65, 'aux_start_um':.35, 'aux_stop_um':.50,
         'progress_exponent':3., 'progress_ratio':5., 'length_um':1500.},
    ] if args.preset == 'deep' else [
        {'id':'F01', 'minimum_gap_um':.65, 'aux_start_um':.35, 'aux_stop_um':.45,
         'progress_exponent':2., 'progress_ratio':2.3333333333, 'length_um':750.},
        {'id':'F02', 'minimum_gap_um':.65, 'aux_start_um':.35, 'aux_stop_um':.50,
         'progress_exponent':2., 'progress_ratio':5., 'length_um':750.},
        {'id':'F03', 'minimum_gap_um':.65, 'aux_start_um':.35, 'aux_stop_um':.55,
         'progress_exponent':3., 'progress_ratio':5., 'length_um':750.},
        {'id':'F04', 'minimum_gap_um':.80, 'aux_start_um':.35, 'aux_stop_um':.45,
         'progress_exponent':2., 'progress_ratio':5., 'length_um':750.},
        {'id':'F05', 'minimum_gap_um':.80, 'aux_start_um':.35, 'aux_stop_um':.50,
         'progress_exponent':3., 'progress_ratio':5., 'length_um':750.},
        {'id':'F06', 'minimum_gap_um':.65, 'aux_start_um':.40, 'aux_stop_um':.55,
         'progress_exponent':3., 'progress_ratio':5., 'length_um':750.},
        {'id':'F07', 'minimum_gap_um':.65, 'aux_start_um':.35, 'aux_stop_um':.45,
         'progress_exponent':2., 'progress_ratio':5., 'trajectory':'slow_crossing', 'length_um':750.},
    ] if args.preset == 'fixed_length' else [
        {'id':'R01', 'minimum_gap_um':.75, 'aux_start_um':.35, 'aux_stop_um':.45,
         'progress_exponent':2., 'progress_ratio':5., 'length_um':750.},
        {'id':'R02', 'minimum_gap_um':.85, 'aux_start_um':.35, 'aux_stop_um':.45,
         'progress_exponent':2., 'progress_ratio':5., 'length_um':750.},
        {'id':'R03', 'minimum_gap_um':.80, 'aux_start_um':.35, 'aux_stop_um':.42,
         'progress_exponent':2., 'progress_ratio':5., 'length_um':750.},
        {'id':'R04', 'minimum_gap_um':.80, 'aux_start_um':.35, 'aux_stop_um':.48,
         'progress_exponent':2., 'progress_ratio':5., 'length_um':750.},
        {'id':'R05', 'minimum_gap_um':.80, 'aux_start_um':.35, 'aux_stop_um':.45,
         'progress_exponent':1.8, 'progress_ratio':5., 'length_um':750.},
        {'id':'R06', 'minimum_gap_um':.80, 'aux_start_um':.35, 'aux_stop_um':.45,
         'progress_exponent':2.2, 'progress_ratio':5., 'length_um':750.},
        {'id':'R07', 'minimum_gap_um':.80, 'aux_start_um':.35, 'aux_stop_um':.45,
         'progress_exponent':2., 'progress_ratio':4., 'length_um':750.},
        {'id':'R08', 'minimum_gap_um':.80, 'aux_start_um':.35, 'aux_stop_um':.45,
         'progress_exponent':2., 'progress_ratio':6., 'length_um':750.},
    ] if args.preset == 'refined' else [
        {'id':'M01', 'minimum_gap_um':.80, 'aux_start_um':.35, 'aux_stop_um':.42,
         'main_start_um':1.33, 'main_stop_um':1.10, 'progress_exponent':2., 'progress_ratio':5., 'length_um':750.},
        {'id':'M02', 'minimum_gap_um':.80, 'aux_start_um':.35, 'aux_stop_um':.42,
         'main_start_um':1.35, 'main_stop_um':1.08, 'progress_exponent':2., 'progress_ratio':5., 'length_um':750.},
        {'id':'M03', 'minimum_gap_um':.80, 'aux_start_um':.35, 'aux_stop_um':.42,
         'main_start_um':1.40, 'main_stop_um':1.10, 'progress_exponent':2., 'progress_ratio':5., 'length_um':750.},
        {'id':'M04', 'minimum_gap_um':.80, 'aux_start_um':.35, 'aux_stop_um':.42,
         'main_start_um':1.45, 'main_stop_um':1.12, 'progress_exponent':2., 'progress_ratio':5., 'length_um':750.},
        {'id':'M05', 'minimum_gap_um':.80, 'aux_start_um':.35, 'aux_stop_um':.42,
         'main_start_um':1.35, 'main_stop_um':1.05, 'progress_exponent':2., 'progress_ratio':5., 'length_um':750.},
        {'id':'M06', 'minimum_gap_um':.80, 'aux_start_um':.35, 'aux_stop_um':.42,
         'main_start_um':1.45, 'main_stop_um':1.05, 'progress_exponent':2., 'progress_ratio':5., 'length_um':750.},
    ]
    state = {'status':'starting', 'created_at':time.strftime('%Y-%m-%d %H:%M:%S'),
             'dimension':'3D EME', 'wavelength_nm':1550., 'candidates':candidates,
             'completed':[], 'full_chain_pass':False,
             'constraint':'minimum core-edge coupling gap > 0.60 um; auxiliary widths >= 0.35 um',
             'screening_note':'粗筛使用1/11/1 EME单元；仅用于选择候选，不替代最终精细收敛',
             'preset':args.preset}
    def save(**kw):
        state.update(kw);state['updated_at']=time.strftime('%Y-%m-%d %H:%M:%S')
        (args.output_dir/'扫描状态.json').write_text(json.dumps(state,ensure_ascii=False,indent=2),encoding='utf-8')
    save()
    for index, candidate in enumerate(candidates, 1):
        tag = f"{candidate['id']}_gap{candidate['minimum_gap_um']:g}_aux{candidate['aux_start_um']:g}to{candidate['aux_stop_um']:g}_main{candidate.get('main_start_um',1.4):g}to{candidate.get('main_stop_um',1.1):g}_L{candidate['length_um']:g}".replace('.', 'p')
        directory=args.output_dir/tag; directory.mkdir()
        geometry=directory/'几何.csv'; smatrix=directory/'S矩阵.csv'; project=directory/'EME.lms'
        command=[sys.executable,'-X','utf8',str(ROOT/'models/optical/mode_mux_eme.py'),
                 '--wavelength-nm','1550','--length-um',str(candidate['length_um']),
                 '--approach-length-um','100','--coupling-length-um',str(candidate['length_um']-200),
                 '--separation-length-um','100','--main-start-um',str(candidate.get('main_start_um',1.4)),'--main-stop-um',str(candidate.get('main_stop_um',1.1)),
                 '--aux-start-um',str(candidate['aux_start_um']),'--aux-stop-um',str(candidate['aux_stop_um']),
                 '--end-gap-um','2.5','--minimum-gap-um',str(candidate['minimum_gap_um']),
                 '--progress-exponent',str(candidate['progress_exponent']),'--progress-ratio',str(candidate['progress_ratio']),
                 '--trajectory',candidate.get('trajectory','asymmetric'),'--geometry-samples','121','--vertical-slices','4',
                 '--y-span-um','12','--mesh-cells-y','180','--mesh-cells-z','90','--modes','8',
                 '--cells','1,11,1','--port1-modes','1,3,4','--port2-modes','1,2,4',
                 '--project',str(project),'--geometry-output',str(geometry),'--s-output',str(smatrix)]
        state.update(status='running',current_index=index,current_case=tag);save()
        log=directory/'控制器输出.log'
        with log.open('x',encoding='utf-8') as stream:
            code=subprocess.run(command,cwd=ROOT,stdout=stream,stderr=subprocess.STDOUT).returncode
        if code or not smatrix.is_file():
            state['completed'].append({'id':candidate['id'],'directory':str(directory),'status':'failed','returncode':code});save();continue
        s=read_s(smatrix)
        row={'id':candidate['id'],'directory':str(directory),'status':'completed',**candidate,
             'port1_main_TE0_to_port2_main_TE0':s.get((4,1)),
             'port1_main_TE1_to_port2_aux_TE0':s.get((5,2)),
             'port1_aux_TE0_to_port2_main_TE1':s.get((6,3)),
             'port2_main_TE0_to_port1_main_TE0':s.get((1,4)),
             'port2_aux_TE0_to_port1_main_TE1':s.get((2,5)),
             'port2_main_TE1_to_port1_aux_TE0':s.get((3,6)),
             's_matrix_shape':[6,6],
             'selected_port_modes':{'port_1':[1,3,4], 'port_2':[1,2,4]},
             'interpretation':'端点FDE已确认port_2辅助TE0为solver mode 2；此处仍需最终端口场形复核'}
        state['completed'].append(row);save()
    valid=[x for x in state['completed'] if x.get('status')=='completed']
    valid.sort(key=lambda x:min(x['port1_main_TE1_to_port2_aux_TE0'] or 0,
                                x['port2_aux_TE0_to_port1_main_TE1'] or 0),reverse=True)
    state.update(status='screening_completed',current_case=None,ranking=valid,
                 best_candidate=valid[0] if valid else None,
                 next_step='对排名最高候选执行精细EME并导出端口场形')
    save()
    (args.output_dir/'筛选排名.json').write_text(json.dumps(valid,ensure_ascii=False,indent=2),encoding='utf-8')


if __name__=='__main__':main()
