"""整理已成功完成但被文件名后处理误报的辅助脊0.35um MUX EME结果；不重跑。"""
from __future__ import annotations
import csv
import hashlib
import json
from pathlib import Path
import re
import time


def read_s(path):
    rows = list(csv.DictReader(path.open(encoding='utf-8-sig', newline='')))
    return {(int(r['输出索引']), int(r['输入索引'])): float(r['功率']) for r in rows}


def main():
    import argparse
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--directory', type=Path, required=True)
    a = p.parse_args(); d = a.directory.resolve()
    log = d / 'MUX_aux0p35_gap1um_EMEplms_p0.log'
    geom = d / 'MUX_aux0p35_gap1um_几何pcsv'
    smat = d / 'MUX_aux0p35_gap1um_S矩阵pcsv'
    project = d / 'MUX_aux0p35_gap1um_EMEplms.lms'
    if not all(x.is_file() for x in (log, geom, smat, project)):
        raise FileNotFoundError('缺少已生成的EME文件')
    log_text = log.read_text(encoding='utf-8', errors='replace')
    if 'Simulation completed successfully' not in log_text:
        raise RuntimeError('EME日志没有成功结束标志')
    with geom.open(encoding='utf-8-sig', newline='') as f:
        rows = list(csv.DictReader(f))
    gaps = [float(x['gap_um']) for x in rows]
    aux = [float(x['auxiliary_width_um']) for x in rows]
    main = [float(x['main_width_um']) for x in rows]
    if abs(min(gaps)-1.0)>1e-9 or abs(min(aux)-.35)>1e-9 or abs(max(aux)-.4)>1e-9:
        raise RuntimeError('生成几何与目标参数不一致')
    s = read_s(smat)
    # 端口选择顺序为[solver mode 1,3,4]；6x6用户S的前3行/列为port_1，
    # 后3行/列为port_2。物理通道按已核对的基线映射取值。
    channels = {
        'port1_main_TE0_to_port2_main_TE0': s[(4,1)],
        'port1_main_TE1_to_port2_aux_TE0': s[(5,2)],
        'port1_aux_TE0_to_port2_main_TE1': s[(6,3)],
        'port2_main_TE0_to_port1_main_TE0': s[(1,4)],
        'port2_aux_TE0_to_port1_main_TE1': s[(2,5)],
        'port2_main_TE1_to_port1_aux_TE0': s[(3,6)],
    }
    baseline_path = Path(__file__).resolve().parents[2] / 'results/optical/模式复用器_70度_EME_精细检查_S矩阵.csv'
    baseline = read_s(baseline_path)
    baseline_channels = {
        'port1_main_TE0_to_port2_main_TE0': baseline[(4,1)],
        'port1_main_TE1_to_port2_aux_TE0': baseline[(5,2)],
        'port1_aux_TE0_to_port2_main_TE1': baseline[(6,3)],
        'port2_main_TE0_to_port1_main_TE0': baseline[(1,4)],
        'port2_aux_TE0_to_port1_main_TE1': baseline[(2,5)],
        'port2_main_TE1_to_port1_aux_TE0': baseline[(3,6)],
    }
    result = {
        'status': 'eme_completed_postprocessed', 'wavelength_nm': 1550.,
        'auxiliary_start_width_um': min(aux), 'auxiliary_end_width_um': max(aux),
        'minimum_core_edge_gap_um': min(gaps), 'maximum_core_edge_gap_um': max(gaps),
        'geometry_rows': len(rows), 'eme_log_success': True,
        'channels_power': channels, 'baseline_channels_power': baseline_channels,
        'channel_change_power': {k: channels[k]-baseline_channels[k] for k in channels},
        'channel_metrics_percent': {k: 100*v for k,v in channels.items()},
        'selected_port_modes': {'port_1': [1,3,4], 'port_2': [1,3,4]},
        'interpretation': '按已核对的基线端口场形映射；完整6x6 S已保存。辅助脊加宽后TE1转换通道明显下降，不作为当前MUX替代版。',
        'solver_project': str(project), 'solver_log': str(log),
        'baseline_s_sha256': hashlib.sha256(baseline_path.read_bytes()).hexdigest(),
        'full_chain_pass': False, 'tolerance_sweep': False,
        'limitations': ['EME模式数/单元数未做独立收敛扫描', '只改辅助脊起始宽度，未重优化渐变轨迹/耦合长度',
                        'PDK标记需代工确认；该结果不能直接证明流片合规'],
    }
    (d/'结果.json').write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding='utf-8')
    state_path = d/'状态.json'; state = json.loads(state_path.read_text(encoding='utf-8'))
    state.update(status='completed_postprocessed', child_pid=None, returncode=0,
                 postprocess_warning='原控制器因扩展名替换误报缺文件；EME日志、几何和S矩阵已核实',
                 result_path=str((d/'结果.json').resolve()), finished_at=time.strftime('%Y-%m-%d %H:%M:%S'))
    state_path.write_text(json.dumps(state, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == '__main__': main()
