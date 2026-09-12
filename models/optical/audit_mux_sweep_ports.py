"""逐例只读核查MUX两端的矢量模场，按物理分支重建通道指标。

不调用run/findmodes/updateportmodes/emepropagate，不保存或修改LMS。
结果只写新目录，旧批次的猜测模式标签保持原样供追溯。
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import numpy as np
from mode_pdk_sweep import DEFAULT_CONFIG, load_config, load_lumapi
from mode_mux_eme import extract_s_matrix


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_new(path, data):
    with path.open('x', encoding='utf-8') as stream:
        json.dump(data, stream, ensure_ascii=False, indent=2)


def sign_changes(field, y, z, center, width):
    mask = abs(y-center) <= width/2 + .25
    height = (z >= .7) & (z <= 1.1)
    if not mask.any() or not height.any():
        return None
    ey = field[:, :, 1]
    zi = np.flatnonzero(height)[np.argmax(np.sum(abs(ey[mask][:, height])**2, axis=0))]
    line = ey[mask, zi]
    if abs(line).max() == 0:
        return None
    phased = line * np.exp(-1j*np.angle(line[np.argmax(abs(line))]))
    values = phased.real[abs(phased.real) > .12 * abs(phased.real).max()]
    return int(np.count_nonzero(np.diff(np.sign(values))))


def port_metrics(data, endpoint, offset):
    wavelengths = np.asarray(data['lambda']).ravel()*1e9
    if wavelengths.size != 1 or not np.allclose(wavelengths, 1550, atol=1e-6, rtol=0):
        raise ValueError('端口模实际频点不是1550nm')
    y, z = [np.asarray(data[k]).ravel()*1e6 for k in ('y', 'z')]
    weights = abs(np.gradient(y))[:, None] * abs(np.gradient(z))[None, :]
    centers = [float(endpoint[k]) for k in ('main_center_um', 'auxiliary_center_um')]
    widths = [float(endpoint[k]) for k in ('main_width_um', 'auxiliary_width_um')]
    gap = float(endpoint['gap_um'])
    padding = min(.45, gap/2-.001)
    modes = sorted(int(k[1:]) for k in data if re.fullmatch(r'E\d+', k))
    rows = []
    for slot, mode in enumerate(modes):
        e = np.asarray(data[f'E{mode}']).squeeze()
        if e.shape != (len(y), len(z), 3):
            raise ValueError('端口场数组维度与坐标不一致')
        intensity = abs(e)**2
        components = np.sum(intensity*weights[:, :, None], axis=(0, 1))
        components /= sum(components)
        lateral = np.sum(intensity*weights[:, :, None], axis=(1, 2))
        branch = [float(sum(lateral[abs(y-center) <= width/2+padding])/sum(lateral))
                  for center, width in zip(centers, widths)]
        nodes = [sign_changes(e, y, z, center, width) for center, width in zip(centers, widths)]
        labels = []
        if components[1] > .5:
            if branch[0] > .6 and nodes[0] == 0:
                labels.append('main_TE0')
            if branch[0] > .6 and nodes[0] == 1:
                labels.append('main_TE1')
            if branch[1] > .6 and nodes[1] == 0:
                labels.append('aux_TE0')
        rows.append({'solver_mode': mode, 'S_index_1based': offset+slot+1,
                     'component_fraction_xyz': components.tolist(),
                     'main_energy_fraction': branch[0], 'aux_energy_fraction': branch[1],
                     'main_Ey_nodes': nodes[0], 'aux_Ey_nodes': nodes[1],
                     'centroid_y_um': float(sum(y*lateral)/sum(lateral)),
                     'peak_y_um': float(y[np.argmax(lateral)]),
                     'label': labels[0] if len(labels) == 1 else 'ambiguous_or_other'})
    mapping = {}
    for label in ('main_TE0', 'main_TE1', 'aux_TE0'):
        found = [row for row in rows if row['label'] == label]
        if len(found) == 1:
            mapping[label] = found[0]['S_index_1based']
    return {'wavelength_nm': float(wavelengths[0]), 'modes': rows, 'physical_S_indices': mapping,
            'all_three_physical_channels_present': len(mapping) == 3}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--batch', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=False)
    state = json.loads((args.batch/'扫描状态.json').read_text(encoding='utf-8'))
    api = load_lumapi(load_config(DEFAULT_CONFIG))
    reports = []
    for entry in state['completed']:
        directory = Path(entry['directory']).resolve()
        project = directory/entry.get('project_file', 'EME.lms')
        before = digest(project)
        rows = list(csv.DictReader((directory/entry.get('geometry_file', '几何.csv')).open(encoding='utf-8-sig')))
        csv_s = list(csv.DictReader((directory/entry.get('s_file', 'S矩阵.csv')).open(encoding='utf-8-sig')))
        power = np.zeros((6, 6))
        for row in csv_s:
            power[int(row['输出索引'])-1, int(row['输入索引'])-1] = float(row['功率'])
        report = {'id': entry['id'], 'source_project': str(project), 'source_sha256': before,
                  'length_total_um': float(rows[-1]['x_um'])-float(rows[0]['x_um']),
                  'min_gap_um': min(float(row['gap_um']) for row in rows), 'ports': {}}
        print('只读两端：'+entry['id'], flush=True)
        with api.MODE(filename=str(project), hide=True) as mode:
            report['solver_version'] = mode.version()
            report['background_material'] = mode.getnamed('EME', 'background material')
            report['solver_type'] = mode.getnamed('EME', 'solver type')
            report['mesh_cells_y'] = float(mode.getnamed('EME', 'mesh cells y'))
            report['mesh_cells_z'] = float(mode.getnamed('EME', 'mesh cells z'))
            report['modes_per_cell'] = float(mode.getnamed('EME', 'number of modes for all cell groups'))
            report['project_wavelength_nm'] = float(mode.getnamed('EME', 'wavelength'))*1e9
            if not np.isclose(report['project_wavelength_nm'], 1550, atol=1e-6, rtol=0):
                raise ValueError('LMS工程波长不一致')
            report['cells'] = np.asarray(mode.getnamed('EME', 'cells')).ravel().tolist()
            report['group_spans_um'] = (np.asarray(mode.getnamed('EME', 'group spans')).ravel()*1e6).tolist()
            matrix = extract_s_matrix(mode.getresult('EME', 'power normalized user s matrix'))
            report['CSV_matches_saved_LMS_S'] = bool(np.allclose(abs(matrix)**2, power, atol=1e-10, rtol=1e-8))
            for index, (port, endpoint) in enumerate((('port_1', rows[0]), ('port_2', rows[-1]))):
                data = mode.getresult('EME::Ports::'+port, 'mode profiles')
                np.savez_compressed(args.output_dir/(entry['id']+'_'+port+'_已保存端口模.npz'),
                                    **{k: v for k, v in data.items() if isinstance(v, np.ndarray)})
                report['ports'][port] = port_metrics(data, endpoint, index*3)
        report['source_unchanged'] = digest(project) == before
        if not report['source_unchanged'] or not report['CSV_matches_saved_LMS_S']:
            raise RuntimeError('源文件变化或CSV与LMS不一致')
        left = report['ports']['port_1']['physical_S_indices']
        right = report['ports']['port_2']['physical_S_indices']
        channels = {}
        for name, input_label, output_label in (
                ('TE0_straight', 'main_TE0', 'main_TE0'),
                ('main_TE1_to_aux_TE0', 'main_TE1', 'aux_TE0'),
                ('aux_TE0_to_main_TE1', 'aux_TE0', 'main_TE1')):
            if input_label in left and output_label in right:
                i, o = left[input_label]-1, right[output_label]-1
                channels[name] = {'forward_power': float(power[o, i]), 'reverse_power': float(power[i, o]),
                                  'input_S_index': i+1, 'output_S_index': o+1}
            else:
                channels[name] = {'unresolved': True, 'reason': '已保存的选模不包含全部目标模，或物理分支不明确'}
        report['verified_channels'] = channels
        report['all_physical_port_labels_resolved'] = len(left) == 3 and len(right) == 3
        report['full_chain_pass'] = False
        write_new(args.output_dir/(entry['id']+'_两端模式核查.json'), report)
        reports.append(report)
        print(json.dumps({'id': entry['id'], 'port1_labels': left, 'port2_labels': right,
                          'channels': channels}, ensure_ascii=False), flush=True)
    write_new(args.output_dir/'批次两端核查汇总.json', {
        'scope': '逐例矢量场局域性和Ey节点检查；无新求解，不改LMS/GDS',
        'classification_note': '分支积分窗互不重叠；TE主分量、芯区局域性及节点联合识别；混合模不强行贴标签',
        'reports': reports, 'prior_fixed_index_ranking_not_validated': True, 'full_chain_pass': False})


if __name__ == '__main__':
    main()
