"""发布前只读核对：回读实际GDS、源保留、间距、规则标记；不生成新版图或启动仿真。"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
from build_15mm_ysj_merge import structures

ROOT = Path(__file__).resolve().parents[2]


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path, data):
    with path.open('x', encoding='utf-8') as stream:
        json.dump(data, stream, ensure_ascii=False, indent=2)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--klayout', type=Path, default=Path('C:/Users/PC/AppData/Roaming/KLayout/klayout_app.exe'))
    args = parser.parse_args()
    source_dir = args.directory.resolve()
    payload = read(source_dir / '回读参数.json')
    manifest = read(source_dir / '合并清单.json')
    original = read(source_dir / '合并回读_完整检查.json')
    sources = {key: Path(payload[key]) for key in ('merged_gds', 'source_gds', 'ysj_source', 'build_payload', 'routes_file')}
    before = {key: digest(path) for key, path in sources.items()}
    if before['merged_gds'] != manifest['sha256']:
        raise ValueError('当前GDS与所选版本清单不同，拒绝用旧报告代表当前版图')
    args.output_dir.mkdir(parents=True, exist_ok=False)
    payload['report'] = str((args.output_dir / '实际GDS最终回读.json').resolve())
    payload['spacing_report'] = str((args.output_dir / '实际GDS最终间距.json').resolve())
    payload.pop('preview_json', None)
    request = args.output_dir / '只读检查参数.json'
    write(request, payload)
    for script in ('audit_merged_block.py', 'audit_route_spacing_gds.py'):
        result = subprocess.run([str(args.klayout), '-b', '-r', str(ROOT / 'models/layout' / script),
                                 '-rd', 'payload_file=' + str(request.resolve())], capture_output=True)
        if result.returncode:
            raise RuntimeError(script + '执行失败：' + result.stderr.decode('utf-8', errors='replace'))
    check = read(Path(payload['report']))
    spacing = read(Path(payload['spacing_report']))
    old_cells, old_units, _ = structures(sources['ysj_source'].read_bytes())
    merged_cells, merged_units, _ = structures(sources['merged_gds'].read_bytes())
    records = {name: merged_cells.get(name) == value for name, value in old_cells.items()}
    after = {key: digest(path) for key, path in sources.items()}
    failures = {k: v for k, v in check['checks'].items() if v is False or ('_um2' in k and v != 0)}
    result = {
        'scope': '用户接受现有仿真的阶段冻结；只读版图核查，不是流片签核',
        'merged_gds_sha256': before['merged_gds'],
        'source_files_unchanged': before == after,
        'YSJ_records_byte_identical': records, 'GDS_units_identical': old_units == merged_units,
        'geometry_failures': failures, 'checks': check['checks'],
        'ordinary_route_spacing_pass': spacing['ordinary_route_spacing_pass'],
        'drc_counts': check['counts'], 'drc_counts_match_previous': check['counts'] == original['counts'],
        'all_checked_drc_zero': all(v == 0 for v in check['counts'].values()),
        'GSG_net_continuity': check['GSG_net_continuity'],
        'combined_bbox_um': check['combined_bbox_um'], 'ours_bbox_um': check['ours_bbox_um'],
        'geometry_stage_pass': not failures and before == after and all(records.values()) and spacing['ordinary_route_spacing_pass'],
        'full_chain_pass': False, 'tapeout_signoff': False,
        'clarifications': [
            '归一化窗口模式下MUX_body_layers_match_frozen_source只比较20/0核心；不表示20/1、10/2从未改变',
            '40条核心线宽、40条刻蚀间距旧标记保留，未被豁免或删除',
            '用户接受现有二维仿真，不等于把status=1改成能量停止或数值收敛',
            '本次不安排MUX外接续、黑盒、射频接口或容差补算；相关性能证据边界如实保留',
        ],
    }
    write(args.output_dir / '最终核查摘要.json', result)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    if not result['geometry_stage_pass']:
        raise SystemExit(2)


if __name__ == '__main__':
    main()
