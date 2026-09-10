"""规划中间版图归档，只生成清单，不移动/删除文件。

逐个检查中间目录是否被目录外的复现脚本、说明或结果元数据引用。
无法确认为独立历史稿的目录留在原位。结果大文件只做哈希，不输出内容。
"""
import argparse
import hashlib
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    base = ROOT / 'results/layout'
    candidates = [p for p in base.iterdir() if p.is_dir() and
                  re.match(r'^YSJ合并15mm_横向接入_v13_(?:07|08|10|12|13|14|15|16|17|18|20)(?:_|$)', p.name)]
    references = {p.name: [] for p in candidates}
    # 限定为项目源文件和结果元数据，不遍历本机个人文件、凭据或PDK。
    for scope in ('models', 'scripts', 'docs', 'results'):
        for file in (ROOT / scope).rglob('*'):
            if not file.is_file() or file.suffix not in ('.py', '.md', '.json'):
                continue
            if file.stat().st_size > 24 * 1024**2 or '_历史归档' in file.parts:
                continue
            text = file.read_text(encoding='utf-8-sig', errors='replace')
            # JSON路径可能用ASCII转义保存；只进行字符串解码以查引用。
            if file.suffix == '.json':
                try:
                    text = json.dumps(json.loads(text), ensure_ascii=False)
                except (ValueError, RecursionError):
                    pass
            for candidate in candidates:
                if candidate.name in text and not file.is_relative_to(candidate):
                    references[candidate.name].append(str(file.relative_to(ROOT)))
    selected, protected = [], []
    archive = base / '_历史归档/20260911_未采用窗口中间稿'
    for candidate in candidates:
        if references[candidate.name]:
            protected.append({'directory': candidate.name, 'references': references[candidate.name]})
            continue
        files = [{'path': str(f.relative_to(candidate)), 'bytes': f.stat().st_size,
                  'sha256': hashlib.sha256(f.read_bytes()).hexdigest()}
                 for f in candidate.rglob('*') if f.is_file()]
        selected.append({'source': str(candidate.resolve()), 'destination': str((archive / candidate.name).resolve()),
                         'files': files, 'bytes': sum(f['bytes'] for f in files)})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x', encoding='utf-8') as stream:
        json.dump({'status': 'planned_not_moved', 'workspace': str(ROOT), 'archive': str(archive.resolve()),
                   'selected': selected, 'protected': protected,
                   'note': '只移动未在其他脚本/说明/元数据中被引用的窗口中间稿；恢复时按原路径移回'},
                  stream, ensure_ascii=False, indent=2)
    print(json.dumps({'selected': [Path(e['source']).name for e in selected],
                      'total_bytes': sum(e['bytes'] for e in selected),
                      'protected': [{'directory': p['directory'], 'reference_count': len(p['references'])} for p in protected]},
                     ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
