"""白名单检查并可选暂存复现源码；不提交、不推送，不读取凭据或结果文件。"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess

ROOT = Path(__file__).resolve().parents[1]
EXTENSIONS = {'.py', '.md', '.json', '.lydrc', '.ps1'}
ROOT_FILES = {'.gitignore', '.gitattributes', 'README.md', 'requirements-hfss.txt',
              'requirements-layout.txt', 'requirements-analysis.txt'}
TOKEN_PATTERNS = [r'sk-proj-[A-Za-z0-9_-]{20,}', r'github_pat_[A-Za-z0-9_]{20,}',
                  r'ghp_[A-Za-z0-9]{20,}', r'-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----']


def allowed(name):
    path = Path(name)
    if name in ROOT_FILES:
        return True
    return bool(path.parts and path.parts[0] in ('models', 'scripts', 'docs') and path.suffix in EXTENSIONS
                and not re.search(r'(?:recovery|credential|\.env|\.pem|\.pfx)', path.name, re.I))


def git(*arguments):
    return subprocess.check_output(['git', *arguments], cwd=ROOT)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--report', type=Path, required=True)
    parser.add_argument('--stage', action='store_true')
    args = parser.parse_args()
    if args.report.exists():
        raise FileExistsError('不覆盖既有检查清单')
    already_staged = git('diff', '--cached', '--name-only', '-z').decode('utf-8').split('\0')
    if any(name for name in already_staged if not allowed(name)):
        raise ValueError('暂存区包含白名单以外文件，不自动修改用户暂存状态')
    tracked = git('ls-files', '-z').decode('utf-8').split('\0')
    new = git('ls-files', '--others', '--exclude-standard', '-z').decode('utf-8').split('\0')
    candidates = sorted({name for name in tracked + new if name and allowed(name)})
    records, issues = [], []
    for name in candidates:
        file = ROOT / name
        if not file.is_file():
            raise ValueError('候选文件被删除，不在本次自动暂存删除：' + name)
        if not file.resolve().is_relative_to(ROOT):
            raise ValueError('文件实际路径越出工作区')
        if file.stat().st_size > 2 * 1024**2:
            raise ValueError('源码候选异常偏大：' + name)
        data = file.read_bytes()
        text = data.decode('utf-8-sig')
        for expression in TOKEN_PATTERNS:
            if re.search(expression, text):
                issues.append({'file': name, 'issue': '疑似凭据，只报告文件名，不输出匹配内容'})
        if re.search(r'xwechat_files[/\\]|wxid_[A-Za-z0-9]{8,}', text):
            issues.append({'file': name, 'issue': '包含与公开复现无关的私人路径'})
        records.append({'path': name, 'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()})
    report = {'scope': '仅源码/配置/中文说明白名单；不含结果、GDS、PDK、黑盒和私人文件',
              'checked_files': records, 'issues': issues, 'staged': False, 'pushed': False}
    if not issues and args.stage:
        for offset in range(0, len(candidates), 35):
            subprocess.run(['git', 'add', '--', *candidates[offset:offset+35]], cwd=ROOT, check=True)
        report['staged'] = True
        report['staged_files'] = [v for v in git('diff', '--cached', '--name-only', '-z').decode('utf-8').split('\0') if v]
        if not all(allowed(name) for name in report['staged_files']):
            raise ValueError('暂存区最终白名单核查失败')
    args.report.parent.mkdir(parents=True, exist_ok=True)
    with args.report.open('x', encoding='utf-8') as stream:
        json.dump(report, stream, ensure_ascii=False, indent=2)
    print(json.dumps({'checked_count': len(records), 'issues': issues,
                      'staged_count': len(report.get('staged_files', [])), 'pushed': False}, ensure_ascii=False, indent=2))
    if issues:
        raise SystemExit(2)


if __name__ == '__main__':
    main()
