"""局部扩展欧拉端口公共窗口，消除新分叉缝；保留LN1核心、MUX本体及M1。"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess

ROOT=Path(__file__).resolve().parents[2]
KLAYOUT=Path('C:/Users/PC/AppData/Roaming/KLayout/klayout_app.exe')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source-dir',type=Path,required=True)
    p.add_argument('--output-dir',type=Path,required=True)
    a=p.parse_args();a.output_dir=a.output_dir.resolve()
    old=json.loads((a.source_dir/'GSG连接与版图检查.json').read_text(encoding='utf-8'))
    source=Path(old['output_gds'])
    if hashlib.sha256(source.read_bytes()).hexdigest()!=old['output_sha256']:raise ValueError('源GDS已变化')
    if a.output_dir.exists():raise FileExistsError('输出目录已存在')
    a.output_dir.mkdir(parents=True,exist_ok=False)
    pitch=old['pitch_um']
    payload={'source':str(source),'source_sha256':old['output_sha256'],
        'output':str(a.output_dir/f'四程10GHz_15mm_GSG{pitch:g}um_分叉窗口修正_待验证工作稿.gds'),
        'top':f'EO4P_10G_15MM_GSG{pitch:g}_WINDOW_FIX_DRAFT','source_report':old,
        'report':str(a.output_dir/'分叉窗口与GSG检查.json'),
        'preview':str(a.output_dir/'版图预览几何.json'),
        'patches':{'20/1':[[785,25.9,851,44],[785,-44,851,-25.9]],
                   '10/2':[[754,21.1,851,64],[754,-64,851,-21.1],
                           [18149,18,18250,57],[18149,-57,18250,-18]]}}
    path=a.output_dir/'窗口修正参数.json';path.write_text(json.dumps(payload,ensure_ascii=False,indent=2),encoding='utf-8')
    run=subprocess.run([str(KLAYOUT),'-b','-r',str(Path(__file__).with_name('klayout_fix_branch_windows.py')),
                        '-rd','payload_file='+str(path)],capture_output=True,text=True,encoding='utf-8',errors='replace')
    report=json.loads(Path(payload['report']).read_text(encoding='utf-8')) if Path(payload['report']).exists() else {}
    print(json.dumps({k:v for k,v in report.items() if k!='source_report'},ensure_ascii=False,indent=2))
    if run.returncode:raise RuntimeError(run.stderr+run.stdout)
    if hashlib.sha256(source.read_bytes()).hexdigest()!=old['output_sha256']:raise RuntimeError('源文件变化')


if __name__=='__main__':main()
