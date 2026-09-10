"""合并GDS的独立回读、时延账目和中文预览；不启动仿真。"""
import argparse,json,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
KLAYOUT=Path('C:/Users/PC/AppData/Roaming/KLayout/klayout_app.exe')
ap=argparse.ArgumentParser(description=__doc__)
ap.add_argument('--directory',type=Path,required=True);ap.add_argument('--routes',type=Path,required=True);ap.add_argument('--ysj',type=Path,required=True)
a=ap.parse_args();d=a.directory.resolve()
bp=json.loads((d/'构建参数.json').read_text(encoding='utf-8'));mf=json.loads((d/'合并清单.json').read_text(encoding='utf-8'))
report=d/'合并回读_完整检查.json'
if report.exists():raise FileExistsError(report)
payload={'merged_gds':mf['merged_gds'],'ours_top':bp['top_name'],'ysj_top':bp['ysj_top'],
 'ysj_source':str(a.ysj.resolve()),'report':str(report),'preview_json':str(d/'合并回读_四程几何.json'),
 'shift_um':bp['shift_um'],'source_gds':bp['source_gds'],'build_payload':str(d/'构建参数.json'),
 'mux_name':f"MUX_TE01_A70_A{bp['active_taper_um']:g}_E{bp['external_taper_um']:g}_BODY750_DRAFT",
 'routes_file':str(a.routes.resolve()),'spacing_report':str(d/'实际GDS路由间距检查.json')}
with (d/'回读参数.json').open('x',encoding='utf-8') as f:json.dump(payload,f,ensure_ascii=False,indent=2)
r=subprocess.run([str(KLAYOUT),'-b','-r',str(ROOT/'models/layout/audit_merged_block.py'),'-rd','payload_file='+str(d/'回读参数.json')],capture_output=True)
if r.returncode:raise RuntimeError(r.stderr.decode('utf-8',errors='replace'))
data=json.loads(report.read_text(encoding='utf-8'))
bad={k:v for k,v in data['checks'].items() if (v is False or ('_um2' in k and v!=0))}
expected=mf.get('native_drc_counts') or {'core_width_030':40,'core_space_030':0,'etch_width_030':0,'etch_space_030':40,
 'trench_width_020':0,'trench_space_020':0,'metal_width_2':0,'metal_space_3':0}
if bad or data['counts']!=expected:raise RuntimeError(json.dumps({'bad':bad,'counts':data['counts']},ensure_ascii=False))
if mf.get('layout_variant')=='horizontal_access':
    s=subprocess.run([str(KLAYOUT),'-b','-r',str(ROOT/'models/layout/audit_route_spacing_gds.py'),'-rd','payload_file='+str(d/'回读参数.json')],capture_output=True)
    if s.returncode:raise RuntimeError(s.stdout.decode('utf-8',errors='replace')[-5000:]+s.stderr.decode('utf-8',errors='replace')[-1500:])
subprocess.run([sys.executable,'-X','utf8',str(ROOT/'scripts/merged_layout_delay_budget.py'),'--routes',str(a.routes),'--output',str(d/'合并版时延预算.json')],check=True)
subprocess.run([sys.executable,'-X','utf8',str(ROOT/'models/layout/preview_merged_layout.py'),'--directory',str(d)],check=True)
print(json.dumps({'checks':data['checks'],'counts':data['counts'],'bbox':data['combined_bbox_um'],'file':mf['merged_gds']},ensure_ascii=False,indent=2))
