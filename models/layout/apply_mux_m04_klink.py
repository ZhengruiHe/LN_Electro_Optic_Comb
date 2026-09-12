import argparse
import hashlib
import json
from pathlib import Path
from klink import KLinkClient

ROOT = Path(__file__).resolve().parents[2]
parser = argparse.ArgumentParser()
parser.add_argument('--source', type=Path, required=True)
parser.add_argument('--output-dir', type=Path, required=True)
parser.add_argument('--port', type=int, default=8765)
args = parser.parse_args()
args.output_dir.mkdir(parents=True, exist_ok=False)
output = (args.output_dir / 'YSJ与四程15mm_M04_MUX_750um_功能候选.gds').resolve()
report = (args.output_dir / 'M04版图几何检查.json').resolve()
payload = {'source_gds': str(args.source.resolve()), 'output_gds': str(output),
           'report': str(report), 'mux_cell': 'MUX_TE01_A70_A100_E50_BODY750_DRAFT'}
before = hashlib.sha256(args.source.read_bytes()).hexdigest()
code = (ROOT / 'models/layout/m04_klink_payload.py').read_text(encoding='utf-8')
client = KLinkClient(host='127.0.0.1', port=args.port)
client.connect()
try:
    response = client.call('exec.python', {'code': 'payload=' + repr(payload) + '\n' + code})
    if response.get('exception'):
        raise RuntimeError(str(response['exception']))
    file_info = client.call('layout.file_info', {'path': str(output), 'detail': 'counts'})
    client.call('layout.show_file', {'path': str(output), 'mode': 'new'})
finally:
    client.close()
if hashlib.sha256(args.source.read_bytes()).hexdigest() != before:
    raise RuntimeError('源GDS发生变化')
result = json.loads(report.read_text(encoding='utf-8'))
result.update(source_unchanged=True, source_sha256=before, file_info=file_info)
report.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps(result, ensure_ascii=False, indent=2))
