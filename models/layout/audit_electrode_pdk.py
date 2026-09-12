"""只读核查GDS中的T形M1电极及手册线宽/间距，不改版图、不启动电磁求解。"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from klink import KLinkClient

ROOT = Path(__file__).resolve().parents[2]
DEFAULT = ROOT / "results/layout/_历史归档/迭代版本_20260912/黑盒交叉_拉锥工作稿_v5/四程10GHz_15mm_PDK黑盒交叉_拉锥接续_待验证工作稿.gds"

CODE = r'''
from collections import Counter
audit_layout = pya.Layout()
audit_layout.read(payload['gds'])
if len(audit_layout.top_cells()) != 1:
    raise RuntimeError('GDS顶层数量不唯一')
audit_top = audit_layout.top_cell()
electrode = audit_layout.cell('T_GSG_15MM_RECTANGULAR_BASELINE')
if electrode is None:
    raise RuntimeError('没有找到当前T形电极单元')
dbu = audit_layout.dbu
metal_layer = audit_layout.layer(42,0)
shapes = list(electrode.shapes(metal_layer).each())
if not all(s.is_box() for s in shapes):
    raise RuntimeError('电极已不是预期矩形组成，请重新审阅而非复用当前识别规则')
boxes = [s.box for s in shapes]
counts = Counter((round(b.width()*dbu,6),round(b.height()*dbu,6)) for b in boxes)
trunks = [b for b in boxes if b.width()*dbu == 15000.]
caps = [b for b in boxes if (round(b.width()*dbu,6),round(b.height()*dbu,6)) == (45.,2.)]
necks = [b for b in boxes if (round(b.width()*dbu,6),round(b.height()*dbu,6)) == (10.,4.)]
centers = sorted(set(round(b.center().x*dbu,6) for b in caps))
cap_tracks = sorted(set((round(b.bottom*dbu,6),round(b.top*dbu,6)) for b in caps))
periods = sorted(set(round(b-a,6) for a,b in zip(centers,centers[1:])))
region = pya.Region(audit_top.begin_shapes_rec(metal_layer)).merged()
pieces = list(region.each())
markers = {'M1_width_below_2um':region.width_check(round(2/dbu)).size(),
           'M1_space_below_3um':region.space_check(round(3/dbu)).size()}
diagnostics = {'width_below_2p001um':region.width_check(round(2.001/dbu)).size(),
               'space_below_4um':region.space_check(round(4/dbu)).size(),
               'space_below_4p001um':region.space_check(round(4.001/dbu)).size()}
checks = {'three_trunks':len(trunks)==3,'1200_caps':len(caps)==1200,'1200_necks':len(necks)==1200,
          '300_period_positions':len(centers)==300,'period_50um':periods==[50.],
          'four_symmetric_cap_tracks':cap_tracks==[(-34.5,-32.5),(-27.5,-25.5),(25.5,27.5),(32.5,34.5)],
          'three_connected_conductors':len(pieces)==3,
          'full_15mm_connected_spans':all(abs(p.bbox().width()*dbu-15000)<1e-6 for p in pieces),
          'all_expected_rectangles':len(boxes)==2403 and sum(counts.values())==2403,
          'no_metal_core_overlap':(region & pya.Region(audit_top.begin_shapes_rec(audit_layout.layer(20,0)))).is_empty(),
          'grid_1nm':dbu==.001,'line_space_rules':not any(markers.values())}
result = {'gds':payload['gds'],'top_cell':audit_top.name,'electrode_cell':electrode.name,'layer':'42/0',
          'box_dimensions_um':[{'length_x':w,'width_y':h,'count':n} for (w,h),n in sorted(counts.items())],
          'period_um':periods,'number_of_periods':len(centers),'cap_tracks_um':cap_tracks,
          'connected_conductor_bboxes_um':[[p.bbox().left*dbu,p.bbox().bottom*dbu,p.bbox().right*dbu,p.bbox().top*dbu] for p in pieces],
          'manual_rule_markers':markers,'threshold_diagnostics':diagnostics,'checks':checks,
          'nominal_geometry_check_pass':all(checks.values())}
result
'''


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--gds',type=Path,default=DEFAULT)
    p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--screenshot',action='store_true',help='仅当活动KLayout文件与被检查文件相同时导出局部截图')
    a = p.parse_args()
    before = hashlib.sha256(a.gds.read_bytes()).hexdigest()
    a.output_dir.mkdir(parents=True,exist_ok=False)
    c = KLinkClient()
    c.connect()
    try:
        r = c.call('exec.python',{'code':'payload = '+repr({'gds':str(a.gds.resolve())})+'\n'+CODE})
        if r.get('exception'):
            raise RuntimeError(str(r['exception']))
        result = r['return_value']
        if a.screenshot:
            active = c.call('layout.info',{'verbosity':'full'})
            if active.get('file') and Path(active['file']).resolve()==a.gds.resolve() and active.get('cell')==result['top_cell']:
                screenshot = (a.output_dir/'T形电极_实际GDS局部.png').resolve()
                c.call('view.screenshot',{'mode':'path','path':str(screenshot),
                                          'bbox_um':[4998,15,5102,43],'width_px':1600,'height_px':520})
                result['screenshot'] = str(screenshot)
            else:
                result['screenshot_skipped'] = '当前活动页不同，不改变用户视图'
    finally:
        c.close()
    after = hashlib.sha256(a.gds.read_bytes()).hexdigest()
    result.update({'gds_sha256':before,'original_gds_unchanged':before==after,
                   'manual_source':'南智光电流片说明-TFLN-on-SiN-V0.2试用版.pdf，4.3节，M1最小线宽2um/间距3um',
                   'scope':'仅当前T形M1名义平面几何，不是全片或电磁性能签核',
                   'foundry_signoff':False,'tolerance_performed':False,
                   'remaining':['帽部2um刚好等于规则下限，未评估制造裕量',
                                '焊盘、RF接入、终端及地连接方案未完成',
                                'GDS不验证金属厚度、电导率、粗糙度及材料堆栈',
                                '未替代完整电磁收敛及官方DRC']})
    (a.output_dir/'T形电极_PDK几何核查.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))
    if not result['nominal_geometry_check_pass'] or before!=after:
        raise SystemExit(2)


if __name__=='__main__':
    main()
