"""通过KLink生成黑盒交叉器工作稿，保留所有历史文件和原生黑盒层级。

只替换交叉器及其四个外接锥形，并补齐两只端面耦合器的宽度过渡。
不宣称截面、光学性能、时延或整片DRC已通过。
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from klink import KLinkClient
from shapely.geometry import Polygon, box

from check_crossing_candidate import build
from route_geometry import check_route_network

ROOT = Path(__file__).resolve().parents[2]
BASELINE = ROOT / "results/layout/_历史归档/迭代版本_20260912/自定义交叉_初稿_v4/四程10GHz_15mm_自定义交叉_待验证工作稿.gds"
BLACKBOX = ROOT / "SiN_TFLN_V0.1_BlackBox/crossing/crossing_ln.gds"
TOP = "EO4P_10G_15MM_BB_TAPER_DRAFT"


def taper_vertices(length):
    if not 0 < length <= 120:
        raise ValueError("当前固定交叉位置只允许0<L<=120um，避免锥形侵入东侧弯曲")
    x = np.linspace(0,length,201)
    u = x/length
    width = .7 + .5*(3*u*u-2*u*u*u)
    return np.vstack((np.column_stack((x,-width/2)),
                      np.column_stack((x[::-1],width[::-1]/2)))).tolist()


KLAYOUT_CODE = r'''
import json
import hashlib

def signature(layout, cell):
    data = []
    for li in layout.layer_indices():
        shapes = sorted(str(s) for s in cell.shapes(li).each())
        if shapes:
            data.append((str(layout.get_info(li)), shapes))
    instances = sorted((i.cell.name,str(i.dcplx_trans)) for i in cell.each_inst())
    return hashlib.sha256(json.dumps([sorted(data),instances],sort_keys=True).encode()).hexdigest()

def region_box(coords):
    return pya.Region(pya.DBox(*coords).to_itype(layout.dbu))

layout = pya.Layout()
layout.read(payload['baseline'])
top = layout.cell('EO4P_10G_15MM_CUSTOM_CROSSING')
if top is None or layout.dbu != .001 or len(layout.top_cells()) != 1:
    raise RuntimeError('基线顶层或绘制格点不匹配')
baseline_cells = {c.name:signature(layout,c) for c in layout.each_cell() if c != top}
baseline_instances = sorted((i.cell.name,str(i.dcplx_trans)) for i in top.each_inst())
top.name = payload['top']
source = pya.Layout()
source.read(payload['blackbox'])
source_cell = source.cell('crossing_ln')
if source_cell is None or source.dbu != layout.dbu:
    raise RuntimeError('原生黑盒单元不匹配')
if str(source_cell.dbbox()) != '(0,-22.5;45,22.5)':
    raise RuntimeError('原生黑盒尺寸与已核对45x45um不一致')
source_pins = sorted(str(s) for s in source_cell.shapes(source.layer(70,30)).each())
expected_pins = sorted(['box (0,-600;300,600)', 'box (21900,-22500;23100,-22200)',
                        'box (44700,-600;45000,600)', 'box (21900,22200;23100,22500)'])
if source_pins != expected_pins:
    raise RuntimeError('原生GDS四个1.2um端口标记已变化')
if layout.cell('crossing_ln') is not None:
    raise RuntimeError('基线已有同名黑盒，拒绝覆盖')
bb = layout.create_cell('crossing_ln')
bb.copy_tree(source_cell)
if signature(source,source_cell) != signature(layout,bb):
    raise RuntimeError('复制后的黑盒内容变化')

length = payload['taper_length_um']
taper = layout.create_cell('LN_TAPER_070_120_L'+str(int(length))+'_DRAFT')
taper.shapes(layout.layer(20,0)).insert(pya.DPolygon([pya.DPoint(*p) for p in payload['taper']]).to_itype(layout.dbu))
for layer,width in (((20,1),7.2),((10,2),17.2)):
    taper.shapes(layout.layer(*layer)).insert(pya.DBox(0,-width/2,length,width/2).to_itype(layout.dbu))
wrapper = layout.create_cell('LN_CROSSING_BB_WITH_TAPERS_DRAFT')
wrapper.insert(pya.DCellInstArray(bb.cell_index(),pya.DCplxTrans(1,0,False,-22.5,0)))
ports = [('west',0,-22.5-length,0),('east',180,22.5+length,0),
         ('north',270,0,22.5+length),('south',90,0,-22.5-length)]
for name,angle,x,y in ports:
    wrapper.insert(pya.DCellInstArray(taper.cell_index(),pya.DCplxTrans(1,angle,False,x,y)))

# 在独立内存副本中只裁掉六处过渡窗口与旧交叉器，绝不清空活动页。
cx,cy = 300.,250.
cross_window = region_box((cx-22.5,cy-22.5,cx+22.5,cy+22.5))
window = cross_window + region_box((cx-22.5-length,cy-8.6,cx+22.5+length,cy+8.6))
window += region_box((cx-8.6,cy-22.5-length,cx+8.6,cy+22.5+length))
window += region_box((0,250-8.6,length,250+8.6))
window += region_box((18700-length,-26.6-8.6,18700,-26.6+8.6))
original = {}
for layer in ((20,0),(20,1),(10,2)):
    li = layout.layer(*layer)
    original[layer] = pya.Region(top.begin_shapes_rec(li)).merged()
    direct = pya.Region(top.shapes(li)).merged()
    new_region = direct-window
    top.shapes(li).clear()
    top.shapes(li).insert(new_region)
top.insert(pya.DCellInstArray(wrapper.cell_index(),pya.DCplxTrans(1,0,False,cx,cy)))
# 耦合器的片内端口也是1.2um；必须同时修正，否则仍是0.7对1.2硬接。
top.insert(pya.DCellInstArray(taper.cell_index(),pya.DCplxTrans(1,180,False,length,250)))
top.insert(pya.DCellInstArray(taper.cell_index(),pya.DCplxTrans(1,0,False,18700-length,-26.6)))
texts = top.shapes(layout.layer(101,0))
texts.insert(pya.DText('BLACKBOX_CROSSING_TAPERS_DRAFT_NOT_TAPEOUT',1800,1010).to_itype(layout.dbu))
texts.insert(pya.DText('LN2_INTERFACE_AND_OPTICAL_DELAY_PENDING',1800,990).to_itype(layout.dbu))

checks = {}
final_regions = {}
for layer in ((20,0),(20,1),(10,2)):
    final_regions[layer] = pya.Region(top.begin_shapes_rec(layout.layer(*layer))).merged()
    checks['unchanged_outside_edit_windows_'+str(layer)] = ((final_regions[layer]^original[layer])-window).is_empty()
    checks['no_custom_mask_in_crossing_blackbox_'+str(layer)] = (final_regions[layer]&cross_window).is_empty()
checks['inherited_cells_unchanged'] = all(signature(layout,layout.cell(name)) == digest for name,digest in baseline_cells.items())
checks['inherited_instances_unchanged'] = all(item in sorted((i.cell.name,str(i.dcplx_trans)) for i in top.each_inst()) for item in baseline_instances)
checks['blackbox_cell_unchanged'] = signature(source,source_cell) == signature(layout,bb)

core = final_regions[(20,0)]
connections = []
for name,angle,x,y in ports:
    tr = pya.DCplxTrans(1,angle,False,cx+x,cy+y)
    connections.append((name,tr))
connections.extend([('edge_input',pya.DCplxTrans(1,180,False,length,250)),
                    ('edge_output',pya.DCplxTrans(1,0,False,18700-length,-26.6))])
port_report = []
for name,tr in connections:
    # 两端向锥形内取1nm以及窄端外取1nm，验证与原路由连续且端面宽度匹配。
    flags = {}
    for key,coords in (('wide_end_1p2', (length-.001,-.6,length,.6)),
                       ('narrow_end_0p7',(0,-.35,.001,.35)),
                       ('narrow_route_connected',(-.001,-.35,0,.35))):
        test = pya.Region(pya.DPolygon(pya.DBox(*coords)).transformed(tr).to_itype(layout.dbu))
        flags[key] = (test-core).is_empty()
    pos = tr*pya.DPoint(length,0)
    port_report.append({'port':name,'wide_end_um':[pos.x,pos.y],'checks':flags})
checks['six_taper_endpoints_connected'] = all(all(r['checks'].values()) for r in port_report)

# 检查实际黑盒边界；排除仅作RF预留框的100/0图形。
bb_regions = []
for inst in top.each_inst():
    if 'edge_coupler_9um_y_1550_ln' in inst.cell.name:
        b = pya.Region(inst.cell.shapes(layout.layer(100,0)))
        bb_regions.append(b.transformed(inst.cplx_trans))
for layer,reg in final_regions.items():
    checks['no_custom_mask_in_edge_blackboxes_'+str(layer)] = all((reg&br).is_empty() for br in bb_regions)
if len(bb_regions) != 2:
    raise RuntimeError('未找到两个原有端面耦合器黑盒')
if not all(checks.values()):
    raise RuntimeError('连接检查未通过：'+json.dumps(checks))

roi = window.sized(1000)
local_core = core&roi
local_etch = (final_regions[(20,1)]-core)&roi
geometry_markers = {'local_core_width_030':local_core.width_check(300).size(),
                    'local_core_space_030':local_core.space_check(300).size(),
                    'local_etch_width_030':local_etch.width_check(300).size(),
                    'local_etch_space_030':local_etch.space_check(300).size()}
if any(geometry_markers.values()):
    raise RuntimeError('新增接续区域手册几何规则未通过：'+json.dumps(geometry_markers))
full_markers = {'core_width_030':core.width_check(300).size(),
               'core_space_030':core.space_check(300).size(),
               'etch1_width_030':(final_regions[(20,1)]-core).width_check(300).size(),
               'etch1_space_030':(final_regions[(20,1)]-core).space_check(300).size()}
layout.write(payload['gds'])

# 独立预览副本仅用于绘图，不写为可提交工艺层文件，不改变GDS层级。
preview = []
for layer in ((10,2),(20,1),(20,0),(42,0),(100,0)):
    reg = pya.Region(top.begin_shapes_rec(layout.layer(*layer))).merged()
    shapes = []
    for polygon in reg.each():
        shapes.append({'exterior':[[p.x*layout.dbu,p.y*layout.dbu] for p in polygon.each_point_hull()],
                       'holes':[[[p.x*layout.dbu,p.y*layout.dbu] for p in polygon.each_point_hole(k)] for k in range(polygon.holes())]})
    preview.append({'layer':list(layer),'polygons':shapes})
with open(payload['preview_json'],'w',encoding='utf-8') as stream:
    json.dump(preview,stream)
result = {'status':'黑盒交叉器拉锥工作稿，非流片最终版','top_cell':top.name,'top_cells':[c.name for c in layout.top_cells()],
          'bbox_um':[top.dbbox().left,top.dbbox().bottom,top.dbbox().right,top.dbbox().top],
          'blackbox_cell':'crossing_ln','blackbox_center_um':[cx,cy],'blackbox_rotation_deg':0,
          'taper_length_um':length,'number_of_new_tapers':len(connections),'ports':port_report,
          'checks':checks,'local_rule_markers':geometry_markers,'inherited_full_core_etch_rule_markers':full_markers,
          'blackbox_signature':signature(source,source_cell),'hierarchy_preserved':True}
result
'''


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir",type=Path,required=True)
    parser.add_argument("--taper-length-um",type=float,default=100.)
    args = parser.parse_args()
    vertices = taper_vertices(args.taper_length_um)
    taper = Polygon(vertices)
    if not taper.is_valid or not box(0,-.6,args.taper_length_um,.6).covers(taper):
        raise ValueError("拉锥轮廓无效")
    budget_file = ROOT/"results/layout/四程电光梳_10GHz_15mm_3T_3.5T_3T_名义_检查.json"
    budget = json.loads(budget_file.read_text(encoding="utf-8"))
    routes = build([r["actual_centerline_length_um"] for r in budget["loop_centerline_checks"]],True)
    topology = check_route_network(routes,[{"routes":["input","loop2"],"center_um":[300.,250.],"angle_deg":90}])
    if not topology["topology_screen_pass"]:
        raise ValueError("原有中心线拓扑检查失败")
    args.output_dir.mkdir(parents=True,exist_ok=False)
    output = (args.output_dir/"四程10GHz_15mm_PDK黑盒交叉_拉锥接续_待验证工作稿.gds").resolve()
    sources = {str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in (BASELINE,BLACKBOX,budget_file)}
    payload = {"baseline":str(BASELINE),"blackbox":str(BLACKBOX),"gds":str(output),"top":TOP,
               "taper_length_um":args.taper_length_um,"taper":vertices,
               "preview_json":str((args.output_dir/"预览几何.json").resolve())}
    client = KLinkClient()
    client.connect()
    try:
        response = client.call("exec.python",{"code":"payload = "+repr(payload)+"\n"+KLAYOUT_CODE})
        if response.get("exception"):
            raise RuntimeError(str(response["exception"]))
        result = response["return_value"]
        result["file_info"] = client.call("layout.file_info",{"path":str(output),"detail":"counts"})
    finally:
        client.close()
    if not all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==h for p,h in sources.items()):
        raise RuntimeError("输入文件变化，请核对；未覆盖原文件")
    result.update({"working_gds":str(output),"source_sha256":sources,"source_files_unchanged":True,
                   "topology":topology,"centerline_route_lengths_unchanged":True,
                   "optical_validation_pass":False,"foundry_signoff":False,
                   "limitations":["黑盒内部无本地可调用S矩阵；手册<0.2dB仅参考，未按晶向分列",
                                  "100um拉锥尚未完成EME/FDTD；0.7um是自定义候选而非PDK默认",
                                  "7.2um为LN1刻蚀包络，不能自动解释为LN2平台；标准端口截面待确认",
                                  "保留原MUX脊拉锥及窄尾端，尾端0.15um相关历史违规未修正",
                                  "原MUX外围窗口接缝未改变，不能按外观判断实际材料连续性",
                                  "没有推测添加LN2全片掩膜；不代表原仿真截面已在工艺上实现",
                                  "保持几何中心线不等于保持群时延；10GHz预算需回填真实过渡与黑盒时延",
                                  "原有RF焊盘/终端、完整级联及官方DRC仍未完成"]})
    (args.output_dir/"版图检查.json").write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({k:v for k,v in result.items() if k in ("working_gds","checks","ports","local_rule_markers","inherited_full_core_etch_rule_markers","bbox_um")},ensure_ascii=False,indent=2))


if __name__ == "__main__":
    main()
