"""平滑对齐四个MUX有源端的20/1与10/2包络，另存GDS，不改光学脊或电极。"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from klink import KLinkClient
from shapely.geometry import Polygon, box

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "results/layout/_历史归档/迭代版本_20260912/尾端线宽修正_工作稿_v6/四程10GHz_15mm_尾端线宽修正_待验证工作稿.gds"
TAPER_LENGTH_UM = 200.0
MAIN_CENTER_UM = -1.95


def shifted_envelope_vertices(width_um: float, samples: int = 201) -> list[list[float]]:
    """在x=-200到0间将包络中心从主脊中心平滑移到MUX共同中心。"""
    if width_um <= 0 or samples < 3:
        raise ValueError("包络宽度和采样数必须有效")
    x = np.linspace(-TAPER_LENGTH_UM, 0.0, samples)
    u = (x + TAPER_LENGTH_UM) / TAPER_LENGTH_UM
    progress = u * u * (3.0 - 2.0 * u)
    center = MAIN_CENTER_UM * (1.0 - progress)
    lower = np.column_stack((x, center - 0.5 * width_um))
    lower = np.vstack((lower, [950.0, -0.5 * width_um]))
    upper = np.vstack(([950.0, 0.5 * width_um],
                       np.column_stack((x[::-1], center[::-1] + 0.5 * width_um))))
    return np.vstack((lower, upper)).tolist()


KLAYOUT_CODE = r'''
import hashlib
import json

def signature(layout, cell):
    shapes = [(str(layout.get_info(li)), sorted(str(s) for s in cell.shapes(li).each()))
              for li in layout.layer_indices() if not cell.shapes(li).is_empty()]
    instances = sorted((i.cell.name, str(i.dcplx_trans)) for i in cell.each_inst())
    return hashlib.sha256(json.dumps([sorted(shapes), instances]).encode()).hexdigest()

layout = pya.Layout()
layout.read(payload['input'])
top = layout.top_cell()
if top is None or layout.dbu != .001 or len(layout.top_cells()) != 1:
    raise RuntimeError('输入GDS顶层或1nm格点不匹配')
muxes = [c for c in layout.each_cell() if c.name.startswith('MUX_TE01_A70_TAIL030')]
if len(muxes) != 1:
    raise RuntimeError('未唯一找到v6尾端修正MUX单元')
mux = muxes[0]
instances = [i for i in top.each_inst() if i.cell == mux]
if len(instances) != 4:
    raise RuntimeError('整片MUX实例数不是四个')
unchanged_cells = {c.name: signature(layout, c) for c in layout.each_cell() if c not in (top, mux)}
unchanged_top_layers = {str(layout.get_info(li)): pya.Region(top.begin_shapes_rec(li)).merged()
                        for li in layout.layer_indices() if (layout.get_info(li).layer, layout.get_info(li).datatype)
                        not in ((20,1),(10,2),(101,0))}
before_local = {}
after_local = {}
for layer_key in ((20,1),(10,2)):
    li = layout.layer(*layer_key)
    shapes = list(mux.shapes(li).each())
    matching = [s for s in shapes if s.bbox() == pya.Box(-200000,
                -payload['widths'][str(layer_key)]*500, 950000, payload['widths'][str(layer_key)]*500)]
    if len(matching) != 1 or len(shapes) != 1:
        raise RuntimeError('MUX包络不是预期单矩形：'+str(layer_key))
    before_local[layer_key] = pya.Region(mux.shapes(li)).merged()
    matching[0].delete()
    polygon = pya.DPolygon([pya.DPoint(*xy) for xy in payload['vertices'][str(layer_key)]])
    mux.shapes(li).insert(polygon.to_itype(layout.dbu))
    after_local[layer_key] = pya.Region(mux.shapes(li)).merged()
mux.name = 'MUX_TE01_A70_TAIL030_ENVELOPE_ALIGNED_UNVERIFIED'
top.name = 'EO4P_10G_15MM_MUX_ENVELOPE_ALIGNED_DRAFT'
text_layer = layout.layer(101,0)
mux.shapes(text_layer).insert(pya.Text('MUX_ENVELOPE_SMOOTH_SHIFT_200UM_UNVERIFIED', pya.Trans(-190000,-6800)))

# 10/2半宽8.6um再偏移1.95um，局部审计窗口需完整覆盖到±12um。
edit_window = pya.Region(pya.Box(-200000,-12000,1,12000))
checks = {
    'four_mux_instances': len(instances) == 4,
    'other_cells_unchanged': all(signature(layout,layout.cell(name)) == digest for name,digest in unchanged_cells.items()),
    'other_top_layers_unchanged': all((pya.Region(top.begin_shapes_rec(li)).merged() ^ original).is_empty()
       for li in layout.layer_indices() if str(layout.get_info(li)) in unchanged_top_layers
       for original in [unchanged_top_layers[str(layout.get_info(li))]]),
}
for layer_key in ((20,1),(10,2)):
    checks['only_input_side_changed_'+str(layer_key)] = ((before_local[layer_key]^after_local[layer_key])-edit_window).is_empty()
    checks['single_valid_polygon_'+str(layer_key)] = after_local[layer_key].size() == 1

joint_rows = []
for index,inst in enumerate(instances,1):
    row = {'index':index}
    main = inst.dcplx_trans * pya.DPoint(-200,-1.95)
    row['active_main_center_um'] = [main.x,main.y]
    for layer_key,width in (((20,1),7.2),((10,2),17.2)):
        li = layout.layer(*layer_key)
        region = pya.Region(top.begin_shapes_rec(li)).merged()
        # 取接头两侧各1nm的完整截面，确认包络没有台阶或缺口。
        local = pya.DBox(-200.001,-1.95-width/2,-199.999,-1.95+width/2)
        test = pya.Region(pya.DPolygon(local).transformed(inst.dcplx_trans).to_itype(layout.dbu))
        covered = (test-region).is_empty()
        row['layer_'+str(layer_key)+'_full_width_connected'] = covered
        checks['join_'+str(index)+'_'+str(layer_key)] = covered
    joint_rows.append(row)

core = pya.Region(top.begin_shapes_rec(layout.layer(20,0))).merged()
clad = pya.Region(top.begin_shapes_rec(layout.layer(20,1))).merged()
trench = pya.Region(top.begin_shapes_rec(layout.layer(10,2))).merged()
etch1 = clad-core
markers = {'core_width_030':core.width_check(300).size(),
           'core_space_030':core.space_check(300).size(),
           'etch1_width_030':etch1.width_check(300).size(),
           'etch1_space_030':etch1.space_check(300).size(),
           'sin_trench_width_020':trench.width_check(200).size(),
           'sin_trench_space_020':trench.space_check(200).size()}
if not all(checks.values()):
    raise RuntimeError('四处包络对齐检查失败：'+json.dumps(checks))
layout.write(payload['output'])
{'top_cell':top.name, 'mux_cell':mux.name, 'instances_updated':4,
 'transition_length_um':200.0, 'mux_body_center_um':0.0, 'active_end_center_um':-1.95,
 'layers_aligned':['20/1','10/2'], 'joint_checks':joint_rows, 'checks':checks,
 'manual_rule_markers_after':markers, 'optical_validation_pass':False,
 'note':'仅平滑对齐刻蚀包络；1.33um LN1脊、电极、圆弯、黑盒和21层均未改变。'}
'''


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=BASE)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    vertices = {'(20, 1)': shifted_envelope_vertices(7.2),
                '(10, 2)': shifted_envelope_vertices(17.2)}
    for width,key in ((7.2,'(20, 1)'),(17.2,'(10, 2)')):
        shape = Polygon(vertices[key])
        if not shape.is_valid or not box(-200,-12,950,12).covers(shape):
            raise ValueError('生成的包络多边形无效：'+key)
        for x,center in ((-200,-1.95),(0,0),(950,0)):
            section = shape.intersection(__import__('shapely').geometry.LineString([(x,-20),(x,20)]))
            if abs(section.length-width)>1e-8 or abs((section.bounds[1]+section.bounds[3])/2-center)>1e-8:
                raise ValueError('包络端点宽度或中心不正确：'+key)
    args.output_dir.mkdir(parents=True, exist_ok=False)
    output = (args.output_dir/'四程10GHz_15mm_MUX包络四处对齐_待验证工作稿.gds').resolve()
    before = hashlib.sha256(args.input.read_bytes()).hexdigest()
    payload = {'input':str(args.input.resolve()), 'output':str(output), 'vertices':vertices,
               'widths':{'(20, 1)':7.2,'(10, 2)':17.2}}
    client = KLinkClient()
    client.connect()
    try:
        response = client.call('exec.python', {'code':'payload = '+repr(payload)+'\n'+KLAYOUT_CODE})
        if response.get('exception'):
            raise RuntimeError(str(response['exception']))
        report = response['return_value']
        report['file_info'] = client.call('layout.file_info', {'path':str(output),'detail':'counts'})
        client.call('layout.show_file', {'path':str(output),'mode':'new'})
        client.call('view.show_cell', {'cell':report['top_cell'],'zoom_fit':True})
        for layer,color in [('20/0','#C62858'),('20/1','#339A90'),('10/2','#65A7CF'),('42/0','#D5A12C')]:
            client.call('layer.set_style', {'layer':layer,'fill_color':color,'frame_color':color,
                        'dither_pattern':0 if layer in ('20/0','42/0') else 1,'line_width':1})
        client.call('layer.set_visible', {'layers':['101/0','100/0','100/30','56/30'],'visible':False})
        client.call('view.hier_levels', {'min':0,'max':10})
        for name,bbox in [('左侧上下两处包络对齐.png',[1750,-44,2050,44]),
                          ('右侧上下两处包络对齐.png',[16950,-44,17250,44])]:
            client.call('view.screenshot', {'mode':'path','path':str((args.output_dir/name).resolve()),
                        'bbox_um':bbox,'width_px':1800,'height_px':620})
        client.call('view.zoom_box', {'bbox_um':[1750,20,2050,42]})
    finally:
        client.close()
    if hashlib.sha256(args.input.read_bytes()).hexdigest()!=before:
        raise RuntimeError('v6源文件哈希发生变化')
    report.update({'working_gds':str(output),'source_gds':str(args.input.resolve()),
                   'source_sha256':before,'source_unchanged':True,'foundry_signoff':False,
                   'remaining':['包络过渡的光学影响和完整MUX级联未验证',
                                'MUX本体辅助脊仍有0.3um法向宽度/格点标记',
                                '未处理21层、欧拉弯、RF焊盘或终端；继承v6其余限制']})
    (args.output_dir/'四处包络对齐检查.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k not in ('file_info',)},ensure_ascii=False,indent=2))


if __name__ == '__main__':
    main()
