"""KLayout批处理只读盘点：不改源GDS，报告顶层、实例、图层和包围盒。

运行示例：klayout_app.exe -b -r 本脚本 -rd input_gds=输入.gds
可选 -rd report_json=报告.json；源图形和层级不写回。
"""
import hashlib
import json
from pathlib import Path

import pya


def box_values(box):
    if box.empty():
        return None
    return [box.left, box.bottom, box.right, box.top]


source = Path(input_gds)
layout = pya.Layout()
layout.read(str(source))
result = {
    "source": str(source),
    "sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
    "dbu_um": layout.dbu,
    "cell_count": layout.cells(),
    "tops": [],
}
for top in layout.top_cells():
    layers = []
    for li in layout.layer_indices():
        bounds = top.dbbox(li)
        if bounds.empty():
            continue
        info = layout.get_info(li)
        layers.append({"layer": info.layer, "datatype": info.datatype,
                       "name": info.name, "bbox_um": box_values(bounds),
                       "direct_shape_count": top.shapes(li).size()})
    instances = [{"cell": inst.cell.name, "transform": str(inst.dcplx_trans),
                  "bbox_um": box_values(inst.dbbox()),
                  "array_na": inst.na, "array_nb": inst.nb}
                 for inst in top.each_inst()]
    result["tops"].append({
        "name": top.name, "bbox_um": box_values(top.dbbox()),
        "layers": layers, "instances": instances,
    })
if "report_json" in globals():
    with Path(report_json).open("x", encoding="utf-8") as stream:
        json.dump(result, stream, ensure_ascii=False, indent=2)
if "geometry_json" in globals():
    if len(layout.top_cells()) != 1:
        raise ValueError("几何预览导出要求只有一个源顶层")
    top = layout.top_cell()
    preview = []
    for li in layout.layer_indices():
        info = layout.get_info(li)
        if info.layer not in (10, 20, 21, 41, 42, 56, 70, 100):
            continue
        region = pya.Region(top.begin_shapes_rec(li)).merged()
        if region.is_empty():
            continue
        preview.append({
            "layer": [info.layer, info.datatype],
            "polygons": [{
                "exterior": [[point.x * layout.dbu, point.y * layout.dbu]
                             for point in polygon.each_point_hull()],
                "holes": [[[point.x * layout.dbu, point.y * layout.dbu]
                           for point in polygon.each_point_hole(hole)]
                          for hole in range(polygon.holes())],
            } for polygon in region.each()],
        })
    with Path(geometry_json).open("x", encoding="utf-8") as stream:
        json.dump(preview, stream, ensure_ascii=False)
print(json.dumps(result, ensure_ascii=False, indent=2))
