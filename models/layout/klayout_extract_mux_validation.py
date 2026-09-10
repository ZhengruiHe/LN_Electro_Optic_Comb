"""只读提取MUX实际绘制层；供仿真/版图一致性检查，不保存GDS。"""
import hashlib
import json
from pathlib import Path
import pya

source=Path(gds_path)
layout=pya.Layout();layout.read(str(source))
cell=layout.cell(mux_name)
if cell is None:raise ValueError('实际GDS缺少指定MUX')
layers=[]
for key in [(20,0),(20,1),(10,2),(21,0),(21,1),(21,2)]:
    region=pya.Region(cell.begin_shapes_rec(layout.layer(*key))).merged()
    layers.append({'layer':list(key),'polygons':[
        {'exterior':[[pt.x*layout.dbu,pt.y*layout.dbu] for pt in poly.each_point_hull()],
         'holes':[[[pt.x*layout.dbu,pt.y*layout.dbu] for pt in poly.each_point_hole(h)] for h in range(poly.holes())]}
        for poly in region.each()]})
result={'source_gds':str(source),'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
        'cell':cell.name,'dbu_um':layout.dbu,'layers':layers}
with open(report_path,'x',encoding='utf-8') as f:json.dump(result,f,ensure_ascii=False,indent=2)
