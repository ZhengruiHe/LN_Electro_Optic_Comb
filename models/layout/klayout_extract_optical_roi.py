"""KLayout只读提取实际GDS芯层的局部区域；不重画、不保存源GDS。"""
import hashlib
import json
from pathlib import Path
import pya


with open(payload_file, encoding="utf-8") as stream:
    payload = json.load(stream)
source = Path(payload["gds"])
before = hashlib.sha256(source.read_bytes()).hexdigest()
layout = pya.Layout()
layout.read(str(source))
reports = []
for job in payload["jobs"]:
    cell = layout.cell(job["cell"])
    if cell is None:
        raise ValueError("源GDS不存在单元：" + job["cell"])
    region = pya.Region(cell.begin_shapes_rec(layout.layer(20, 0))).merged()
    clip = pya.Region(pya.DBox(*job["bounds_um"]).to_itype(layout.dbu))
    region = (region & clip).merged()
    reports.append({
        "name": job["name"], "cell": cell.name, "bounds_um": job["bounds_um"],
        "core_polygons": [{
            "exterior": [[pt.x * layout.dbu, pt.y * layout.dbu] for pt in poly.each_point_hull()],
            "holes": [[[pt.x * layout.dbu, pt.y * layout.dbu] for pt in poly.each_point_hole(h)]
                      for h in range(poly.holes())],
        } for poly in region.each()],
    })
if hashlib.sha256(source.read_bytes()).hexdigest() != before:
    raise RuntimeError("只读提取期间GDS发生变化")
with open(payload["output"], "x", encoding="utf-8") as stream:
    json.dump({"gds": str(source), "sha256": before, "dbu_um": layout.dbu,
               "source_unchanged": True, "jobs": reports}, stream, ensure_ascii=False, indent=2)
