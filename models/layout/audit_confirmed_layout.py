"""KLayout批处理只读复核已有四程GDS，不生成或修改几何。

参数：-rd source_gds=... -rd baseline_json=... -rd report_json=新报告.json
LN2/21层仅记录，不作为本轮判据；原有MUX标记不会被豁免或清除。
"""
import hashlib
import json
from pathlib import Path
import pya


source = Path(source_gds)
baseline = json.loads(Path(baseline_json).read_text(encoding="utf-8"))
digest = hashlib.sha256(source.read_bytes()).hexdigest()
layout = pya.Layout()
layout.read(str(source))
top = layout.top_cell()
dbu = layout.dbu


def region(key):
    return pya.Region(top.begin_shapes_rec(layout.layer(*key))).merged()


def box(bounds):
    return pya.Region(pya.DBox(*bounds).to_itype(dbu))


core, window, trench, metal = [region(key) for key in
                             [(20, 0), (20, 1), (10, 2), (42, 0)]]
etch = (window - core).merged()
markers = {
    "core_width_030": core.width_check(round(.30/dbu)).size(),
    "core_space_030": core.space_check(round(.30/dbu)).size(),
    "etch1_width_030": etch.width_check(round(.30/dbu)).size(),
    "etch1_space_030": etch.space_check(round(.30/dbu)).size(),
    "sin_width_020": trench.width_check(round(.20/dbu)).size(),
    "sin_space_020": trench.space_check(round(.20/dbu)).size(),
    "metal_width_2": metal.width_check(round(2/dbu)).size(),
    "metal_space_3": metal.space_check(round(3/dbu)).size(),
    "core_not_in_clad": (core-window).size(),
    "core_not_in_trench": (core-trench).size(),
}
conductors = [pya.Region(poly) for poly in metal.each()]
gsg = baseline["source_report"]
landings = {}
for side, points in gsg["probe_contact_centers_um"].items():
    for net, (x, y) in points.items():
        landing = box([x-25, y-25, x+25, y+25])
        landings[f"{side}_{net}"] = [i for i, conductor in enumerate(conductors)
                                     if (landing-conductor).is_empty()]
net_ids = {}
for net, y in [("S", 0), ("Gupper", 88.5), ("Glower", -88.5)]:
    middle = box([9499, y-1, 9501, y+1])
    net_ids[net] = [i for i, conductor in enumerate(conductors)
                    if (middle-conductor).is_empty()]
connections = {net: len(ids) == 1 and all(landings[f"{side}_{net}"] == ids
               for side in ("RF_IN", "RF_OUT")) for net, ids in net_ids.items()}
active_core = {str(y): (box([2000, y-.665, 17000, y+.665])-core).is_empty()
               for y in (30., -30.)}
flags = {
    "source_matches_last_verified_sha256": digest == baseline["output_sha256"],
    "one_top_and_1nm_grid": len(layout.top_cells()) == 1 and abs(dbu-.001) < 1e-12,
    "actual_markers_unchanged": markers == baseline["actual_saved_gds_markers"],
    "three_metal_nets": len(conductors) == 3,
    "both_gsg_ends_connected": all(connections.values()),
    "active_two_15mm_1p33um_ribs_present": all(active_core.values()),
    "no_metal_core_planar_overlap": (core & metal).is_empty(),
}
result = {
    "scope": "本轮采用版复核，非流片签核；LN2/21层暂缓",
    "source_gds": str(source), "sha256": digest,
    "top_cell": top.name, "bbox_um": str(top.dbbox()),
    "checks": flags, "nominal_version_integrity_pass": all(flags.values()),
    "manual_rule_markers_excluding_ln2": markers,
    "gsg_net_connections": connections, "gsg_landing_conductor_ids": landings,
    "active_core_rectangles": active_core,
    "ln2_present_layers": [str(layout.get_info(li)) for li in layout.layer_indices()
                           if layout.get_info(li).layer == 21 and not top.dbbox(li).empty()],
    "remaining_not_waived": ["MUX原有40条线宽及40条刻蚀间距标记", "GSG接口S参数尚未求解"],
    "optical_bend_decision": "采用已得90度粗网格时间窗初筛；按用户要求不再补查，不等同于全部弯收敛通过",
    "source_unchanged": hashlib.sha256(source.read_bytes()).hexdigest() == digest,
}
with Path(report_json).open("x", encoding="utf-8") as stream:
    json.dump(result, stream, ensure_ascii=False, indent=2)
print(json.dumps(result, ensure_ascii=False, indent=2))
if not all(flags.values()):
    raise RuntimeError("存在需要说明的版本或连接差异；未修改GDS")
