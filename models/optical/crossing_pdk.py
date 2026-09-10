"""自定义交叉器的PDK前置检查；仅使用本地资料，不是官方签核DRC。

按流片手册第2.2/4.3节检查，而不是IPKISS通用的0.15/0.20 um默认值。
同时反向验证两次刻蚀掩膜能否得到400/200/0 nm三个LN高度区域。
"""
from __future__ import annotations

import argparse
import hashlib
from itertools import product
import json
import math
import subprocess
from pathlib import Path
from types import SimpleNamespace

from shapely.geometry import box

ROOT = Path(__file__).resolve().parents[2]
RULE_VERSION = "manual_ln_030_030_two_etch_v1"


def validate_nominal(args, cfg, core, slab):
    """轻量前置门槛；不改变仿真几何，失败时禁止启动求解器。"""
    slab_width = getattr(args, "slab_center_width_um", 7.2)
    values = [args.width_x_um, args.width_y_um, args.port_width_um, args.half_length_um, slab_width]
    if not all(math.isfinite(v) and v > 0 for v in values):
        raise ValueError("交叉器尺寸必须为有限正数")
    if min(values[:3]) < .3:
        raise ValueError("违反手册LN最小线宽0.30 um")
    if slab_width < .3:
        raise ValueError("LN2平台违反手册最小线宽0.30 um")
    if args.reference_axis and slab_width != 7.2:
        raise ValueError("同长直波导基准必须使用7.2um等宽LN2平台")
    if min(values[:2]) < args.port_width_um:
        raise ValueError("当前已验证参数族为渐扩交叉器，中心不得比端口窄")
    stack, wg = cfg["stack_um"], cfg["waveguide_um"]
    if not (math.isclose(stack["ln"], .4) and wg["etch_depths"] == [.2]
            and math.isclose(wg["sidewall_angle_deg_from_horizontal"], 70)
            and wg["active_region_sin_fully_removed"]
            and math.isclose(wg["ln_isolation_width"], 7.2)
            and math.isclose(stack["sin"] + stack["ln_sin_interlayer_oxide"], .7)
            and math.isclose(stack["top_cladding"], 1.0)):
        raise ValueError("交叉器固定截面与项目名义工艺配置不一致；不能静默忽略配置变化")
    for shape in (core, slab):
        if not shape.is_valid or shape.geom_type != "Polygon" or shape.interiors:
            raise ValueError("当前交叉器必须是有效、无孔、连通图形")
    # 每200nm台阶在单侧横向外扩72.794nm，不能把顶宽当底宽。
    offset = .2 / math.tan(math.radians(70))
    roi = box(-args.half_length_um-3, -args.half_length_um-3,
              args.half_length_um+3, args.half_length_um+3)
    bottom = core.buffer(offset, join_style=2).intersection(roi)
    if bottom.difference(slab).area > 1e-9:
        raise ValueError("LN1脊底面超出LN2残余平台；两次套刻截面不可按当前模型实现")
    return {"rule_version": RULE_VERSION, "parameter_precheck_pass": True,
            "ln_min_line_um": .3, "ln_min_space_um": .3,
            "sidewall_lateral_offset_per_200nm_um": offset,
            "port_top_width_um": args.port_width_um,
            "port_rib_bottom_width_um": args.port_width_um+2*offset,
            "ln2_center_top_width_um": slab_width,
            "ln2_port_top_width_um": 7.2,
            "provenance": {"400nm_LN": "本地PDK手册典型值",
                           "200nm_plus_200nm_etch": "用户转述的代工工艺说明",
                           "70deg_both_steps": "用户指定名义值，与水平面夹角",
                           "top_cladding_1um": "当前项目名义假设，非手册保证值",
                           "no_SiN_below_LN": "手册4.4节；SiO2回填按用户确认"},
            "foundry_signoff": False}


def main():
    from crossing_fdtd import GEOMETRY_VERSION, geometry
    from mode_pdk_sweep import DEFAULT_CONFIG, load_config
    from klink import KLinkClient
    p = argparse.ArgumentParser()
    p.add_argument("--widths", default="1.5,2,2.5,3")
    p.add_argument("--half-lengths", default="12,8,16")
    p.add_argument("--slab-widths", default="7.2", help="LN2单臂中心顶宽列表，端口平台固定7.2um")
    p.add_argument("--slab-flat-half-length-um", type=float, default=6.4)
    p.add_argument("--klayout-exe", type=Path, help="使用已安装KLayout批处理，不要求打开GUI或KLink")
    p.add_argument("--output-dir", type=Path, required=True)
    a = p.parse_args()
    a.output_dir.mkdir(parents=True, exist_ok=True)
    report_path = a.output_dir / "PDK前置检查.json"
    if report_path.exists():
        raise FileExistsError("请使用新输出目录，保留历史PDK核查报告")
    cfg = load_config(DEFAULT_CONFIG)
    cases = []
    for slab_width, length, width in product(map(float,a.slab_widths.split(',')),
                                             map(float,a.half_lengths.split(',')),
                                             map(float,a.widths.split(','))):
        args = SimpleNamespace(width_x_um=width, width_y_um=width,
                               port_width_um=.7, half_length_um=length, reference_axis=None,
                               slab_center_width_um=slab_width,
                               slab_flat_half_length_um=a.slab_flat_half_length_um)
        extent = length+6
        core, slab = geometry(args, extent)
        precheck = validate_nominal(args, cfg, core, slab)
        case = {"width_x_um": width, "width_y_um": width, "half_length_um": length,
                "port_width_um": .7, "slab_center_width_um": slab_width, "precheck": precheck,
                "geometry_version": GEOMETRY_VERSION,
                "slab_flat_half_length_um": a.slab_flat_half_length_um,
                "core": list(core.exterior.coords)[:-1],
                "slab": list(slab.exterior.coords)[:-1], "extent": extent,
                "cell": (f"CrossW{width:g}L{length:g}S{slab_width:g}").replace('.', 'p')}
        cases.append(case)
    output = (a.output_dir / "交叉器候选_两次刻蚀映射核查.gds").resolve()
    payload = {"cases": cases, "output": str(output)}
    code = "payload = " + repr(payload) + "\n" + '''
audit_layout = pya.Layout()
audit_layout.dbu = 0.001
audit_top = audit_layout.create_cell("CrossingPDKAudit")
reports = []
for k, case in enumerate(payload["cases"]):
    cell = audit_layout.create_cell(case["cell"])
    audit_top.insert(pya.CellInstArray(cell.cell_index(), pya.Trans((k%4)*70000, (k//4)*70000)))
    def region(points):
        return pya.Region(pya.DPolygon([pya.DPoint(*xy) for xy in points]).to_itype(audit_layout.dbu))
    core, slab = region(case["core"]), region(case["slab"])
    e = case["extent"]
    window = pya.Region(pya.DBox(-e,-e,e,e).to_itype(audit_layout.dbu))
    # 两次套刻：第一遍窗口减脊；第二遍窗口减残余平台。
    # 20/1本身不是材料平台，必须经过手册布尔逻辑还原。
    etch1, etch2 = window-core, window-slab
    layers = {(20,0):core, (20,1):window, (21,0):slab, (21,1):window, (10,2):window}
    for layer, reg in layers.items():
        cell.shapes(audit_layout.layer(*layer)).insert(reg)
    checks = {}
    for name, reg in (("LN1_core",core),("LN2_keep",slab),("LN1_etch",etch1),("LN2_etch",etch2)):
        checks[name+"_width_030"] = reg.width_check(300).size()
        checks[name+"_space_030"] = reg.space_check(300).size()
    checks["etch2_outside_etch1"] = (etch2-etch1).size()
    checks["reconstructed_400nm_mismatch"] = ((window-etch1)^core).size()
    checks["reconstructed_200nm_mismatch"] = ((etch1-etch2)^(slab-core)).size()
    checks["reconstructed_0nm_mismatch"] = (etch2^(window-slab)).size()
    checks["core_outside_slab"] = (core-slab).size()
    reports.append({"cell":case["cell"],"checks":checks,"known_geometry_rules_pass":all(v==0 for v in checks.values())})
audit_layout.write(payload["output"])
reports
'''
    if a.klayout_exe:
        if not a.klayout_exe.is_file():
            raise FileNotFoundError(a.klayout_exe)
        numeric_path = (a.output_dir / "规则数值.json").resolve()
        job_path = (a.output_dir / "生成并检查掩膜.py").resolve()
        job_path.write_text("import json\nimport pya\n" + code +
                            "\nwith open(" + repr(str(numeric_path)) + ", 'w', encoding='utf-8') as stream:\n"
                            "    json.dump(reports, stream, ensure_ascii=False, indent=2)\n", encoding="utf-8")
        subprocess.run([str(a.klayout_exe.resolve()), "-b", "-r", str(job_path)],
                       check=True, creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
        response = {"return_value": json.loads(numeric_path.read_text(encoding="utf-8"))}
    else:
        client = KLinkClient()
        client.connect()
        try:
            response = client.call("exec.python", {"code": code})
            if response.get("exception"):
                raise RuntimeError(str(response["exception"]))
        finally:
            client.close()
    for case, checks in zip(cases, response["return_value"], strict=True):
        case.update(checks)
        case.pop("core"); case.pop("slab")
    report = {"rule_version": RULE_VERSION, "cases": cases,
              "all_known_geometry_rules_pass": all(c["known_geometry_rules_pass"] for c in cases),
              "gds": str(output), "gds_sha256": hashlib.sha256(output.read_bytes()).hexdigest(),
              "manual": "南智光电流片说明-TFLN-on-SiN-V0.2试用版.pdf；2.2节及4.3节",
              "not_checked": ["完整芯片官方DRC、黑盒禁布与端口连接", "真实工艺掩膜偏置与套刻容差",
                              "局部仿真域之外未刻LN/SiN材料对光场的影响", "完整光学级联"],
              "foundry_signoff": False,
              "note": "独立候选掩膜映射示例；不得用旧版20/1包络直接假定LN2刻除已经实现"}
    report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({"report": str(report_path), "all_known_geometry_rules_pass": report["all_known_geometry_rules_pass"],
                      "cases": [{"cell": c["cell"], "checks": c["checks"]} for c in cases]}, ensure_ascii=False, indent=2))
    if not report["all_known_geometry_rules_pass"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
