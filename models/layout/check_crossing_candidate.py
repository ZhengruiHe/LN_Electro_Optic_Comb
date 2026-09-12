"""复核旧输入交叉，并生成有意正交交叉的中心线候选。

不写 GDS、不修改 PDK、不把黑盒当成已通过的光学器件。
保持历史回路中心线长度；实际交叉器/锥形/转向的群时延仍需回填。
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from generate_four_pass_gds import append_arc, polyline_length, return_route, semicircle_chord_length
from route_geometry import check_route_network

ROOT = Path(__file__).resolve().parents[2]
BASELINE = ROOT / "results/layout/四程电光梳_10GHz_15mm_3T_3.5T_3T_名义_检查.json"


def orthogonal_loop2(length: float) -> list[list[float]]:
    xr, xl, ys, ye = 18150.0, 850.0, 26.6, -30.15
    launch, r0, rm, rf = 400.0, 270.0, 80.0, 80.0
    xc, xt = 300.0, 2150.0
    lane1 = ys + 2 * r0
    lane2, lane3 = lane1 + 2 * rm, lane1 + 4 * rm
    arc_total = sum(semicircle_chord_length(r) for r in (r0, rm, rm))
    # 两个 90 度弯各取 32 段，与 180 度 64 段离散一致。
    arc_total += semicircle_chord_length(rf)
    vertical = lane3 - ye - 2 * rf
    fixed = xr - xl + 2 * launch + 2 * (xl - xc - rf) + vertical + arc_total
    backtrack = (length - fixed) / 2
    if backtrack <= 0 or xt + backtrack + rm >= xr + launch - r0:
        raise ValueError("当前正交回路无法容纳指定长度")
    p = [[xr, ys], [xr + launch, ys]]
    append_arc(p, xr + launch, ys + r0, r0, -math.pi / 2, math.pi / 2)
    p.append([xt, lane1])
    append_arc(p, xt, lane1 + rm, rm, -math.pi / 2, -3 * math.pi / 2)
    p.append([xt + backtrack, lane2])
    append_arc(p, xt + backtrack, lane2 + rm, rm, -math.pi / 2, math.pi / 2)
    p.append([xc + rf, lane3])
    append_arc(p, xc + rf, lane3 - rf, rf, math.pi / 2, math.pi, segments=32)
    p.append([xc, ye + rf])
    append_arc(p, xc + rf, ye + rf, rf, math.pi, 3 * math.pi / 2, segments=32)
    p.append([xl, ye])
    if not math.isclose(polyline_length(p), length, abs_tol=1e-6):
        raise RuntimeError("正交回路长度回标失败")
    return p


def candidate_input() -> list[list[float]]:
    p = [[0.0, 250.0], [445.0, 250.0]]
    append_arc(p, 445, 170, 80, math.pi / 2, 0, segments=32)
    p.append([525.0, 110.15])
    append_arc(p, 605, 110.15, 80, math.pi, 3 * math.pi / 2, segments=32)
    p.append([850.0, 30.15])
    return p


def build(lengths: list[float], candidate: bool) -> dict:
    common = dict(x_right=18150, x_left=850, start_radius=80)
    loop1, _ = return_route(**common, y_start=29.85, y_end=33.4,
                            target_length=lengths[0], outward_sign=1)
    loop3, _ = return_route(**common, y_start=-29.85, y_end=-33.4,
                            target_length=lengths[2], outward_sign=-1)
    if candidate:
        loop2 = orthogonal_loop2(lengths[1])
        incoming = candidate_input()
    else:
        loop2, _ = return_route(x_right=18150, x_left=850, y_start=26.6, y_end=-30.15,
                               target_length=lengths[1], outward_sign=1,
                               start_radius=250, launch_extension=400)
        incoming = [[300, 30.15], [850, 30.15]]
    return {"loop1": loop1, "loop2": loop2, "loop3": loop3,
            "input": incoming, "output": [[18150, -26.6], [18700, -26.6]]}


def svg_preview(routes: dict, crossing: dict) -> str:
    # 非掩膜图：局部横纵同倍率，箭头标识正交穿行，空白框明确为黑盒预留。
    colors = {"loop1": "#cc802a", "loop2": "#306ab4", "loop3": "#46a176",
              "input": "#c22f58", "output": "#8055aa"}
    pieces = ['<svg xmlns="http://www.w3.org/2000/svg" width="1200" height="980" viewBox="0 0 1200 980">',
              '<rect width="1200" height="980" fill="white"/>',
              '<style>text{font-family:Microsoft YaHei,Arial,sans-serif;font-size:18px}</style>',
              '<text x="35" y="35">正交交叉候选：左侧光路中心线（不是可流片掩膜）</text>',
              '<svg x="45" y="65" width="1110" height="775" viewBox="-120 -930 1200 1100">',
              '<g transform="scale(1,-1)" fill="none">']
    for name, points in routes.items():
        xy = " ".join(f"{x:.4f},{y:.4f}" for x, y in points)
        pieces.append(f'<polyline points="{xy}" stroke="{colors[name]}" stroke-width="2.6"/>')
    pieces += ['<rect x="277.5" y="227.5" width="45" height="45" fill="white" stroke="#171717" stroke-width="2"/>',
               '<path d="M280,250 h40 M300,230 v40" stroke="#777" stroke-width="1.8" stroke-dasharray="3 2"/>',
               '<path d="M255,250 l-9,5 v-10 z" fill="#c22f58"/>',
               '<path d="M300,295 l-5,9 h10 z" fill="#306ab4"/>', '</g>',
               '<text x="365" y="-255">45 × 45 µm 交叉器预留</text>',
               '<text x="365" y="-230">横向：光输入 →；纵向：回路 2 ↓</text>',
               '<text x="680" y="-8">上方输入 / MUX 端口</text>',
               '<text x="680" y="65">下方 MUX 端口</text>', '</svg>']
    for index, (name, label) in enumerate((("input", "光输入"), ("loop1", "回路 1"),
                                          ("loop2", "回路 2"), ("loop3", "回路 3"))):
        x = 50 + index * 235
        pieces.append(f'<path d="M{x},870 h35" stroke="{colors[name]}" stroke-width="3"/>'
                      f'<text x="{x+45}" y="876">{label}</text>')
    pieces += ['<text x="45" y="920">检查范围：输入、输出、三条回路；器件内部及全层 DRC 另行检查。</text>',
               '<text x="45" y="950">0.7 ↔ 1.2 µm 过渡、交叉器内部及晶向相关损耗/时延尚未验证。</text>', '</svg>']
    return "\n".join(pieces)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline-report", type=Path, default=BASELINE)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "results/layout/_历史归档/迭代版本_20260912/正交交叉候选_v1")
    args = parser.parse_args()
    baseline = json.loads(args.baseline_report.read_text(encoding="utf-8"))
    if baseline["loop_periods"] != [3.0, 3.5, 3.0] or baseline["electrode"]["active_length_um"] != 15000:
        raise ValueError("此候选仅适用于已确认端口坐标的 15 mm、3T/3.5T/3T 基线")
    lengths = [row["actual_centerline_length_um"] for row in baseline["loop_centerline_checks"]]
    old, new = build(lengths, False), build(lengths, True)
    crossing = {"routes": ["input", "loop2"], "center_um": [300.0, 250.0], "angle_deg": 90,
                "bbox_um": [277.5, 227.5, 322.5, 272.5],
                "status": "PDK黑盒接口预留；非光学通过"}
    result = {
        "baseline_report": str(args.baseline_report.resolve()),
        "old_layout": check_route_network(old),
        "candidate": check_route_network(new, [crossing]),
        "candidate_routes_um": new,
        "length_checks": {name: {"baseline_um": polyline_length(old[name]),
                                  "candidate_um": polyline_length(new[name])}
                          for name in new},
        "crossing_port_taper_plan": {"route_width_um": 0.7, "pdk_template_width_um": 1.2,
                                      "length_each_um": 100, "status": "仅预留，未仿真"},
        "limitations": ["只通过中心线拓扑筛查时仍不能称为版图DRC通过",
                        "新增交叉器不具有本地真实掩膜或实测复S矩阵，不能采用论文损耗代替",
                        "中心线长度保持不等于真实群时延保持；交叉/渐变/晶向变化需重新回填",
                        "现有EME快速相位扫频的固定LN折射率未计入材料色散，不能用于最终时延签核",
                        "有源电极不缩短，10 GHz目标不改变，暂不做容差"],
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "光路拓扑检查.json").write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    (args.output_dir / "正交交叉_局部中心线.svg").write_text(svg_preview(new, crossing), encoding="utf-8")
    print(json.dumps({"old_unexpected": result["old_layout"]["unexpected_intersections"],
                      "candidate_pass": result["candidate"]["topology_screen_pass"],
                      "candidate_unexpected": result["candidate"]["unexpected_intersections"],
                      "output_dir": str(args.output_dir)}, ensure_ascii=False, indent=2))
    if not result["candidate"]["topology_screen_pass"]:
        raise RuntimeError("候选仍有未经设计的交叉，不能进入GDS更新")


if __name__ == "__main__":
    main()
