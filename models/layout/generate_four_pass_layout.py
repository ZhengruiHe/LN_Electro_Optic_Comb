"""生成候选A四程调制器的参数化顶视版图骨架（SVG与JSON）。"""

from __future__ import annotations

import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "results" / "layout"


def mm(value_um: float) -> float:
    return value_um / 1000.0


def serpentine(
    start_x: float,
    rail_y: float,
    first_row_y: float,
    target_length_mm: float,
    direction: int,
    runs: int = 4,
    radius_mm: float = 0.08,
) -> tuple[str, float]:
    """返回同侧进出的四横段蛇形路径与用于回标的横段长度。"""
    last_row_y = first_row_y + direction * 2.0 * radius_mm * (runs - 1)
    connector_length = abs(rail_y - first_row_y) + abs(rail_y - last_row_y)
    bend_length = (runs - 1) * math.pi * radius_mm
    span = (target_length_mm - connector_length - bend_length) / runs
    if span <= 0:
        raise ValueError("目标回路长度不足以容纳指定弯曲半径")

    other_x = start_x - span if start_x > 6.0 else start_x + span
    commands = [f"M {start_x:.6f} {rail_y:.6f}", f"L {start_x:.6f} {first_row_y:.6f}"]
    current_x = start_x
    current_y = first_row_y
    for run in range(runs):
        next_x = other_x if current_x == start_x else start_x
        commands.append(f"L {next_x:.6f} {current_y:.6f}")
        current_x = next_x
        if run + 1 < runs:
            next_y = current_y + direction * 2.0 * radius_mm
            sweep = 1 if (direction > 0) == (current_x == other_x) else 0
            commands.append(
                f"A {radius_mm:.6f} {radius_mm:.6f} 0 0 {sweep} "
                f"{current_x:.6f} {next_y:.6f}"
            )
            current_y = next_y
    commands.append(f"L {start_x:.6f} {rail_y:.6f}")
    return " ".join(commands), span


def rect(x: float, y: float, width: float, height: float, css: str) -> str:
    return (
        f'<rect x="{x:.6f}" y="{y:.6f}" width="{width:.6f}" '
        f'height="{height:.6f}" class="{css}"/>'
    )


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)

    active_x0 = 1.0
    active_x1 = 11.0
    active_length = active_x1 - active_x0
    center_y = 1.60
    signal_width = mm(43.0)
    ground_width = mm(100.0)
    trunk_gap = mm(17.0)
    inner_gap = mm(5.0)
    cap_width = mm(2.0)
    neck_length = mm(4.0)
    neck_longitudinal_width = mm(10.0)
    cap_longitudinal_length = mm(45.0)
    period = mm(50.0)
    waveguide_width = mm(1.33)
    bend_radius = mm(80.0)
    target_frequency_ghz = 10.0
    preliminary_group_indices = [2.22, 2.22, 2.22]
    loop_multipliers = [2.0, 2.5, 2.0]
    speed_of_light_mm_per_s = 299_792_458.0 * 1e3
    total_loop_lengths_mm = [
        speed_of_light_mm_per_s * periods
        / (group_index * target_frequency_ghz * 1e9)
        for periods, group_index in zip(loop_multipliers, preliminary_group_indices)
    ]
    # 论文补充材料中的“回路总长度”包含上一程调制臂、模式
    # 复用器、锥形过渡、弯曲和外部回环。本骨架已单独画出10 mm
    # 调制臂，因此蛇形回环只能占用剩余的被动过渡长度。
    passive_transition_lengths_mm = [
        total_length_mm - active_length for total_length_mm in total_loop_lengths_mm
    ]
    if min(passive_transition_lengths_mm) <= 0:
        raise ValueError("回路总长度不足以容纳已画出的调制臂")

    signal_top = center_y - 0.5 * signal_width
    signal_bottom = center_y + 0.5 * signal_width
    upper_ground_bottom = signal_top - trunk_gap
    upper_ground_top = upper_ground_bottom - ground_width
    lower_ground_top = signal_bottom + trunk_gap
    upper_rail_y = center_y - (
        0.5 * signal_width + neck_length + cap_width + 0.5 * inner_gap
    )
    lower_rail_y = center_y + (
        0.5 * signal_width + neck_length + cap_width + 0.5 * inner_gap
    )

    loop1, loop1_span = serpentine(
        active_x1,
        upper_rail_y,
        1.30,
        passive_transition_lengths_mm[0],
        -1,
        radius_mm=bend_radius,
    )
    loop2, loop2_span = serpentine(
        active_x0,
        upper_rail_y,
        0.70,
        passive_transition_lengths_mm[1],
        -1,
        radius_mm=bend_radius,
    )
    # 回路2的末端从上调制区跨到下调制区，因此最后连接点替换为下光轨。
    loop2 = loop2.rsplit("L", 1)[0] + f"L {active_x0:.6f} {lower_rail_y:.6f}"
    loop3, loop3_span = serpentine(
        active_x1,
        lower_rail_y,
        1.90,
        passive_transition_lengths_mm[2],
        1,
        radius_mm=bend_radius,
    )

    metal_parts = [
        rect(active_x0, signal_top, active_length, signal_width, "signal"),
        rect(active_x0, upper_ground_top, active_length, ground_width, "ground"),
        rect(active_x0, lower_ground_top, active_length, ground_width, "ground"),
    ]
    x = active_x0 + 0.5 * period
    while x < active_x1:
        # 信号电极上下T形支节。
        metal_parts.extend(
            [
                rect(
                    x - 0.5 * neck_longitudinal_width,
                    signal_top - neck_length,
                    neck_longitudinal_width,
                    neck_length,
                    "signal",
                ),
                rect(
                    x - 0.5 * cap_longitudinal_length,
                    signal_top - neck_length - cap_width,
                    cap_longitudinal_length,
                    cap_width,
                    "signal",
                ),
                rect(
                    x - 0.5 * neck_longitudinal_width,
                    signal_bottom,
                    neck_longitudinal_width,
                    neck_length,
                    "signal",
                ),
                rect(
                    x - 0.5 * cap_longitudinal_length,
                    signal_bottom + neck_length,
                    cap_longitudinal_length,
                    cap_width,
                    "signal",
                ),
                # 两侧地电极向内的T形支节。
                rect(
                    x - 0.5 * neck_longitudinal_width,
                    upper_ground_bottom,
                    neck_longitudinal_width,
                    neck_length,
                    "ground",
                ),
                rect(
                    x - 0.5 * cap_longitudinal_length,
                    upper_ground_bottom + neck_length,
                    cap_longitudinal_length,
                    cap_width,
                    "ground",
                ),
                rect(
                    x - 0.5 * neck_longitudinal_width,
                    lower_ground_top - neck_length,
                    neck_longitudinal_width,
                    neck_length,
                    "ground",
                ),
                rect(
                    x - 0.5 * cap_longitudinal_length,
                    lower_ground_top - neck_length - cap_width,
                    cap_longitudinal_length,
                    cap_width,
                    "ground",
                ),
            ]
        )
        x += period

    svg = f'''<svg xmlns="http://www.w3.org/2000/svg" width="1400" height="470" viewBox="0 0 12 3.2">
  <defs>
    <marker id="arrow1" markerWidth="0.10" markerHeight="0.10" refX="0.08" refY="0.05" orient="auto"><path d="M0,0 L0.10,0.05 L0,0.10 Z" class="fill-p1"/></marker>
    <marker id="arrow2" markerWidth="0.10" markerHeight="0.10" refX="0.08" refY="0.05" orient="auto"><path d="M0,0 L0.10,0.05 L0,0.10 Z" class="fill-p2"/></marker>
    <marker id="arrow3" markerWidth="0.10" markerHeight="0.10" refX="0.08" refY="0.05" orient="auto"><path d="M0,0 L0.10,0.05 L0,0.10 Z" class="fill-p3"/></marker>
    <marker id="arrow4" markerWidth="0.10" markerHeight="0.10" refX="0.08" refY="0.05" orient="auto"><path d="M0,0 L0.10,0.05 L0,0.10 Z" class="fill-p4"/></marker>
    <style>
      .die {{ fill:#f7f8fb; stroke:#27324a; stroke-width:.012; }}
      .signal {{ fill:#c58a1b; }} .ground {{ fill:#8a6b32; }}
      .p1,.p2,.p3,.p4,.loop {{ fill:none; stroke-width:.020; }}
      .p1 {{ stroke:#d12c58; }} .p2 {{ stroke:#6548c7; }} .p3 {{ stroke:#16845b; }} .p4 {{ stroke:#d56b1f; }}
      .loop {{ stroke:#526177; stroke-width:.012; }}
      .fill-p1 {{ fill:#d12c58; }} .fill-p2 {{ fill:#6548c7; }} .fill-p3 {{ fill:#16845b; }} .fill-p4 {{ fill:#d56b1f; }}
      .rail {{ fill:#d12c58; opacity:.18; }}
      .mux {{ fill:#d8deea; stroke:#526177; stroke-width:.008; }}
      .label {{ font: .12px sans-serif; fill:#172033; }}
      .small {{ font: .09px sans-serif; fill:#34425c; }}
      .title {{ font: 600 .16px sans-serif; fill:#172033; }}
    </style>
  </defs>
  <rect x="0.02" y="0.02" width="11.96" height="3.16" rx=".04" class="die"/>
  <text x="0.18" y="3.03" class="title">候选A四程T形行波调制器：按12 mm × 3.2 mm版图骨架显示</text>
  <text x="0.18" y="3.17" class="small">金属有源区10 mm；LN最小弯曲半径80 µm；回路长度按10 GHz群时延回标</text>
  {''.join(metal_parts)}
  <rect x="{active_x0}" y="{upper_rail_y - waveguide_width/2:.6f}" width="{active_length}" height="{waveguide_width:.6f}" class="rail"/>
  <rect x="{active_x0}" y="{lower_rail_y - waveguide_width/2:.6f}" width="{active_length}" height="{waveguide_width:.6f}" class="rail"/>
  <path d="M {active_x0} {upper_rail_y} L {active_x1} {upper_rail_y}" class="p1" marker-end="url(#arrow1)"/>
  <path d="{loop1}" class="loop"/>
  <path d="M {active_x1} {upper_rail_y + .008} L {active_x0} {upper_rail_y + .008}" class="p2" marker-end="url(#arrow2)"/>
  <path d="{loop2}" class="loop"/>
  <path d="M {active_x0} {lower_rail_y} L {active_x1} {lower_rail_y}" class="p3" marker-end="url(#arrow3)"/>
  <path d="{loop3}" class="loop"/>
  <path d="M {active_x1} {lower_rail_y + .008} L {active_x0} {lower_rail_y + .008}" class="p4" marker-end="url(#arrow4)"/>
  <rect x="10.82" y="1.245" width=".18" height=".12" rx=".02" class="mux"/><text x="10.825" y="1.235" class="small">MUX1</text>
  <rect x="0.82" y="1.30" width=".18" height=".12" rx=".02" class="mux"/><text x="0.825" y="1.29" class="small">MUX2</text>
  <rect x="10.82" y="1.79" width=".18" height=".12" rx=".02" class="mux"/><text x="10.825" y="1.78" class="small">MUX3</text>
  <text x="5.05" y="1.51" class="small">上调制区：第1程TE0 →　第2程TE1 ←</text>
  <text x="5.05" y="1.82" class="small">下调制区：第3程TE0 →　第4程TE1 ←</text>
  <text x="4.45" y="1.27" class="small">回路1总长：{total_loop_lengths_mm[0]:.3f} mm（{loop_multipliers[0]:g}T，待实际ng回标）</text>
  <text x="1.20" y="0.13" class="small">回路2总长：{total_loop_lengths_mm[1]:.3f} mm（{loop_multipliers[1]:g}T，补偿π极性差）</text>
  <text x="6.10" y="2.46" class="small">回路3总长：{total_loop_lengths_mm[2]:.3f} mm（{loop_multipliers[2]:g}T，待实际ng回标）</text>
  <text x="1.02" y="1.70" class="small">光输入</text><text x="0.36" y="1.92" class="small">光输出</text>
</svg>'''

    svg_path = OUTPUT / "候选A四程实际结构_10GHz.svg"
    svg_path.write_text(svg, encoding="utf-8")

    metadata = {
        "status": "10ghz_layout_skeleton_pending_a70_group_delay_not_mask_ready",
        "target_frequency_ghz": target_frequency_ghz,
        "footprint_mm": [12.0, 3.2],
        "active_electrode_length_mm": active_length,
        "pdk_ln_min_bend_radius_mm": bend_radius,
        "upper_rail_y_mm": upper_rail_y,
        "lower_rail_y_mm": lower_rail_y,
        "four_passes": [
            {"pass": 1, "mode": "TE0", "region": "upper", "direction": "left_to_right"},
            {"pass": 2, "mode": "TE1", "region": "upper", "direction": "right_to_left"},
            {"pass": 3, "mode": "TE0", "region": "lower", "direction": "left_to_right"},
            {"pass": 4, "mode": "TE1", "region": "lower", "direction": "right_to_left"},
        ],
        "loop_period_multipliers": loop_multipliers,
        "preliminary_loop_group_indices": preliminary_group_indices,
        "total_loop_lengths_mm": total_loop_lengths_mm,
        "passive_transition_lengths_after_subtracting_10mm_modulation_arm_mm": passive_transition_lengths_mm,
        "serpentine_horizontal_spans_mm": [loop1_span, loop2_span, loop3_span],
        "electrode_um": {
            "signal_width": 60.0,
            "ground_width": 100.0,
            "inner_gap": 5.0,
            "trunk_gap": 21.0,
            "t_cap_width": 2.0,
            "t_neck_length": 6.0,
            "t_neck_longitudinal_width": 2.0,
            "t_cap_longitudinal_length": 20.0,
            "period": 50.0,
        },
        "limitations": [
            "Each total loop length includes the preceding 10 mm modulation arm; the drawn serpentine uses only the remaining passive transition length",
            "TE0/TE1 multiplexers and crossing are placeholders pending EME/FDTD design",
            "serpentine corner and transition lengths require final routed path integration",
            "M1 pads, termination, optical couplers, and official PDK DRC are not yet included",
        ],
    }
    json_path = OUTPUT / "候选A四程实际结构参数_10GHz.json"
    json_path.write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"svg={svg_path}")
    print(f"metadata={json_path}")


if __name__ == "__main__":
    main()
