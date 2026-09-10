"""通过 KLink 在 KLayout 中生成四程电光梳调制器的名义 GDS。

本脚本只生成用于功能级联和长度回标的名义版图，不代表通过代工方 DRC。
四次有源调制均沿 +x 方向传播，与射频行波同向；三个回路把光从右端
送回左端。当前 15 mm、10 GHz 基线采用 3T/3.5T/3T，以满足真实折返
几何和 80 um 最小弯曲半径。
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

from klink import KLinkClient
from shapely.geometry import LineString

from route_geometry import check_route_network


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "results" / "layout"
PDK_ROOT = ROOT / "ioptee_sin_tfln_v1.0_beta" / "ioptee_sin_tfln_v1.0_beta"
PDK_LYP = PDK_ROOT / "techfiles" / "klayout" / "ioptee_sin_tfln.lyp"
EDGE_COUPLER_GDS = (
    ROOT
    / "SiN_TFLN_V0.1_BlackBox"
    / "edge_coupler"
    / "edge_coupler_9um_y_1550_ln.gds"
)
MUX_GROUP_DELAY_JSON = (
    ROOT / "results" / "optical" / "模式复用器_70度_EME_群时延摘要.json"
)

TOP_CELL = "EO4P_10G_15MM_3T_3P5T_3T_NOMINAL"
MUX_CELL = "MUX_TE01_A70_EME_NOMINAL"
ELECTRODE_CELL = "T_GSG_15MM_RECTANGULAR_BASELINE"
ACTIVE_CELL = "DUAL_RAIL_ACTIVE_LN"

LAYERS = {
    "TRENCH_SIN": (10, 2),
    "WGCORE_LN1": (20, 0),
    "WGCLAD_LN1": (20, 1),
    "METAL": (42, 0),
    "BB_BND": (100, 0),
    "BB_PARAM": (101, 0),
}


def parse_periods(value: str) -> list[float]:
    periods = [float(item.strip()) for item in value.split(",")]
    if len(periods) != 3 or min(periods) <= 0:
        raise argparse.ArgumentTypeError("回路周期必须是三个正数，例如3,3.5,3")
    if not math.isclose(periods[0] % 1.0, 0.0, abs_tol=1e-9):
        raise argparse.ArgumentTypeError("回路1必须是整数周期")
    if not math.isclose(periods[1] % 1.0, 0.5, abs_tol=1e-9):
        raise argparse.ArgumentTypeError("回路2必须是半整数周期，用于补偿上下电场反号")
    if not math.isclose(periods[2] % 1.0, 0.0, abs_tol=1e-9):
        raise argparse.ArgumentTypeError("回路3必须是整数周期")
    return periods


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--active-length-mm", type=float, default=15.0)
    parser.add_argument("--frequency-ghz", type=float, default=10.0)
    parser.add_argument("--loop-periods", type=parse_periods, default=[3.0, 3.5, 3.0])
    parser.add_argument(
        "--replace-active-layout",
        action="store_true",
        help="清空当前KLayout活动页；运行前应先保存其他版图",
    )
    parser.add_argument(
        "--skip-coupler-blackboxes",
        action="store_true",
        help="不放置PDK端面耦合器黑盒",
    )
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8765)
    parser.add_argument(
        "--allow-unverified-crossings", action="store_true",
        help="仅供复现旧版交叉反例；允许导出未验证交叉，绝不表示光学通过",
    )
    return parser.parse_args()


def polygon_from_center_width(
    x_values: list[float], centers: list[float], widths: list[float]
) -> list[list[float]]:
    lower = [
        [x_value, center - 0.5 * width]
        for x_value, center, width in zip(x_values, centers, widths)
    ]
    upper = [
        [x_value, center + 0.5 * width]
        for x_value, center, width in zip(x_values, centers, widths)
    ]
    return lower + list(reversed(upper))


def taper_polygon(
    x0: float,
    x1: float,
    center0: float,
    center1: float,
    width0: float,
    width1: float,
) -> list[list[float]]:
    return [
        [x0, center0 - 0.5 * width0],
        [x1, center1 - 0.5 * width1],
        [x1, center1 + 0.5 * width1],
        [x0, center0 + 0.5 * width0],
    ]


def mode_mux_paths(samples: int = 121) -> dict[str, list[float]]:
    """复现已完成1550 nm名义EME检查的渐变轨迹。"""
    length_um = 750.0
    main_start_um = 1.40
    main_stop_um = 1.10
    aux_start_um = 0.30
    aux_stop_um = 0.40
    end_gap_um = 2.50
    minimum_gap_um = 1.00
    exponent = 2.0
    ratio = 2.3333333333

    x_values: list[float] = []
    main_widths: list[float] = []
    aux_widths: list[float] = []
    main_centers: list[float] = []
    aux_centers: list[float] = []
    for index in range(samples):
        u = index / (samples - 1)
        numerator = u**exponent
        denominator = numerator + ratio * (1.0 - u) ** exponent
        progress = numerator / denominator if denominator > 0 else 0.0
        if index == samples - 1:
            progress = 1.0
        main_width = main_start_um + (main_stop_um - main_start_um) * progress
        aux_width = aux_start_um + (aux_stop_um - aux_start_um) * progress
        gap = minimum_gap_um + (end_gap_um - minimum_gap_um) * (
            0.5 + 0.5 * math.cos(2.0 * math.pi * u)
        )
        x_values.append(length_um * u)
        main_widths.append(main_width)
        aux_widths.append(aux_width)
        main_centers.append(-0.5 * gap - 0.5 * main_width)
        aux_centers.append(0.5 * gap + 0.5 * aux_width)
    return {
        "x": x_values,
        "main_width": main_widths,
        "aux_width": aux_widths,
        "main_center": main_centers,
        "aux_center": aux_centers,
    }


def polyline_length(points: list[list[float]]) -> float:
    return sum(
        math.hypot(p1[0] - p0[0], p1[1] - p0[1])
        for p0, p1 in zip(points, points[1:])
    )


def append_arc(
    points: list[list[float]],
    center_x: float,
    center_y: float,
    radius: float,
    angle_start: float,
    angle_stop: float,
    segments: int = 64,
) -> None:
    for index in range(1, segments + 1):
        angle = angle_start + (angle_stop - angle_start) * index / segments
        points.append(
            [
                center_x + radius * math.cos(angle),
                center_y + radius * math.sin(angle),
            ]
        )


def semicircle_chord_length(radius: float, segments: int = 64) -> float:
    return 2.0 * radius * segments * math.sin(math.pi / (2.0 * segments))


def return_route(
    *,
    x_right: float,
    y_start: float,
    x_left: float,
    y_end: float,
    target_length: float,
    outward_sign: int,
    start_radius: float,
    middle_radius: float = 80.0,
    launch_extension: float = 5.0,
    end_extension: float = 5.0,
) -> tuple[list[list[float]], dict[str, float]]:
    """生成带两级蛇形的右端到左端回路，并精确回标折线中心线长度。"""
    if outward_sign not in (-1, 1):
        raise ValueError("outward_sign只能取-1或1")
    sign = float(outward_sign)
    lane1_y = y_start + sign * 2.0 * start_radius
    lane2_y = lane1_y + sign * 2.0 * middle_radius
    lane3_y = lane2_y + sign * 2.0 * middle_radius
    end_radius = sign * (lane3_y - y_end) / 2.0
    if min(start_radius, middle_radius, end_radius) < 80.0 - 1e-9:
        raise ValueError("回路中出现小于80 um的弯曲半径")

    arc_length = (
        semicircle_chord_length(start_radius)
        + 2.0 * semicircle_chord_length(middle_radius)
        + semicircle_chord_length(end_radius)
    )
    direct_span = x_right - x_left
    horizontal_backtrack = 0.5 * (
        target_length
        - direct_span
        - 2.0 * launch_extension
        - 2.0 * end_extension
        - arc_length
    )
    if horizontal_backtrack <= 0:
        minimum_length = direct_span + 2.0 * launch_extension + arc_length
        raise ValueError(
            f"目标回路{target_length:.3f} um短于当前折返几何最小值"
            f"{minimum_length:.3f} um"
        )

    x_turn = x_left + max(1300.0, end_radius + 300.0)
    x_second = x_turn + horizontal_backtrack
    if x_second + middle_radius >= x_right + launch_extension - start_radius:
        raise ValueError("蛇形横向回折过长，需要增加版图长度或增加蛇形级数")

    points: list[list[float]] = [[x_right, y_start]]
    if launch_extension > 0:
        points.append([x_right + launch_extension, y_start])
    start_x = x_right + launch_extension

    append_arc(
        points,
        start_x,
        y_start + sign * start_radius,
        start_radius,
        -sign * math.pi / 2.0,
        sign * math.pi / 2.0,
    )
    points.append([x_turn, lane1_y])
    append_arc(
        points,
        x_turn,
        lane1_y + sign * middle_radius,
        middle_radius,
        -sign * math.pi / 2.0,
        -sign * 3.0 * math.pi / 2.0,
    )
    points.append([x_second, lane2_y])
    append_arc(
        points,
        x_second,
        lane2_y + sign * middle_radius,
        middle_radius,
        -sign * math.pi / 2.0,
        sign * math.pi / 2.0,
    )
    final_arc_x = x_left - end_extension
    points.append([final_arc_x, lane3_y])
    append_arc(
        points,
        final_arc_x,
        0.5 * (lane3_y + y_end),
        end_radius,
        sign * math.pi / 2.0,
        sign * 3.0 * math.pi / 2.0,
    )
    if end_extension > 0:
        points.append([x_left, y_end])

    actual_length = polyline_length(points)
    if abs(actual_length - target_length) > 1e-6:
        raise RuntimeError(
            f"回路中心线长度回标失败：目标{target_length:.9f} um，"
            f"实际{actual_length:.9f} um"
        )
    return points, {
        "target_length_um": target_length,
        "actual_centerline_length_um": actual_length,
        "start_radius_um": start_radius,
        "middle_radius_um": middle_radius,
        "end_radius_um": end_radius,
        "horizontal_backtrack_um": horizontal_backtrack,
        "launch_extension_um": launch_extension,
        "end_extension_um": end_extension,
        "outer_lane_y_um": lane3_y,
    }


def call(client: KLinkClient, name: str, params: dict[str, Any]) -> Any:
    return client.call(name, params)


def insert_path_family(
    client: KLinkClient,
    cell: str,
    layer_indices: dict[str, int],
    points: list[list[float]],
    core_width: float,
) -> None:
    for layer_name, width in (
        ("TRENCH_SIN", 17.2),
        ("WGCLAD_LN1", 7.2),
        ("WGCORE_LN1", core_width),
    ):
        buffered = LineString(points).buffer(
            0.5 * width,
            quad_segs=8,
            cap_style=2,
            join_style=1,
        )
        if buffered.geom_type != "Polygon" or not buffered.is_valid:
            raise RuntimeError(f"{layer_name}路径缓冲后不是有效单连通多边形")
        call(
            client,
            "shape.insert_polygon",
            {
                "cell": cell,
                "layer_index": layer_indices[layer_name],
                "points_um": [
                    [float(x_value), float(y_value)]
                    for x_value, y_value in list(buffered.exterior.coords)[:-1]
                ],
            },
        )


def add_mux_cell(
    client: KLinkClient, layer_indices: dict[str, int]
) -> dict[str, list[float]]:
    paths = mode_mux_paths()
    call(
        client,
        "shape.insert_boxes",
        {
            "cell": MUX_CELL,
            "layer_index": layer_indices["TRENCH_SIN"],
            "boxes_um": [[-200.0, -8.6, 950.0, 8.6]],
        },
    )
    call(
        client,
        "shape.insert_boxes",
        {
            "cell": MUX_CELL,
            "layer_index": layer_indices["WGCLAD_LN1"],
            "boxes_um": [[-200.0, -3.6, 950.0, 3.6]],
        },
    )
    core_polygons = [
        polygon_from_center_width(
            paths["x"], paths["main_center"], paths["main_width"]
        ),
        polygon_from_center_width(
            paths["x"], paths["aux_center"], paths["aux_width"]
        ),
        taper_polygon(-200.0, 0.0, -1.95, -1.95, 1.33, 1.40),
        taper_polygon(-100.0, 0.0, 1.40, 1.40, 0.15, 0.30),
        taper_polygon(750.0, 950.0, -1.80, -1.80, 1.10, 0.70),
        taper_polygon(750.0, 950.0, 1.45, 1.45, 0.40, 0.70),
    ]
    for polygon in core_polygons:
        call(
            client,
            "shape.insert_polygon",
            {
                "cell": MUX_CELL,
                "layer_index": layer_indices["WGCORE_LN1"],
                "points_um": polygon,
            },
        )
    for label, position in (
        ("ACTIVE_MAIN_1P33", [-200.0, -1.95]),
        ("EXT_MAIN_0P70", [950.0, -1.80]),
        ("EXT_AUX_0P70", [950.0, 1.45]),
        ("AUX_TAIL_0P15_PENDING_EME", [-100.0, 1.40]),
    ):
        call(
            client,
            "shape.insert_text",
            {
                "cell": MUX_CELL,
                "layer_index": layer_indices["BB_PARAM"],
                "string": label,
                "position_um": position,
                "size_um": 5.0,
            },
        )
    return paths


def add_electrode_cell(
    client: KLinkClient,
    layer_indices: dict[str, int],
    active_length_um: float,
) -> int:
    signal_width = 43.0
    ground_width = 100.0
    trunk_gap = 17.0
    neck_length = 4.0
    neck_x_width = 10.0
    cap_y_width = 2.0
    cap_x_length = 45.0
    period = 50.0
    x0 = 0.0
    x1 = active_length_um
    signal_top = 0.5 * signal_width
    signal_bottom = -0.5 * signal_width
    upper_ground_bottom = signal_top + trunk_gap
    upper_ground_top = upper_ground_bottom + ground_width
    lower_ground_top = signal_bottom - trunk_gap
    lower_ground_bottom = lower_ground_top - ground_width

    boxes: list[list[float]] = [
        [x0, signal_bottom, x1, signal_top],
        [x0, upper_ground_bottom, x1, upper_ground_top],
        [x0, lower_ground_bottom, x1, lower_ground_top],
    ]
    x_value = 0.5 * period
    while x_value < active_length_um:
        x_neck0 = x_value - 0.5 * neck_x_width
        x_neck1 = x_value + 0.5 * neck_x_width
        x_cap0 = x_value - 0.5 * cap_x_length
        x_cap1 = x_value + 0.5 * cap_x_length
        boxes.extend(
            [
                [x_neck0, signal_top, x_neck1, signal_top + neck_length],
                [x_cap0, signal_top + neck_length, x_cap1, signal_top + neck_length + cap_y_width],
                [x_neck0, upper_ground_bottom - neck_length, x_neck1, upper_ground_bottom],
                [x_cap0, upper_ground_bottom - neck_length - cap_y_width, x_cap1, upper_ground_bottom - neck_length],
                [x_neck0, signal_bottom - neck_length, x_neck1, signal_bottom],
                [x_cap0, signal_bottom - neck_length - cap_y_width, x_cap1, signal_bottom - neck_length],
                [x_neck0, lower_ground_top, x_neck1, lower_ground_top + neck_length],
                [x_cap0, lower_ground_top + neck_length, x_cap1, lower_ground_top + neck_length + cap_y_width],
            ]
        )
        x_value += period
    call(
        client,
        "shape.insert_boxes",
        {
            "cell": ELECTRODE_CELL,
            "layer_index": layer_indices["METAL"],
            "boxes_um": boxes,
        },
    )
    return len(boxes)


def main() -> None:
    args = parse_args()
    if args.active_length_mm <= 0 or args.frequency_ghz <= 0:
        raise ValueError("有源长度和射频频率必须为正数")
    active_length_um = args.active_length_mm * 1000.0
    if active_length_um >= 21_800.0:
        raise ValueError("有源长度必须小于PDK有效设计长度21.8 mm")

    client = KLinkClient(host=args.host, port=args.port)
    client.connect()
    initial = call(client, "layout.info", {"verbosity": "full"})
    is_blank = (
        initial.get("file") is None
        and not initial.get("layers")
        and initial.get("top_cells") in ([], ["TOP"])
    )
    if not is_blank and not args.replace_active_layout:
        raise RuntimeError(
            "当前KLayout活动页不是空白。请先保存其他版图，再加"
            "--replace-active-layout运行。"
        )
    call(client, "layout.clear", {})
    blank = call(client, "layout.info", {"verbosity": "normal"})
    if blank.get("top_cells") == ["TOP"]:
        call(
            client,
            "cell.rename",
            {"cell": "TOP", "new_name": TOP_CELL},
        )
    else:
        call(client, "cell.create", {"name": TOP_CELL})
    for cell_name in (MUX_CELL, ELECTRODE_CELL, ACTIVE_CELL):
        call(client, "cell.create", {"name": cell_name})

    layer_indices: dict[str, int] = {}
    for name, (layer, datatype) in LAYERS.items():
        result = call(
            client,
            "layer.ensure",
            {"layer": layer, "datatype": datatype, "name": name},
        )
        layer_indices[name] = int(result["layer_index"])

    mux_paths = add_mux_cell(client, layer_indices)
    electrode_box_count = add_electrode_cell(
        client, layer_indices, active_length_um
    )

    active_x0 = 2000.0
    active_x1 = active_x0 + active_length_um
    upper_rail_y = 30.0
    lower_rail_y = -30.0
    for rail_y in (upper_rail_y, lower_rail_y):
        insert_path_family(
            client,
            ACTIVE_CELL,
            layer_indices,
            [[active_x0, rail_y], [active_x1, rail_y]],
            1.33,
        )

    call(
        client,
        "instance.insert_many",
        {
            "parent": TOP_CELL,
            "items": [
                {"child": ELECTRODE_CELL, "position_um": [active_x0, 0.0]},
                {"child": ACTIVE_CELL, "position_um": [0.0, 0.0]},
                {
                    "child": MUX_CELL,
                    "position_um": [active_x0 - 200.0, upper_rail_y + 1.95],
                    "rotation": 180,
                    "mirror": True,
                    "klink_id": "MUX_LT",
                },
                {
                    "child": MUX_CELL,
                    "position_um": [active_x1 + 200.0, upper_rail_y - 1.95],
                    "mirror": True,
                    "klink_id": "MUX_RT",
                },
                {
                    "child": MUX_CELL,
                    "position_um": [active_x0 - 200.0, lower_rail_y - 1.95],
                    "rotation": 180,
                    "klink_id": "MUX_LB",
                },
                {
                    "child": MUX_CELL,
                    "position_um": [active_x1 + 200.0, lower_rail_y + 1.95],
                    "klink_id": "MUX_RB",
                },
            ],
        },
    )

    # MUX局部坐标的外部0.7 um端口为x=950，主/辅中心为-1.80/+1.45。
    # 这一镜像只保证三条回路彼此不交叉；输入/输出仍必须参与全路由检查。
    x_left = active_x0 - 200.0 - 950.0
    x_right = active_x1 + 200.0 + 950.0
    left_upper_main_y = upper_rail_y + 0.15
    left_upper_aux_y = upper_rail_y + 3.40
    right_upper_main_y = upper_rail_y - 0.15
    right_upper_aux_y = upper_rail_y - 3.40
    left_lower_main_y = lower_rail_y - 0.15
    left_lower_aux_y = lower_rail_y - 3.40
    right_lower_main_y = lower_rail_y + 0.15
    right_lower_aux_y = lower_rail_y + 3.40

    speed_of_light_um_per_s = 299_792_458.0 * 1e6
    ng_te0 = 2.22129949419841
    ng_te1 = 2.22060778246259
    ng_passive = 2.20234880605680
    ng_mux_estimate = 2.20
    ng_narrow_port = 2.12461650319239
    active_group_indices = [ng_te0, ng_te1, ng_te0]
    mux_body_length_per_loop_um = 2.0 * 750.0
    active_tapers_per_loop_um = 2.0 * 200.0
    external_tapers_per_loop_um = 2.0 * 200.0
    default_mux_delay_s = (
        mux_body_length_per_loop_um * ng_mux_estimate / speed_of_light_um_per_s
    )
    mux_delays_s = [default_mux_delay_s] * 3
    mux_delay_source = "以750um本体群折射率2.20估算"
    if MUX_GROUP_DELAY_JSON.is_file():
        mux_delay_data = json.loads(MUX_GROUP_DELAY_JSON.read_text(encoding="utf-8"))
        loop_delay_values = list(mux_delay_data["loop_mux_delay_ps"].values())
        if len(loop_delay_values) == 3 and min(loop_delay_values) > 0:
            mux_delays_s = [float(value) * 1e-12 for value in loop_delay_values]
            mux_delay_source = (
                "MODE 2023 R2复传输相位斜率："
                + str(MUX_GROUP_DELAY_JSON.resolve())
            )
    external_taper_ng = 0.5 * (ng_narrow_port + ng_passive)

    route_targets: list[float] = []
    delay_rows: list[dict[str, float]] = []
    for loop_index, (periods, active_ng) in enumerate(
        zip(args.loop_periods, active_group_indices)
    ):
        target_delay_s = periods / (args.frequency_ghz * 1e9)
        active_delay_s = (
            active_length_um * active_ng / speed_of_light_um_per_s
        )
        active_taper_ng = 0.5 * (active_ng + ng_mux_estimate)
        active_taper_delay_s = (
            active_tapers_per_loop_um
            * active_taper_ng
            / speed_of_light_um_per_s
        )
        external_taper_delay_s = (
            external_tapers_per_loop_um
            * external_taper_ng
            / speed_of_light_um_per_s
        )
        passive_route_length_um = (
            speed_of_light_um_per_s
            * (
                target_delay_s
                - active_delay_s
                - mux_delays_s[loop_index]
                - active_taper_delay_s
                - external_taper_delay_s
            )
            / ng_passive
        )
        route_targets.append(passive_route_length_um)
        delay_rows.append(
            {
                "periods": periods,
                "target_delay_ps": target_delay_s * 1e12,
                "active_delay_ps": active_delay_s * 1e12,
                "two_mux_bodies_delay_ps": mux_delays_s[loop_index] * 1e12,
                "two_mux_bodies_delay_source": mux_delay_source,
                "two_active_tapers_delay_ps": active_taper_delay_s * 1e12,
                "two_external_tapers_delay_ps": external_taper_delay_s * 1e12,
                "passive_route_target_um": passive_route_length_um,
            }
        )

    loop1, check1 = return_route(
        x_right=x_right,
        y_start=right_upper_main_y,
        x_left=x_left,
        y_end=left_upper_aux_y,
        target_length=route_targets[0],
        outward_sign=1,
        start_radius=80.0,
    )
    loop2, check2 = return_route(
        x_right=x_right,
        y_start=right_upper_aux_y,
        x_left=x_left,
        y_end=left_lower_main_y,
        target_length=route_targets[1],
        outward_sign=1,
        start_radius=250.0,
        launch_extension=400.0,
    )
    loop3, check3 = return_route(
        x_right=x_right,
        y_start=right_lower_main_y,
        x_left=x_left,
        y_end=left_lower_aux_y,
        target_length=route_targets[2],
        outward_sign=-1,
        start_radius=80.0,
    )
    for loop in (loop1, loop2, loop3):
        insert_path_family(
            client, TOP_CELL, layer_indices, loop, core_width=0.70
        )

    input_port_x = 300.0
    output_port_x = x_right + 550.0
    input_route = [
        [input_port_x, left_upper_main_y],
        [x_left, left_upper_main_y],
    ]
    output_route = [
        [x_right, right_lower_aux_y],
        [output_port_x, right_lower_aux_y],
    ]
    all_route_checks = check_route_network({
        "loop1": loop1, "loop2": loop2, "loop3": loop3,
        "input": input_route, "output": output_route,
    })
    if not all_route_checks["topology_screen_pass"] and not args.allow_unverified_crossings:
        client.close()
        raise RuntimeError(
            "检测到未经设计的输入/输出交叉，已拒绝写入GDS，磁盘历史结果保持不变。"
            "先运行check_crossing_candidate.py检查正交交叉候选；"
            "旧版反例仅可显式使用--allow-unverified-crossings导出。"
        )
    insert_path_family(
        client, TOP_CELL, layer_indices, input_route, core_width=0.70
    )
    insert_path_family(
        client, TOP_CELL, layer_indices, output_route, core_width=0.70
    )

    coupler_status = "未放置"
    if not args.skip_coupler_blackboxes and EDGE_COUPLER_GDS.exists():
        import_result = call(
            client,
            "layout.import_file",
            {
                "path": str(EDGE_COUPLER_GDS.resolve()),
                "on_conflict": "rename",
            },
        )
        imported_top = import_result["file_info"]["top_cells"][0]
        call(
            client,
            "instance.insert_many",
            {
                "parent": TOP_CELL,
                "items": [
                    {
                        "child": imported_top,
                        "position_um": [input_port_x, left_upper_main_y],
                        "klink_id": "PDK_EC_INPUT_BB",
                    },
                    {
                        "child": imported_top,
                        "position_um": [output_port_x, right_lower_aux_y],
                        "rotation": 180,
                        "klink_id": "PDK_EC_OUTPUT_BB",
                    },
                ],
            },
        )
        coupler_status = "已放置PDK的Y方向LN端面耦合器黑盒；黑盒不含真实掩膜"

    annotations = [
        ("RF_IN_REFERENCE_PLANE_NO_PAD", [active_x0, -170.0]),
        ("RF_OUT_REFERENCE_PLANE_NO_TERMINATION", [active_x1 - 800.0, -170.0]),
        ("PASS1_TE0_L2R", [0.45 * active_x0 + 0.55 * active_x1, 42.0]),
        ("PASS2_TE1_L2R", [0.45 * active_x0 + 0.55 * active_x1, 52.0]),
        ("PASS3_TE0_L2R", [0.45 * active_x0 + 0.55 * active_x1, -42.0]),
        ("PASS4_TE1_L2R", [0.45 * active_x0 + 0.55 * active_x1, -52.0]),
        ("LOOP1_3T", [x_left + 500.0, check1["outer_lane_y_um"] + 15.0]),
        ("LOOP2_3P5T", [x_left + 500.0, check2["outer_lane_y_um"] + 15.0]),
        ("LOOP3_3T", [x_left + 500.0, check3["outer_lane_y_um"] - 15.0]),
        ("NOMINAL_ONLY_NOT_DRC_CLOSED", [active_x0, check2["outer_lane_y_um"] + 70.0]),
    ]
    for label, position in annotations:
        call(
            client,
            "shape.insert_text",
            {
                "cell": TOP_CELL,
                "layer_index": layer_indices["BB_PARAM"],
                "string": label,
                "position_um": position,
                "size_um": 18.0,
            },
        )

    # 为尚未设计的GSG焊盘和50 ohm终端保留非工艺层窗口。
    call(
        client,
        "shape.insert_boxes",
        {
            "cell": TOP_CELL,
            "layer_index": layer_indices["BB_BND"],
            "boxes_um": [
                [active_x0 - 160.0, -180.0, active_x0 + 160.0, 180.0],
                [active_x1 - 160.0, -180.0, active_x1 + 160.0, 180.0],
            ],
        },
    )

    OUTPUT.mkdir(parents=True, exist_ok=True)
    stem = (
        f"四程电光梳_{args.frequency_ghz:g}GHz_{args.active_length_mm:g}mm_"
        f"{args.loop_periods[0]:g}T_{args.loop_periods[1]:g}T_"
        f"{args.loop_periods[2]:g}T_名义"
    )
    if not all_route_checks["topology_screen_pass"]:
        stem += "_交叉未验证反例"
    working_gds = (OUTPUT / f"{stem}_工作版.gds").resolve()
    clean_gds = (OUTPUT / f"{stem}_核心工艺层_未DRC.gds").resolve()
    screenshot = (OUTPUT / f"{stem}_KLayout预览.png").resolve()
    electrode_screenshot = (OUTPUT / f"{stem}_电极周期放大.png").resolve()
    mux_screenshot = (OUTPUT / f"{stem}_左端复用器放大.png").resolve()
    report_path = (OUTPUT / f"{stem}_检查.json").resolve()

    call(client, "view.show_cell", {"cell": TOP_CELL, "zoom_fit": True})
    if PDK_LYP.exists():
        call(client, "layer.load_lyp", {"path": str(PDK_LYP.resolve())})
    display_styles = {
        "10/2": ("#BFE8FF", "#75BFE8", 2, 1),
        "20/1": ("#98E4D1", "#36A78A", 3, 1),
        "20/0": ("#E62958", "#A80D36", 0, 2),
        "42/0": ("#E0A52B", "#8A5C00", 0, 1),
        "100/0": ("#FF9F43", "#C96A00", 4, 2),
        "101/0": ("#2CBF6E", "#128A45", 0, 1),
    }
    for layer, (fill, frame, dither, width) in display_styles.items():
        call(
            client,
            "layer.set_style",
            {
                "layer": layer,
                "fill_color": fill,
                "frame_color": frame,
                "dither_pattern": dither,
                "line_width": width,
            },
        )
    call(client, "view.hier_levels", {"min": 0, "max": 10})
    save_result = call(
        client, "layout.save_file", {"path": str(working_gds)}
    )
    clean_result = call(
        client,
        "layout.export_clean",
        {
            "path": str(clean_gds),
            "allowlist_layers": ["10/2", "20/0", "20/1", "42/0"],
            "cells": [TOP_CELL],
        },
    )
    # 端口名和检查文字只用于版图调试。截图中隐藏该注释层，避免遮挡
    # 波导、模式复用器及电极的真实几何；GDS文件本身仍保留注释。
    call(client, "layer.set_visible", {"layers": ["101/0"], "visible": False})
    call(client, "view.zoom_fit", {})
    shot_result = call(
        client,
        "view.screenshot",
        {
            "mode": "path",
            "path": str(screenshot),
            "width_px": 2200,
            "height_px": 900,
        },
    )
    electrode_shot_result = call(
        client,
        "view.screenshot",
        {
            "mode": "path",
            "path": str(electrode_screenshot),
            "bbox_um": [
                active_x0 + 2900.0,
                -150.0,
                active_x0 + 3250.0,
                150.0,
            ],
            "width_px": 1400,
            "height_px": 1200,
        },
    )
    mux_shot_result = call(
        client,
        "view.screenshot",
        {
            "mode": "path",
            "path": str(mux_screenshot),
            "bbox_um": [x_left - 150.0, -80.0, active_x0 + 100.0, 100.0],
            "width_px": 2400,
            "height_px": 500,
        },
    )

    file_info = call(
        client,
        "layout.file_info",
        {"path": str(working_gds), "detail": "counts"},
    )
    clean_info = call(
        client,
        "layout.file_info",
        {"path": str(clean_gds), "detail": "counts"},
    )
    core_metal_overlap = call(
        client,
        "geometry.boolean",
        {
            "a": {"cell": TOP_CELL, "layer": "20/0"},
            "b": {"cell": TOP_CELL, "layer": "42/0"},
            "op": "and",
        },
    )
    loop_checks = [check1, check2, check3]
    loop_lines = [LineString(loop) for loop in (loop1, loop2, loop3)]
    pair_names = [(0, 1, "loop1_loop2"), (0, 2, "loop1_loop3"), (1, 2, "loop2_loop3")]
    route_geometry_checks = {
        "all_centerlines_simple": all(line.is_simple for line in loop_lines),
        "pairwise_intersections": {
            name: loop_lines[first].intersects(loop_lines[second])
            for first, second, name in pair_names
        },
        "pairwise_centerline_distance_um": {
            name: loop_lines[first].distance(loop_lines[second])
            for first, second, name in pair_names
        },
    }
    if not route_geometry_checks["all_centerlines_simple"] or any(
        route_geometry_checks["pairwise_intersections"].values()
    ):
        raise RuntimeError("回路中心线存在自交或回路间相交，拒绝输出通过状态")
    if core_metal_overlap["area_um2"] != 0.0:
        raise RuntimeError("LN芯层与金属层发生几何重叠，拒绝输出通过状态")
    phase_tolerance_deg = 5.0
    delay_tolerance_s = (
        phase_tolerance_deg
        / 360.0
        / (args.frequency_ghz * 1e9)
    )
    route_tolerance_um = (
        speed_of_light_um_per_s * delay_tolerance_s / ng_passive
    )
    all_points = loop1 + loop2 + loop3 + input_route + output_route
    bbox_um = [
        min(point[0] for point in all_points) - 520.0,
        min(point[1] for point in all_points) - 20.0,
        max(point[0] for point in all_points) + 520.0,
        max(point[1] for point in all_points) + 20.0,
    ]
    report = {
        "status": ("存在未验证交叉：反例，不能用于功能通过结论"
                   if not all_route_checks["topology_screen_pass"]
                   else "名义版图已生成；未完成完整级联、终端、焊盘、容差和官方DRC"),
        "working_gds": str(working_gds),
        "clean_process_gds": str(clean_gds),
        "screenshot": str(screenshot),
        "electrode_screenshot": str(electrode_screenshot),
        "mux_screenshot": str(mux_screenshot),
        "top_cell": TOP_CELL,
        "pdk_layer_map": {
            name: f"{layer}/{datatype}"
            for name, (layer, datatype) in LAYERS.items()
        },
        "pdk_local_cross_section": {
            "ln_total_nm": 400,
            "ln1_rib_etch_nm": 200,
            "ln_residual_nm": 200,
            "sidewall_angle_deg_from_horizontal": 70,
            "sin_under_ln": "TRENCH_SIN局部去除，原区域按SiO2填充的名义模型",
        },
        "optical_sequence": [
            {"pass": 1, "rail": "upper", "mode": "TE0", "direction": "+x"},
            {"pass": 2, "rail": "upper", "mode": "TE1", "direction": "+x"},
            {"pass": 3, "rail": "lower", "mode": "TE0", "direction": "+x"},
            {"pass": 4, "rail": "lower", "mode": "TE1", "direction": "+x"},
        ],
        "loop_periods": args.loop_periods,
        "delay_budget": delay_rows,
        "loop_centerline_checks": loop_checks,
        "route_geometry_checks": route_geometry_checks,
        "all_optical_routes_geometry_checks": all_route_checks,
        "core_metal_overlap_check": core_metal_overlap,
        "phase_tolerance_at_target": {
            "phase_deg": phase_tolerance_deg,
            "time_ps": delay_tolerance_s * 1e12,
            "equivalent_passive_route_length_um": route_tolerance_um,
        },
        "nominal_bbox_um_estimate": bbox_um,
        "nominal_footprint_mm_estimate": [
            (bbox_um[2] - bbox_um[0]) / 1000.0,
            (bbox_um[3] - bbox_um[1]) / 1000.0,
        ],
        "mode_mux": {
            "body_length_um": 750.0,
            "eme_checked_body_only": True,
            "group_delay_source": mux_delay_source,
            "main_width_um": [1.40, 1.10],
            "aux_width_um": [0.30, 0.40],
            "gap_um": [2.50, 1.00, 2.50],
            "active_taper_um": 200.0,
            "external_taper_um": 200.0,
            "aux_tail_tip_um": 0.15,
            "aux_tail_status": "待EME/FDTD功能验证",
            "instances": 4,
            "traversals_in_three_loops": 6,
            "mirrored_instances": ["MUX_RT", "MUX_LB"],
            "mirror_note": "仅关于传播轴镜像，不旋转有源传播方向；被动等价性待级联复核",
        },
        "electrode": {
            "active_length_um": active_length_um,
            "metal_layer": "42/0",
            "signal_width_um": 43.0,
            "ground_width_um": 100.0,
            "trunk_gap_um": 17.0,
            "inner_gap_um": 5.0,
            "period_um": 50.0,
            "duty_percent": 90.0,
            "inserted_box_count": electrode_box_count,
            "rf_pads": "仅保留窗口，未画正式探针焊盘",
            "termination": "未画50 ohm终端",
        },
        "edge_coupler": coupler_status,
        "klink_outputs": {
            "save": save_result,
            "clean_export": clean_result,
            "screenshot": shot_result,
            "electrode_screenshot": electrode_shot_result,
            "mux_screenshot": mux_shot_result,
        },
        "file_checks": {
            "working": file_info,
            "clean": clean_info,
        },
        "limitations": [
            "模式复用器750 um本体继承名义EME轨迹，但新增端口锥形和辅助尾端尚未级联验证",
            "回路弯曲当前为圆弧离散中心线；Euler弯曲损耗和模式串扰尚未仿真",
            "复用器本体已回填快速扫频复S相位斜率；固定LN折射率缺少材料色散，不可作为最终群时延签核",
            "未加入正式GSG焊盘、射频过渡、50 ohm终端、切割禁区和官方DRC",
            "本轮只做名义功能版图，不执行工艺容差扫描",
        ],
    }
    report_path.write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(f"working_gds={working_gds}")
    print(f"clean_gds={clean_gds}")
    print(f"screenshot={screenshot}")
    print(f"report={report_path}")
    print(
        "loop_lengths_um="
        + json.dumps(
            [item["actual_centerline_length_um"] for item in loop_checks]
        )
    )


if __name__ == "__main__":
    main()
