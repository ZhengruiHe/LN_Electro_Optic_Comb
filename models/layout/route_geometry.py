"""光路中心线检查：输入/输出引线必须与回路一起检查。

这是名义几何筛查，不是光学性能验证，也不代替掩膜 DRC。
允许交叉必须明确指定两条光路、位置和夹角；不能全局忽略交叉。
"""

from __future__ import annotations

import itertools
import math

from shapely.geometry import LineString, Point


def crossing_angle_deg(line_a: LineString, line_b: LineString, point: Point) -> float:
    """交点两侧各取 0.1 um，估计局部无向锐夹角。"""
    vectors = []
    for line in (line_a, line_b):
        s = line.project(point)
        p0 = line.interpolate(max(0.0, s - 0.1))
        p1 = line.interpolate(min(line.length, s + 0.1))
        vectors.append((p1.x - p0.x, p1.y - p0.y))
    a, b = vectors
    norm = math.hypot(*a) * math.hypot(*b)
    if norm == 0:
        raise ValueError("交点处无法确定切向")
    cosine = abs(a[0] * b[0] + a[1] * b[1]) / norm
    return math.degrees(math.acos(min(1.0, cosine)))


def check_route_network(routes: dict[str, list[list[float]]],
                        allowed_crossings: list[dict] | None = None) -> dict:
    lines = {name: LineString(points) for name, points in routes.items()}
    allowed = allowed_crossings or []
    found_allowed = set()
    pairs = {}
    unexpected = []
    for first, second in itertools.combinations(lines, 2):
        a, b = lines[first], lines[second]
        intersection = a.intersection(b)
        record = {
            "centerline_distance_um": a.distance(b),
            "intersection_type": intersection.geom_type if not intersection.is_empty else None,
            "intersections": [],
        }
        if not intersection.is_empty:
            points = list(intersection.geoms) if intersection.geom_type == "MultiPoint" else [intersection]
            for point in points:
                item = {"wkt": point.wkt, "allowed_reserved_crossing": False}
                if point.geom_type == "Point":
                    item.update(x_um=point.x, y_um=point.y,
                                angle_deg=crossing_angle_deg(a, b, point))
                    for index, spec in enumerate(allowed):
                        if (set(spec["routes"]) == {first, second}
                            and point.distance(Point(spec["center_um"])) <= 0.001
                            and abs(item["angle_deg"] - spec.get("angle_deg", 90)) <= 0.01):
                            item["allowed_reserved_crossing"] = True
                            found_allowed.add(index)
                record["intersections"].append(item)
                if not item["allowed_reserved_crossing"]:
                    unexpected.append({"routes": [first, second], **item})
        pairs[f"{first}__{second}"] = record
    self_crossing = [name for name, line in lines.items() if not line.is_simple]
    missing = [spec for i, spec in enumerate(allowed) if i not in found_allowed]
    return {
        "scope": "输入、输出、三条回路；中心线级；不含器件内部/包络/黑盒真实掩膜",
        "self_crossing_routes": self_crossing,
        "pairwise_checks": pairs,
        "unexpected_intersections": unexpected,
        "missing_reserved_crossings": missing,
        "reserved_crossings": allowed,
        "topology_screen_pass": not (self_crossing or unexpected or missing),
        "optical_function_pass": False,
    }
