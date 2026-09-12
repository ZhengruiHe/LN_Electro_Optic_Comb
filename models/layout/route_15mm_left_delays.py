"""15mm上方布局的左侧延时路线试算；有交叉时只报告，不输出GDS。

固定YSJ、15mm电极、750um MUX本体，采用用户同意的100um外接续。
沿已有ng插值计算工程时延，不启动电磁求解或改变源版图。
"""
import argparse
import json
import math
import sys
from pathlib import Path
import numpy as np
from shapely.geometry import LineString, box, Point
from shapely.ops import substring
from shapely.affinity import translate, scale
from euler_routes import append_euler, bend_displacement, radius_for_semicircle_displacement
from route_geometry import check_route_network
from search_fixed_block_placement import geometry, preview_one_placement
from preview_15mm_upper_candidate import polygons

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/"scripts"))
from delay_budget import route_delay


class PathBuilder:
    def __init__(self, points, heading):
        self.points = [list(p) for p in points]
        self.heading = heading
        self.bends = []

    def straight(self, x=None, y=None):
        p = self.points[-1]
        q = [p[0] if x is None else x, p[1] if y is None else y]
        dx, dy = q[0]-p[0], q[1]-p[1]
        forward = dx*math.cos(self.heading)+dy*math.sin(self.heading)
        lateral = -dx*math.sin(self.heading)+dy*math.cos(self.heading)
        if forward < -1e-6 or abs(lateral) > 1e-6:
            raise ValueError(f"直段反向/斜接：{p}->{q}, heading={self.heading}")
        if forward > 1e-8:
            self.points.append(q)

    def bend(self, angle, radius, label):
        self.heading = append_euler(self.points, self.heading, angle, radius, self.bends, label)


def build_left(prefix, index, xr, x_left, y_target, cfg):
    q = bend_displacement(math.pi/2, 80)[0]
    radius_end = cfg["final_radius"]
    y_final_start = (cfg["final_start_y"] if cfg["final_rectangular"] else
                     y_target-bend_displacement(math.pi, radius_end)[1])
    p = PathBuilder(prefix, math.pi)
    p.straight(x=cfg["down_col"]+q)
    p.bend(math.pi/2, 80, "左侧转下")
    p.straight(y=cfg["entry_y"]+q)
    p.bend(math.pi/2, 80, "进入延时直段")
    p.straight(x=xr)
    p.bend(math.pi, 80, "延时折返1")
    if cfg["u_count"] == 3:
        p.straight(x=cfg["fold_left"])
        p.bend(-math.pi, 80, "延时折返2")
        p.straight(x=xr)
        p.bend(math.pi, 80, "延时折返3")
    p.straight(x=cfg["return_col"]+q)
    p.bend(-math.pi/2, 80, "返回竖直段")
    p.straight(y=y_final_start-q)
    p.bend(math.pi/2, 80, "转向MUX左侧")
    p.straight(x=cfg["final_x"])
    if cfg["final_rectangular"]:
        q_end = bend_displacement(math.pi/2, radius_end)[0]
        p.bend(-math.pi/2, radius_end, "MUX前上转1")
        p.straight(y=y_target-q_end)
        p.bend(-math.pi/2, radius_end, "MUX前上转2")
    else:
        p.bend(-math.pi, radius_end, "MUX前上折返")
    if abs(p.points[-1][1]-y_target) > 1e-6:
        raise ValueError("MUX端口高度不匹配")
    p.straight(x=x_left)
    return p.points, p.bends


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output-dir", type=Path, required=True)
    ap.add_argument("--include-io", action="store_true")
    a = ap.parse_args()
    a.output_dir.mkdir(parents=True, exist_ok=False)
    base = ROOT/"results/layout/_历史归档/迭代版本_20260912/合并YSJ_版图确认_20260909_v1"
    study = json.loads((base/"15mm重新布局可行性_v9_四出口完整局部检查/15mm上方模块可行性.json").read_text(encoding="utf-8"))
    budget = json.loads((base/"15mm重新布局可行性_v9_四出口完整局部检查/局部候选与时延保留检查.json").read_text(encoding="utf-8"))
    ng = json.loads((ROOT/"results/optical/欧拉时延输入_0p7um方向群折射率_v1/0p7um_rib_两个晶向群折射率.json").read_text(encoding="utf-8"))
    ny, nz = ng["axes"]["Y"]["group_index"], ng["axes"]["Z"]["group_index"]
    dx, dy = -11620., 1620.
    x_left = 950+dx
    ends = [33.4+dy, -30.15+dy, -33.4+dy]
    d80 = bend_displacement(math.pi, 80)[1]
    configs = [
        dict(down_col=-10845., entry_y=650., return_col=-10200., u_count=1,
             fold_left=None, final_rectangular=True, final_start_y=1050.,
             final_radius=83., final_x=-10728.638053930058),
        dict(down_col=-10865., entry_y=300., return_col=-8700., u_count=3,
             fold_left=-8300., final_rectangular=True, final_start_y=1180.,
             final_radius=80., final_x=-10720.),
        dict(down_col=-10890., entry_y=-400., return_col=-10530., u_count=1,
             fold_left=None, final_rectangular=False, final_start_y=None,
             final_radius=80., final_x=-10671.),
    ]
    routes, metas = {}, {}
    for i, (cfg, xr, target_y) in enumerate(zip(configs, [-8800., -6500., -10100.], ends), 1):
        prefix = [[x+dx, y+dy] for x, y in study["right_escape_local_points_um"][f"loop{i}"]]
        target = budget["delay_targets_retained"][i-1]["new_required_passive_route_delay_ps_estimate"]
        history = []
        for _ in range(4):
            points, bends = build_left(prefix, i, xr, x_left, target_y, cfg)
            d = route_delay(points, ny, nz)
            error = target-d["directional_delay_ps"]
            history.append({"right_fold_x_um": xr, "delay_ps": d["directional_delay_ps"], "error_ps": error})
            xr += error*299.792458/ny/(cfg["u_count"]+1)
        points, bends = build_left(prefix, i, xr, x_left, target_y, cfg)
        routes[f"loop{i}"] = points
        d = route_delay(points, ny, nz)
        metas[f"loop{i}"] = {"target_passive_delay_ps": target, **d,
                              "residual_ps": d["directional_delay_ps"]-target,
                              "right_fold_x_um": xr, "config": cfg, "bends": bends, "solve_history": history}
    if a.include_io:
        q = bend_displacement(math.pi/2, 80)[0]
        incoming = PathBuilder([[-10380., -1700.]], 0.)
        incoming.straight(x=-10280.)  # 预留端面耦合器1.2/0.7um接续拉锥
        incoming.straight(x=-8300.)
        incoming.bend(math.pi, 80., "输入底部折返")
        incoming.straight(x=-8800.+q)
        incoming.bend(-math.pi/2, 80., "输入转北")
        incoming.straight(y=1120.-q)
        incoming.bend(math.pi/2, 80., "输入转西")
        incoming.straight(x=-10729.638053930058)
        incoming.bend(-math.pi/2, 80., "输入MUX前北转")
        incoming.straight(y=dy+30.15-q)
        incoming.bend(-math.pi/2, 80., "输入MUX前东转")
        incoming.straight(x=x_left)
        routes['input'] = incoming.points
        metas['input'] = {'bends': incoming.bends, **route_delay(incoming.points, ny, nz)}
        prefix = [[x+dx, y+dy] for x,y in study['right_escape_local_points_um']['output']]
        outgoing = PathBuilder(prefix, math.pi)
        outgoing.straight(x=-10877.+q)
        outgoing.bend(math.pi/2, 80., "输出转下")
        outgoing.straight(y=100.+q)
        outgoing.bend(math.pi/2, 80., "输出向东绕行")
        outgoing.straight(x=-9500.)
        output_radius = radius_for_semicircle_displacement(600., 80.)
        outgoing.bend(-math.pi, output_radius, "输出底部折返")
        outgoing.straight(x=-10280.)
        outgoing.straight(x=-10380.)
        routes['output'] = outgoing.points
        metas['output'] = {'bends': outgoing.bends, **route_delay(outgoing.points, ny, nz)}
    topology = check_route_network(routes)
    crossing_details = []
    for cross in topology["unexpected_intersections"]:
        rec = dict(cross)
        if "x_um" in cross:
            point = Point(cross["x_um"], cross["y_um"])
            rec["straight_245um_each_route"] = {}
            for name in cross["routes"]:
                line = LineString(routes[name]); distance = line.project(point)
                sub = substring(line, distance-122.5, distance+122.5)
                rec["straight_245um_each_route"][name] = (
                    distance >= 122.5 and distance+122.5 <= line.length
                    and abs(sub.length-math.dist(sub.coords[0], sub.coords[-1])) < 1e-5)
        crossing_details.append(rec)
    old = json.loads((base/"YSJ源版图几何.json").read_text(encoding="utf-8"))
    old_core = geometry(old, (20, 0))
    shape_checks = {}
    for name, points in routes.items():
        core = LineString(points).buffer(.35, cap_style=2, join_style=2)
        window = LineString(points).buffer(8.6, cap_style=2, join_style=2)
        shape_checks[name] = {"YSJ_core_overlap_um2": old_core.intersection(core).area,
                              "window_inside_block": box(-10900,-1900,10900,1900).covers(window),
                              "window_bbox_um": list(window.bounds)}
    result = {"status": "left_routes_trial_no_GDS", "electrode_length_um": 15000,
              "shift_um": [dx,dy], "routes": routes, "delays": metas,
              "topology": topology, "crossing_details": crossing_details,
              "shape_checks": shape_checks, "input_and_final_output_included": a.include_io,
              "full_layout_pass": False}
    with (a.output_dir/"左侧延时路线试算.json").open("x",encoding="utf-8") as f:
        json.dump(result,f,ensure_ascii=False,indent=2)
    print(json.dumps({k:v for k,v in result.items() if k not in ("routes","delays","topology")},ensure_ascii=False,indent=2))
    print(json.dumps({k:{q:r[q] for q in ("right_fold_x_um","directional_delay_ps","residual_ps") if q in r} for k,r in metas.items()},indent=2))


if __name__ == "__main__":
    main()
