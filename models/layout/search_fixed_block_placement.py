"""固定YSJ原图，在21.8×3.8mm内搜索四程器件的整体平移；不写GDS。

只以20/0对20/0、42/0对42/0判重。栅格仅产生候选，再用原轮廓
做平面布尔检查。没有找到栅格候选不等于证明所有连续平移都不可行。
"""
import argparse
import json
import math
from pathlib import Path
import numpy as np
from PIL import Image, ImageDraw
from scipy.signal import fftconvolve
from shapely import make_valid
from shapely.affinity import translate
from shapely.geometry import Polygon, box
from shapely.ops import unary_union
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.patches import Rectangle


BLOCK = [-10900., -1900., 10900., 1900.]
KEYS = [(20, 0), (42, 0)]


def records(data, key):
    return next((layer["polygons"] for layer in data if tuple(layer["layer"]) == key), [])


def geometry(data, key):
    # 自接触轮廓只在分析副本中处理；所有源GDS均保持不变。
    return unary_union([make_valid(Polygon(row["exterior"], row["holes"]))
                        for row in records(data, key)])


def raster(data, key, bounds, step):
    width = math.ceil((bounds[2]-bounds[0])/step)+1
    height = math.ceil((bounds[3]-bounds[1])/step)+1
    im = Image.new("L", (width, height), 0)
    draw = ImageDraw.Draw(im)
    def coords(points):
        return [((x-bounds[0])/step, (y-bounds[1])/step) for x, y in points]
    for row in records(data, key):
        draw.polygon(coords(row["exterior"]), fill=1)
        for hole in row["holes"]:
            draw.polygon(coords(hole), fill=0)
    return np.asarray(im, dtype=np.float32)


def draw_layer(ax, data, key, dx, dy, color):
    lines = []
    for row in records(data, key):
        for ring in [row["exterior"], *row["holes"]]:
            points = [[(x+dx)/1000., (y+dy)/1000.] for x, y in ring]
            lines.append(points+[points[0]])
    ax.add_collection(LineCollection(lines, colors=color, linewidths=.42, alpha=.9))


def preview_one_placement(old, new, dx, dy, output, title):
    """仅生成指定平移的预览图，不搜索、不生成GDS。"""
    output = Path(output)
    if output.exists():
        raise FileExistsError(output)
    plt.rcParams.update({"font.sans-serif": ["Microsoft YaHei", "SimHei", "DejaVu Sans"],
                         "axes.unicode_minus": False})
    fig, ax = plt.subplots(figsize=(17, 4.6), layout="constrained")
    intersections = {}
    for key, c0, c1 in [((20, 0), "#357955", "#d3296a"),
                        ((42, 0), "#aa792a", "#276dc0")]:
        draw_layer(ax, old, key, 0, 0, c0)
        draw_layer(ax, new, key, dx, dy, c1)
        inter = geometry(old, key).intersection(translate(geometry(new, key), dx, dy))
        intersections[str(key)] = {"area_um2": inter.area,
                                   "bbox_um": list(inter.bounds) if not inter.is_empty else None}
        if not inter.is_empty:
            x0, y0, x1, y1 = inter.bounds
            ax.add_patch(Rectangle((x0/1000, y0/1000), (x1-x0)/1000, (y1-y0)/1000,
                                   fill=False, edgecolor="red", linewidth=1.2, linestyle=":"))
    ax.add_patch(Rectangle((-10.9, -1.9), 21.8, 3.8, fill=False,
                           linestyle="--", edgecolor="#333"))
    ax.set(xlim=(-11.2, 11.2), ylim=(-2.05, 2.05), aspect="equal",
           xlabel="x（mm）", ylabel="y（mm）", title=title)
    ax.grid(alpha=.12)
    fig.savefig(output, dpi=180)
    plt.close(fig)
    return {"dx_um": dx, "dy_um": dy, "intersections": intersections,
            "note": "红色虚框仅表示冲突总体范围，不是完整填充禁布区；未写GDS"}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ysj-geometry", type=Path, required=True)
    ap.add_argument("--ours-geometry", type=Path, required=True)
    ap.add_argument("--ours-inventory", type=Path, required=True)
    ap.add_argument("--output-dir", type=Path, required=True)
    ap.add_argument("--step-um", type=float, default=10.)
    a = ap.parse_args()
    if a.step_um <= 0:
        raise ValueError("步长必须为正")
    a.output_dir.mkdir(parents=True, exist_ok=False)
    old = json.loads(a.ysj_geometry.read_text(encoding="utf-8"))
    new = json.loads(a.ours_geometry.read_text(encoding="utf-8"))
    inventory = json.loads(a.ours_inventory.read_text(encoding="utf-8"))
    included = [row["bbox_um"] for row in inventory["tops"][0]["layers"]
                if row["layer"] != 101]
    bounds = [min(b[0] for b in included), min(b[1] for b in included),
              max(b[2] for b in included), max(b[3] for b in included)]
    anchor = [math.floor(bounds[0]/a.step_um)*a.step_um,
              math.floor(bounds[1]/a.step_um)*a.step_um,
              math.ceil(bounds[2]/a.step_um)*a.step_um,
              math.ceil(bounds[3]/a.step_um)*a.step_um]
    shifts = [BLOCK[0]-bounds[0], BLOCK[1]-bounds[1],
              BLOCK[2]-bounds[2], BLOCK[3]-bounds[3]]
    scores = {}
    for key in KEYS:
        obstacle = raster(old, key, BLOCK, a.step_um)
        moving = raster(new, key, anchor, a.step_um)
        scores[key] = np.rint(fftconvolve(obstacle, moving[::-1, ::-1], mode="valid")).astype(np.int32)
    metal_zero = scores[(42, 0)] == 0
    both_zero = metal_zero & (scores[(20, 0)] == 0)
    rank = scores[(42, 0)].astype(float)*1e6 + scores[(20, 0)]
    flat = np.argsort(rank, axis=None, kind="stable")[:120]
    old_geom = {key: geometry(old, key) for key in KEYS}
    new_geom = {key: geometry(new, key) for key in KEYS}
    inspected = []
    for idx in flat:
        row, col = np.unravel_index(idx, rank.shape)
        dx = BLOCK[0]+col*a.step_um-anchor[0]
        dy = BLOCK[1]+row*a.step_um-anchor[1]
        checks = {}
        for key in KEYS:
            moved = translate(new_geom[key], dx, dy)
            inter = old_geom[key].intersection(moved)
            checks[str(key)] = {"intersection_area_um2": inter.area,
                                "intersects": old_geom[key].intersects(moved),
                                "distance_um": old_geom[key].distance(moved),
                                "intersection_bounds_um": list(inter.bounds) if not inter.is_empty else None}
        inspected.append({"dx_um": dx, "dy_um": dy,
                          "raster_overlap_pixels": {str(key): int(scores[key][row,col]) for key in KEYS},
                          "layers": checks,
                          "same_layer_clear": not any(c["intersects"] for c in checks.values())})
        if inspected[-1]["same_layer_clear"]:
            break
    result = {
        "status": "analysis_only_no_GDS_written", "block_um": BLOCK,
        "moving_effective_bbox_um": bounds, "translation_range_um": shifts,
        "rotation_deg": 0, "scale": 1, "raster_step_um": a.step_um,
        "raster_positions": int(rank.size), "metal_clear_raster_positions": int(metal_zero.sum()),
        "both_clear_raster_positions": int(both_zero.sum()),
        "best_raster_core_overlap_with_metal_clear_pixels": int(scores[(20,0)][metal_zero].min()) if metal_zero.any() else None,
        "inspected_candidates": inspected,
        "found_exact_same_layer_clear": any(r["same_layer_clear"] for r in inspected),
        "scope": "只检20/0核心与42/0金属同层冲突；不以跨层投影重叠判失败；窗口与黑盒后续分别记录",
        "limitation": "粗栅格产生候选，未覆盖全部连续平移；未经KLayout原生复核不输出合并GDS",
    }
    with (a.output_dir/"同层放置检查.json").open("x", encoding="utf-8") as stream:
        json.dump(result, stream, ensure_ascii=False, indent=2)
    print(json.dumps({k:v for k,v in result.items() if k != "inspected_candidates"}, ensure_ascii=False, indent=2), flush=True)
    print(json.dumps(inspected[:3], ensure_ascii=False, indent=2), flush=True)
    best = next((r for r in inspected if r["same_layer_clear"]), inspected[0])
    plt.rcParams.update({"font.sans-serif":["Microsoft YaHei","SimHei","DejaVu Sans"],"axes.unicode_minus":False})
    fig, ax = plt.subplots(figsize=(17,4.6), layout="constrained")
    for key, c0, c1 in [((20,0),"#357955","#d3296a"),((42,0),"#aa792a","#276dc0")]:
        draw_layer(ax,old,key,0,0,c0)
        draw_layer(ax,new,key,best["dx_um"],best["dy_um"],c1)
    ax.add_patch(Rectangle((-10.9,-1.9),21.8,3.8,fill=False,linestyle="--",edgecolor="#333"))
    for key in KEYS:
        inter=old_geom[key].intersection(translate(new_geom[key],best["dx_um"],best["dy_um"]))
        parts=list(inter.geoms) if hasattr(inter,"geoms") else [inter]
        for part in parts:
            if part.is_empty: continue
            point=part.representative_point()
            ax.plot(point.x/1000,point.y/1000,"o",mfc="none",mec="red",ms=6)
    ax.set(xlim=(-11.2,11.2),ylim=(-2.05,2.05),aspect="equal",xlabel="x（mm）",ylabel="y（mm）",
           title="固定原YSJ，仅平移四程：绿/棕为原图，粉/蓝为四程；红圈为同层交叠")
    ax.grid(alpha=.12)
    fig.savefig(a.output_dir/"候选放置与同层冲突.png",dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()
