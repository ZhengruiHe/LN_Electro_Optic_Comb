"""从KLayout只读导出的多边形绘制合并位置预览，不写GDS。"""
import argparse
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.patches import Rectangle


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--geometry", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--title", default="YSJ 原始版图：仅查看，没有修改结构")
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    data = json.loads(args.geometry.read_text(encoding="utf-8"))
    plt.rcParams.update({"font.sans-serif": ["Microsoft YaHei", "SimHei", "DejaVu Sans"],
                         "axes.unicode_minus": False})
    fig, ax = plt.subplots(figsize=(17, 4.4), layout="constrained")
    colors = {(10, 2): "#d2e8f0", (20, 1): "#91c8b0", (20, 0): "#225e48",
              (42, 0): "#b9882c", (41, 0): "#c7656c", (10, 0): "#705cb4",
              (10, 1): "#b8adde", (20, 2): "#a6b2c3", (21, 2): "#6f7f93",
              (100, 0): "#d1683e"}
    for key in colors:
        layer = next((row for row in data if tuple(row["layer"]) == key), None)
        if layer is None:
            continue
        for record in layer["polygons"]:
            vertices = [[x / 1000, y / 1000] for x, y in record["exterior"]]
            ax.add_collection(PolyCollection([vertices], facecolors=colors[key],
                                            edgecolors=colors[key], linewidths=.25,
                                            alpha=.75 if key != (100, 0) else .15))
            for hole in record["holes"]:
                ax.add_collection(PolyCollection([[[x/1000, y/1000] for x, y in hole]],
                                                facecolors="white", edgecolors=colors[key],
                                                linewidths=.2))
    ax.add_patch(Rectangle((-10.9, -1.9), 21.8, 3.8, fill=False, edgecolor="#394b68",
                           linestyle="--", linewidth=1.2))
    ax.text(-5.5, 0, "左侧大部分为空\n约 10.9 mm 横向空间", ha="center", va="center",
            fontsize=12, color="#576780")
    ax.set(xlim=(-11.2, 11.2), ylim=(-2.05, 2.05), aspect="equal",
           xlabel="x（mm）", ylabel="y（mm）", title=args.title)
    ax.grid(alpha=.13)
    fig.savefig(args.output, dpi=170)
    plt.close(fig)
    print(args.output.resolve())


if __name__ == "__main__":
    main()
