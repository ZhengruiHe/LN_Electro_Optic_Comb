"""画实际GDS二维近邻仿真结构，不是传播场/模式结果。"""
import argparse
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--preparation", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError("不覆盖旧结构图")
    report = json.loads((a.preparation / "几何准备核对.json").read_text(encoding="utf-8"))
    plt.rcParams.update({"font.sans-serif": ["Microsoft YaHei", "SimHei", "DejaVu Sans"], "axes.unicode_minus": False})
    fig, axes = plt.subplots(2, 2, figsize=(11, 9), layout="constrained")
    names = {"left_upper": "左上", "left_lower": "左下", "right_upper": "右上", "right_lower": "右下"}
    for ax, entry in zip(axes.flat, report['ready']):
        c = json.loads(Path(entry['config']).read_text(encoding='utf-8'))
        for vertices in c['slab_polygons_um']:
            q = np.array(vertices)
            ax.fill(q[:, 0], q[:, 1], color='#e4edf4', label=None)
        for vertices in c['core_polygons_um']:
            q = np.array(vertices)
            ax.fill(q[:, 0], q[:, 1], color='#1263b3')
            ax.plot(q[:, 0], q[:, 1], color='#1263b3', lw=.7)
        for name, (normal, pos, center, direction) in c['ports'].items():
            span = c['port_spans_um'][name]
            ax.plot([pos, pos], [center-span/2, center+span/2], color='#cc472a', lw=1.8)
            ax.annotate({'near_pair':'近端双波导', 'far_a':'输入 A', 'far_b':'输入 B'}[name],
                        (pos, center), xytext=(-8 if c['pair'].startswith('left') else 8, 8),
                        textcoords='offset points', fontsize=8, ha='right' if c['pair'].startswith('left') else 'left')
        b = c['bounds_um']
        ax.set(xlim=(b[0]-4, b[2]+4), ylim=(b[1]-4, b[3]+4), aspect='equal',
               xlabel='局部 x（µm）', ylabel='局部 y（µm）',
               title=names[c['pair']]+'：'+ ' / '.join(c['route_names']))
        ax.grid(alpha=.15)
    fig.suptitle('实际 GDS 提取的四组二维近邻结构\n蓝色：实际 LN1 芯层；浅蓝：既有等效平台假设；红色：端口', fontsize=12)
    fig.savefig(a.output, dpi=160)
    plt.close(fig)


if __name__ == '__main__':
    main()
