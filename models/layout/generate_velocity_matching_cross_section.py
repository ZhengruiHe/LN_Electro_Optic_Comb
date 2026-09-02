"""绘制候选A的PDK截面与微波减速机理示意图。"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Ellipse, PathPatch, Polygon, Rectangle
from matplotlib.path import Path as MplPath


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "results" / "layout"


def configure_font() -> None:
    candidates = [
        Path(r"C:\Windows\Fonts\msyh.ttc"),
        Path(r"C:\Windows\Fonts\simhei.ttf"),
    ]
    for candidate in candidates:
        if candidate.is_file():
            font_manager.fontManager.addfont(str(candidate))
            plt.rcParams["font.family"] = font_manager.FontProperties(
                fname=str(candidate)
            ).get_name()
            break
    plt.rcParams["axes.unicode_minus"] = False


def add_field_line(ax: plt.Axes, x0: float, x1: float, y: float, height: float) -> None:
    path = MplPath(
        [(x0, y), (x0 + 0.25 * (x1 - x0), y - height),
         (x0 + 0.75 * (x1 - x0), y - height), (x1, y)],
        [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4],
    )
    ax.add_patch(
        PathPatch(path, facecolor="none", edgecolor="#cc3d5a", lw=1.4, alpha=0.78)
    )


def main() -> None:
    configure_font()
    OUTPUT.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(12.8, 7.2), dpi=180)
    ax.set_xlim(-175, 175)
    ax.set_ylim(-6.5, 14.0)
    ax.axis("off")

    # 垂直方向为清晰展示而压缩了硅衬底厚度。
    ax.add_patch(Rectangle((-170, -6), 340, 6, color="#b8bcc5"))
    ax.add_patch(Rectangle((-170, 0), 340, 8.2, color="#bfe2f3"))
    ax.add_patch(Rectangle((-170, 8.2), 340, 0.3, color="#bfe2f3"))

    rail_centers = (-40.5, 40.5)
    ax.add_patch(Rectangle((-170, 8.5), 340, 0.4, color="#d9eef7"))
    ax.add_patch(Rectangle((-170, 8.9), 340, 1.4, color="#d9eef7"))

    # 两级LN刻蚀：200-nm残余薄层 + 200-nm高的梯形脊，两层侧壁均为60°。
    for center in rail_centers:
        ax.add_patch(Polygon([
            (center - 3.7155, 8.9), (center + 3.7155, 8.9),
            (center + 3.6, 9.1), (center - 3.6, 9.1),
        ], closed=True, color="#8f79bd"))
        ax.add_patch(Polygon([
            (center - 0.9155, 9.1), (center + 0.9155, 9.1),
            (center + 0.8, 9.3), (center - 0.8, 9.3),
        ], closed=True, color="#7155a8"))
        ax.add_patch(
            Ellipse(
                (center, 9.12), 4.8, 0.43,
                facecolor="#f04f73", edgecolor="none", alpha=0.34, zorder=5,
            )
        )

    # T形加载最强截面：60-um信号主干、100-um地、5-um内间隙。
    metal_y = 10.3
    metal_h = 1.0
    metal = "#c58b2a"
    ax.add_patch(Rectangle((-30, metal_y), 60, metal_h, color=metal))
    ax.add_patch(Rectangle((-151, metal_y), 100, metal_h, color=metal))
    ax.add_patch(Rectangle((51, metal_y), 100, metal_h, color=metal))

    # 信号和地的T帽/颈，使两个局部调制间隙均为5 um。
    for x, width in [(-38, 8), (30, 8), (-51, 8), (43, 8)]:
        ax.add_patch(Rectangle((x, metal_y), width, metal_h, color=metal))

    add_field_line(ax, -38, -43, 10.3, 1.0)
    add_field_line(ax, 38, 43, 10.3, 1.0)
    add_field_line(ax, -37.2, -43.8, 10.3, 0.65)
    add_field_line(ax, 37.2, 43.8, 10.3, 0.65)

    # 主要结构标注。
    ax.text(0, -3.2, "高阻硅衬底  ρ≈10000 Ω·cm\n（实际约716.8 µm，图中压缩显示）",
            ha="center", va="center", fontsize=12, color="#30343b")
    ax.text(-166, 4.1, "底部SiO2\n8.2 µm", ha="left", va="center", fontsize=11)
    ax.text(-166, 8.35, "SiO₂填充300 nm\n（SiN全部移除）", ha="left", va="center", fontsize=9)
    ax.text(-166, 8.70, "层间SiO2 400 nm", ha="left", va="center", fontsize=10)
    ax.text(-166, 9.18, "X切LN 400 nm", ha="left", va="center", fontsize=10)
    ax.text(-166, 9.75, "顶部SiO2 1 µm", ha="left", va="center", fontsize=10)
    ax.text(0, 10.80, "信号电极 60 µm", ha="center", va="center", fontsize=11)
    ax.text(-101, 10.80, "地 100 µm", ha="center", va="center", fontsize=10)
    ax.text(101, 10.80, "地 100 µm", ha="center", va="center", fontsize=10)

    for center in rail_centers:
        ax.annotate(
            "1.6 µm LN顶宽\n下方SiN全部移除",
            xy=(center, 9.08), xytext=(center, 12.0),
            ha="center", va="bottom", fontsize=10,
            arrowprops={"arrowstyle": "->", "color": "#465267", "lw": 1.0},
        )
    ax.annotate(
        "5 µm调制间隙\nT形电容加载",
        xy=(40.5, 10.35), xytext=(72, 12.7),
        ha="center", va="bottom", fontsize=10,
        arrowprops={"arrowstyle": "->", "color": "#465267", "lw": 1.0},
    )

    ax.text(
        168, 5.5,
        "为什么Si衬底仍能匹配？\n\n"
        "① 8.2 µm厚SiO2把RF场与高ε硅隔开\n"
        "② 高阻硅降低衬底导电损耗\n"
        "③ 有源区完全去除SiN，减少寄生介质加载\n"
        "④ T形电极增加单位长度电容，使微波减速\n"
        "⑤ 60 µm主干与间隙同时回调阻抗和nRF\n\n"
        "目标：nRF ≈ ng ≈ 2.2–2.3\n"
        "等效介电常数约 nRF² ≈ 5",
        ha="right", va="center", fontsize=10.5,
        bbox={"boxstyle": "round,pad=0.55", "fc": "#f2f4f7", "ec": "#adb4c0"},
    )

    ax.text(
        0, 13.65,
        "候选A横截面：利用厚氧化层隔离硅衬底，再用T形电容加载精调微波群速度",
        ha="center", va="top", fontsize=15, fontweight="bold",
    )
    ax.text(
        0, -6.25,
        "注意：所有含有限宽SiN隔离窗口的旧HFSS数值均已作废；当前全SiO₂填充截面需要重新求解。",
        ha="center", va="bottom", fontsize=10, color="#9b253c",
    )

    fig.tight_layout(pad=0.4)
    png = OUTPUT / "候选A速度匹配截面图.png"
    svg = OUTPUT / "候选A速度匹配截面图.svg"
    fig.savefig(png, dpi=220, bbox_inches="tight", facecolor="white")
    fig.savefig(svg, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"png={png}")
    print(f"svg={svg}")


if __name__ == "__main__":
    main()
