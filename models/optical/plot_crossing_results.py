"""画交叉器场分布和端口模式；仅展示真实FDTD输出，不绘制假想光场。"""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def main():
    p = argparse.ArgumentParser()
    p.add_argument("directory", type=Path)
    a = p.parse_args()
    data = json.loads((a.directory / "结果.json").read_text(encoding="utf-8"))
    plt.rcParams.update({"font.sans-serif": ["Microsoft YaHei", "SimHei", "DejaVu Sans"],
                         "axes.unicode_minus": False})
    field = np.load(a.directory / "场分布.npz")
    actual_field_nm = float(field["lambda"].ravel()[0])*1e9
    e2 = np.sum(np.abs(field["E"].squeeze())**2, axis=-1)
    e2 /= e2.max()
    fig, ax = plt.subplots(figsize=(8, 7), layout="constrained")
    im = ax.pcolormesh(field["x"].ravel()*1e6, field["y"].ravel()*1e6,
                       np.maximum(10*np.log10(e2.T+1e-30), -40),
                       cmap="inferno", vmin=-40, vmax=0, shading="auto")
    pts = np.asarray(data["geometry"]["core_top_polygon_um"])
    ax.plot(pts[:, 0], pts[:, 1], color="#70d6ff", lw=.65, label="LN脊顶部边界")
    ax.set(xlim=(field["x"].min()*1e6, field["x"].max()*1e6),
           ylim=(field["y"].min()*1e6, field["y"].max()*1e6),
           aspect="equal", xlabel="x / 晶体 Y (µm)", ylabel="y / 晶体 Z (µm)",
           title=f"自定义交叉器 · {data['source']} 入射 · 实际场采样 {actual_field_nm:.3f} nm\n"
                 "真实三维 FDTD 场切片 z = 1.0 µm")
    ax.legend(loc="upper left")
    fig.colorbar(im, ax=ax, label="归一化 |E|² (dB；非光功率积分)")
    fig.savefig(a.directory / "交叉器_真实场分布.png", dpi=180)
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), layout="constrained")
    metrics = {}
    for ax, port in zip(axes, ("west", "north")):
        path = a.directory / f"端口模式_{port}.npz"
        if not path.exists():
            ax.set_visible(False)
            continue
        m = np.load(path)
        e = m["E1"].squeeze()
        dens = np.sum(abs(e)**2, axis=-1)
        transverse = "y" if port == "west" else "x"
        in_plane_component = 1 if port == "west" else 0
        t, z = m[transverse].ravel()*1e6, m["z"].ravel()*1e6
        # 用非均匀网格积分，检查入射模式是面内TE，而不是仅相信mode编号。
        integral = lambda v: float(np.trapezoid(np.trapezoid(v, z, axis=1), t))
        te_fraction = integral(abs(e[..., in_plane_component])**2) / integral(dens)
        metrics[port] = {"in_plane_electric_fraction": te_fraction}
        ax.pcolormesh(t, z, dens.T/dens.max(), cmap="inferno", shading="auto")
        ax.set(xlim=(-2, 2), ylim=(0.3, 1.7), xlabel="横向 (µm)", ylabel="高度 (µm)",
               title=f"{port} 端口：面内电场占比 {te_fraction:.2%}")
    fig.savefig(a.directory / "交叉器_端口模式核对.png", dpi=180)
    plt.close(fig)
    (a.directory / "端口模式核对.json").write_text(json.dumps(metrics, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(metrics, ensure_ascii=False))


if __name__ == "__main__":
    main()
