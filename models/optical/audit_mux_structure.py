"""只读加载已保存EME工程，导出真实几何并绘制结构图，不运行模式或传播求解。"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as PlotPolygon, Rectangle, Patch
import numpy as np
from shapely.geometry import Polygon, LineString
from shapely.ops import unary_union

from mode_pdk_sweep import DEFAULT_CONFIG, load_config, load_lumapi

ROOT = Path(__file__).resolve().parents[2]
FAMILIES = ("LN2_Residual_Platform", "LN1_Main_Rib", "LN1_Auxiliary_Rib")
COLORS = {FAMILIES[0]: "#aedadd", FAMILIES[1]: "#d89537", FAMILIES[2]: "#d89537"}


def serial(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    return value


def section_interval(obj, x):
    polygon = Polygon(obj["vertices_um"])
    section = polygon.intersection(LineString([(x,-100.),(x,100.)]))
    if section.is_empty or section.geom_type != "LineString":
        raise ValueError(f"截面不连续：{obj['name']} x={x}")
    return section.bounds[1], section.bounds[3]


def inspect_sections(objects):
    sections = []
    for x in (0.,375.,750.):
        entry = {"x_um":x, "families":{}}
        for family in FAMILIES:
            slabs = sorted([o for o in objects if o["family"] == family], key=lambda o:o["z_min_um"])
            zmid, widths, centers = [], [], []
            for obj in slabs:
                left,right = section_interval(obj,x)
                zmid.append((obj["z_min_um"]+obj["z_max_um"])/2)
                widths.append(right-left)
                centers.append((left+right)/2)
            slope, intercept = np.polyfit(zmid,widths,1)
            z0,z1 = slabs[0]["z_min_um"],slabs[-1]["z_max_um"]
            entry["families"][family] = {
                "number_of_slices":len(slabs), "z_bottom_um":z0,"z_top_um":z1,
                "height_um":z1-z0,"nominal_top_width_from_slice_fit_um":float(slope*z1+intercept),
                "center_um":float(np.mean(centers)),
                "fitted_angle_deg_from_horizontal":math.degrees(math.atan2(2.,-slope)),
                "slice_widths_um":widths}
        main,aux = (entry["families"][f] for f in FAMILIES[1:])
        entry["nominal_top_gap_um"] = (aux["center_um"]-aux["nominal_top_width_from_slice_fit_um"]/2
                                       -main["center_um"]-main["nominal_top_width_from_slice_fit_um"]/2)
        entry["height_and_angle_match"] = all(abs(f["height_um"]-.2)<1e-8 and
            abs(f["fitted_angle_deg_from_horizontal"]-70.)<1e-6 for f in entry["families"].values())
        sections.append(entry)
    return sections


def extract(project):
    api = load_lumapi(load_config(DEFAULT_CONFIG))
    before = hashlib.sha256(project.read_bytes()).hexdigest()
    info = {"project":str(project.resolve()), "sha256_before":before, "objects":[]}
    print("只读加载已保存EME工程，不调用run或emepropagate",flush=True)
    try:
        with api.MODE(filename=str(project.resolve()),hide=True) as mode:
            info["software_version"] = mode.version()
            for family in FAMILIES:
                for index in range(1,33):
                    name = f"{family}_slice_{index}"
                    if not mode.getnamednumber(name):
                        break
                    vertices = np.asarray(mode.getnamed(name,"vertices"),dtype=float)*1e6
                    vertices += np.array([float(mode.getnamed(name,"x")),float(mode.getnamed(name,"y"))])*1e6
                    obj = {"family":family,"name":name,"vertices_um":vertices.tolist()}
                    for key in ("z min","z max"):
                        obj[key.replace(" ","_")+"_um"] = float(mode.getnamed(name,key))*1e6
                    for key in ("material","index","mesh order"):
                        obj[key.replace(" ","_")] = serial(mode.getnamed(name,key))
                    info["objects"].append(obj)
            info["eme"] = {}
            for key in ("x min","y","y span","z min","z max","wavelength","background material",
                        "group spans","cells","number of modes for all cell groups","energy conservation"):
                info["eme"][key] = serial(mode.getnamed("EME",key))
    finally:
        info["sha256_after"] = hashlib.sha256(project.read_bytes()).hexdigest()
        info["original_project_unchanged"] = info["sha256_after"] == before
        if not info["original_project_unchanged"]:
            raise RuntimeError("原工程哈希变化，需用户核对")
    if not all(any(o["family"] == f for o in info["objects"]) for f in FAMILIES):
        raise ValueError("工程缺少预期的LN几何对象")
    info["sections"] = inspect_sections(info["objects"])
    info["scope"] = "直接读回原生LMS；图中显示真实分层几何，不是模式图或新求解结果"
    info["limitations"] = ["只核对几何，不代表光学收敛或PDK签核", "上部脊与平台颜色不同但同为连续LN材料",
                           "EME均匀SiO2背景，不显式包含硅衬底、空气界面或金属"]
    return info


def plot(info, output):
    plt.rcParams.update({"font.sans-serif":["Microsoft YaHei","SimHei","DejaVu Sans"],
                         "axes.unicode_minus":False,"font.size":11})
    fig, axes = plt.subplots(4,1,figsize=(11.5,12),gridspec_kw={"height_ratios":[1.5,1,1,1]})
    ax = axes[0]
    for family in FAMILIES:
        top = max((o for o in info["objects"] if o["family"] == family),key=lambda o:o["z_max_um"])
        ax.add_patch(PlotPolygon(top["vertices_um"],facecolor=COLORS[family],edgecolor="#40494d",lw=.8))
    for x,tag in ((0,"A"),(375,"B"),(750,"C")):
        ax.axvline(x,color="#7a4760",ls="--",lw=1)
        ax.text(x,4.15,tag,ha="center",color="#54334c")
    ax.set(xlim=(-15,765),ylim=(-4.1,4.6),xlabel="传播方向 x (µm)",ylabel="横向 y (µm)",
           title="俯视：两条 LN1 脊共用连续 LN2 平台")
    ax.legend(handles=[Patch(facecolor=COLORS[FAMILIES[1]],label="LN1 脊"),
                       Patch(facecolor=COLORS[FAMILIES[0]],label="LN2 残余平台")],
              loc="lower center",bbox_to_anchor=(.5,1.12),ncol=2,frameon=False)
    for ax,section,tag in zip(axes[1:],info["sections"],"ABC"):
        ax.set_facecolor("#edf1f3")
        for family in FAMILIES:
            shapes = []
            for obj in info["objects"]:
                if obj["family"] != family:
                    continue
                left,right = section_interval(obj,section["x_um"])
                rect = Rectangle((left,obj["z_min_um"]),right-left,obj["z_max_um"]-obj["z_min_um"],
                                 facecolor=COLORS[family],edgecolor="none")
                ax.add_patch(rect)
                shapes.append(Polygon([(left,obj["z_min_um"]),(right,obj["z_min_um"]),
                                       (right,obj["z_max_um"]),(left,obj["z_max_um"])]))
            outline = np.asarray(unary_union(shapes).exterior.coords)
            ax.plot(outline[:,0],outline[:,1],color="#40494d",lw=.9)
        main,aux = (section["families"][f] for f in FAMILIES[1:])
        w1,w2 = (r["nominal_top_width_from_slice_fit_um"] for r in (main,aux))
        ax.text(main["center_um"],1.17,f"主脊顶宽 {w1:.2f} µm",ha="center")
        ax.annotate(f"辅助脊顶宽 {w2:.2f} µm",xy=(aux["center_um"],1.1),xytext=(2.7,1.20),
                    ha="center",arrowprops={"arrowstyle":"-","color":"#40494d","lw":.8})
        ax.text(0,.79,"共用 LN2 平台：名义顶宽 7.20 µm",ha="center",va="center")
        ax.text(-3.7,1.2,"SiO2",color="#525e65")
        ax.set(xlim=(-4,4),ylim=(.6,1.34),yticks=[.7,.9,1.1],ylabel="高度 z (µm)",
               title=f"{tag} 截面  x={section['x_um']:.0f} µm；名义脊顶间隙 {section['nominal_top_gap_um']:.2f} µm")
        ax.set_xlabel("横向 y (µm)")
        ax.grid(axis="y",color="#b2bcc2",lw=.5,alpha=.55)
        ax.set_axisbelow(True)
    fig.suptitle("MUX EME 工程结构核查：直接读取保存的 .lms 几何",fontsize=16,y=.986)
    fig.subplots_adjust(left=.085,right=.97,bottom=.075,top=.90,hspace=.73)
    fig.text(.085,.026,"两级各 200 nm，每级 4 个 50 nm 薄层近似 70° 侧壁。各图纵横比例不同；颜色不代表不同材料。",fontsize=10)
    fig.text(.085,.009,"仅显示已保存 EME 几何；均匀 SiO2 背景，未显式建模顶包层/空气界面和硅衬底。",fontsize=10)
    fig.savefig(output,dpi=180,facecolor="white")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project",type=Path,default=ROOT/"results/optical/模式复用器_70度_EME_精细检查.lms")
    parser.add_argument("--output-dir",type=Path,required=True)
    parser.add_argument("--reuse-audit",type=Path,help="仅用已有只读审计JSON重新画图，不再打开MODE")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True,exist_ok=False)
    info = json.loads(args.reuse_audit.read_text(encoding="utf-8")) if args.reuse_audit else extract(args.project)
    audit = args.output_dir/"EME原生工程_结构核查.json"
    audit.write_text(json.dumps(info,ensure_ascii=False,indent=2),encoding="utf-8")
    output = args.output_dir/"MUX_EME_实际几何_俯视与截面.png"
    plot(info,output)
    print(json.dumps({"image":str(output.resolve()),"audit":str(audit.resolve()),
                      "objects":len(info["objects"]),"sections":info["sections"],"eme":info["eme"],
                      "project_unchanged":info["original_project_unchanged"]},ensure_ascii=False,indent=2),flush=True)


if __name__ == "__main__":
    main()
