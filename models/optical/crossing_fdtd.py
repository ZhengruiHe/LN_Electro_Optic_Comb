"""按当前PDK截面建立自定义渐扩LN交叉器，运行三维矢量FDTD。

借鉴Zhang 2023交叉器的渐扩轮廓与双晶向验证方法；余弦轮廓及尺寸
是本项目设计，不是论文原始尺寸，也不是代工方黑盒真实几何。
所有端口直接接0.7 um路由；参数扫描是名义设计优化，不是制造容差。
"""

from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path

import numpy as np
from shapely.affinity import rotate
from shapely.geometry import Polygon, box
from shapely.ops import unary_union

from mode_pdk_sweep import DEFAULT_CONFIG, load_config, load_lumapi, zelmon_ln_indices
from crossing_pdk import validate_nominal

ROOT = Path(__file__).resolve().parents[2]
MODEL_VERSION = "crossing_exact_frequency_v2"
GEOMETRY_VERSION = "crossing_local_ln2_flat_center_v2"
C0 = 299792458.0


def arm_polygon(half_length: float, center_width: float, port_width: float,
                extent: float, samples: int = 161) -> Polygon:
    s = np.linspace(-half_length, half_length, samples)
    width = port_width + (center_width - port_width) * np.cos(np.pi * s / (2 * half_length)) ** 2
    s = np.concatenate(([-extent], s, [extent]))
    width = np.concatenate(([port_width], width, [port_width]))
    return Polygon(np.vstack((np.column_stack((s, -width / 2)),
                              np.column_stack((s[::-1], width[::-1] / 2)))))


def slab_arm_polygon(args, extent):
    center = getattr(args, "slab_center_width_um", 7.2)
    if not math.isfinite(center) or center <= 0:
        raise ValueError("LN2平台尺寸必须为有限正数")
    if center == 7.2:
        return box(-extent, -3.6, extent, 3.6)
    flat = getattr(args, "slab_flat_half_length_um", 6.4)
    if not math.isfinite(flat) or not max(center, 7.2)/2 < flat < args.half_length_um:
        raise ValueError("LN2等宽中心必须覆盖交汇区且小于渐变半长")
    # 交汇处保持正交直边，避免直接交叉两条渐缩曲线形成尖锐刻蚀缝。
    samples = np.unique(np.r_[np.linspace(-args.half_length_um,args.half_length_um,161), -flat, flat])
    phase = np.clip((np.abs(samples)-flat)/(args.half_length_um-flat),0,1)
    widths = 7.2 + (center-7.2)*np.cos(np.pi*phase/2)**2
    s = np.r_[-extent, samples, extent]
    widths = np.r_[7.2,widths,7.2]
    return Polygon(np.vstack((np.column_stack((s,-widths/2)),
                              np.column_stack((s[::-1],widths[::-1]/2)))))


def geometry(args, extent: float) -> tuple[Polygon, Polygon]:
    horizontal = arm_polygon(args.half_length_um, args.width_x_um, args.port_width_um, extent)
    vertical = rotate(arm_polygon(args.half_length_um, args.width_y_um, args.port_width_um, extent),
                      90, origin=(0, 0))
    # 端口维持原7.2 um平台；只在交叉区平滑改变LN2，不改变端口本征模。
    slab_x = slab_arm_polygon(args, extent)
    slab_y = rotate(slab_x, 90, origin=(0, 0))
    if args.reference_axis == "x":
        return (box(-extent, -args.port_width_um/2, extent, args.port_width_um/2),
                box(-extent, -3.6, extent, 3.6))
    if args.reference_axis == "y":
        return (box(-args.port_width_um/2, -extent, args.port_width_um/2, extent),
                box(-3.6, -extent, 3.6, extent))
    core, slab = unary_union([horizontal, vertical]), unary_union([slab_x, slab_y])
    grid_nm = getattr(args, "mask_grid_nm", 0.)
    if grid_nm:
        grid = grid_nm*1e-3
        core = Polygon(np.round(np.asarray(core.exterior.coords)/grid)*grid)
        slab = Polygon(np.round(np.asarray(slab.exterior.coords)/grid)*grid)
    return core, slab


def assert_wavelength(dataset, target_nm, label):
    actual_nm = np.asarray(dataset["lambda"], dtype=float).ravel()*1e9
    if actual_nm.size != 1 or not np.all(np.isfinite(actual_nm)) or not np.allclose(actual_nm, target_nm, rtol=0, atol=1e-6):
        raise ValueError(f"{label}实际波长{actual_nm.tolist()}nm不等于目标{target_nm}nm")
    return float(actual_nm[0])


def configure_exact_frequency(f, target_nm):
    """端口沿用源频率范围；源必须以频率对称，而非波长对称。

    在所有对象创建后设置，避免波长跨度转换/对象默认设置改写采样中心。
    """
    fc = C0/(target_nm*1e-9)
    f.setglobalsource("set frequency", True)
    f.setglobalsource("frequency start", fc-10e12)
    f.setglobalsource("frequency stop", fc+10e12)
    f.setglobalmonitor("use source limits", False)
    f.setglobalmonitor("frequency span", 0.)
    f.setglobalmonitor("frequency center", fc)
    f.setglobalmonitor("frequency points", 1)
    f.setnamed("FDTD::ports", "monitor frequency points", 1)
    for label, actual in (("source", f.getglobalsource("center frequency")),
                          ("global_monitor", f.getglobalmonitor("frequency center"))):
        if not math.isclose(float(actual), fc, rel_tol=1e-12):
            raise ValueError(f"{label}频率设置回读不一致")
    return {"target_nm": target_nm, "source_center_hz": fc, "monitor_center_hz": fc}


def add_layer(f, polygon: Polygon, z0: float, height: float, slices: int,
              angle: float, index: str, name: str) -> None:
    for i in range(slices):
        a = z0 + height*i/slices
        b = z0 + height*(i+1)/slices
        offset = (z0 + height - (a+b)/2) / math.tan(math.radians(angle))
        shape = polygon.buffer(offset, join_style=2)
        if shape.geom_type != "Polygon" or shape.interiors:
            raise ValueError("当前侧壁分层需要无孔单连通多边形")
        f.addpoly()
        f.set("name", f"{name}_{i+1}")
        f.set("vertices", np.asarray(shape.exterior.coords[:-1])*1e-6)
        f.set("z min", a*1e-6)
        f.set("z max", b*1e-6)
        f.set("material", "<Object defined dielectric>")
        f.set("index", index)
        f.set("override mesh order from material database", True)
        f.set("mesh order", 1)


def build(f, args, cfg):
    half_span = args.half_length_um + args.xy_padding_um
    extent = half_span + 2.0
    core, slab = geometry(args, extent)
    no, ne = zelmon_ln_indices(np.asarray([args.wavelength_nm*1e-3]))
    # 固定全局晶轴：仿真x=晶体Y，仿真y=晶体Z，仿真z=晶体X。
    index = f"{no[0]:.12g};{ne[0]:.12g};{no[0]:.12g}"
    f.switchtolayout()
    f.deleteall()
    f.addfdtd()
    for key, value in {"dimension": "3D", "x": 0., "y": 0.,
                       "x span": 2*half_span*1e-6, "y span": 2*half_span*1e-6,
                       "z min": args.z_min_um*1e-6, "z max": args.z_max_um*1e-6,
                       "simulation time": args.time_fs*1e-15,
                       "auto shutoff min": args.shutoff,
                       "mesh accuracy": 2}.items():
        f.set(key, value)
    for axis in "xyz":
        for edge in ("min", "max"):
            f.set(f"{axis} {edge} bc", "PML")
    f.addrect()
    for key, value in {"name": "SiO2_fill_and_cladding", "x": 0., "y": 0.,
                       "x span": 2*extent*1e-6, "y span": 2*extent*1e-6,
                       "z min": min(-3., args.z_min_um-2)*1e-6, "z max": 2.1e-6,
                       "material": cfg["material_models"]["sio2"],
                       "override mesh order from material database": True, "mesh order": 3}.items():
        f.set(key, value)
    add_layer(f, slab, .7, .2, args.slices, 70, index, "LN2_residual")
    add_layer(f, core, .9, .2, args.slices, 70, index, "LN1_rib")
    f.addmesh()
    for key, value in {"name": "LN_mesh", "x": 0., "y": 0.,
                       "x span": 2*half_span*1e-6, "y span": 2*half_span*1e-6,
                       "z min": .5e-6, "z max": 1.3e-6,
                       "dx": args.mesh_nm*1e-9, "dy": args.mesh_nm*1e-9,
                       "dz": args.mesh_z_nm*1e-9}.items():
        f.set(key, value)
    # 改PML距离时保持端口参考面固定，避免把额外传播长度误当边界误差。
    port_offset = args.half_length_um+2.
    ports = {"west": ("x", -port_offset, "Forward"),
             "east": ("x", port_offset, "Backward"),
             "south": ("y", -port_offset, "Forward"),
             "north": ("y", port_offset, "Backward")}
    if args.reference_axis:
        ports = {k: v for k, v in ports.items() if v[0] == args.reference_axis}
    for name, (axis, position, direction) in ports.items():
        f.addport()
        f.set("name", name)
        f.set("injection axis", axis+"-axis")
        f.set("direction", direction)
        f.set(axis, position*1e-6)
        other = "y" if axis == "x" else "x"
        f.set(other, 0.)
        f.set(other+" span", 10e-6)
        f.set("z min", -1e-6)
        f.set("z max", 2.9e-6)
        f.set("mode selection", "fundamental TE mode")
        f.seteigensolver("number of trial modes", 12)
    f.select("FDTD::ports")
    f.set("source port", args.source)
    f.addpower()
    for key, value in {"name": "field_xy", "monitor type": "2D Z-normal", "z": 1e-6,
                       "x": 0., "y": 0., "x span": 2*(half_span-1)*1e-6,
                       "y span": 2*(half_span-1)*1e-6}.items():
        f.set(key, value)
    frequency_settings = configure_exact_frequency(f, args.wavelength_nm)
    return {"model_version": MODEL_VERSION, "geometry_version": GEOMETRY_VERSION,
            "slab_center_width_um": getattr(args, "slab_center_width_um", 7.2),
            "slab_port_width_um": 7.2,
            "slab_flat_half_length_um": getattr(args, "slab_flat_half_length_um", 6.4),
            "slab_profile": "每臂中央等宽平台，外接余弦平方过渡；±L之外接7.2um等宽平台",
            "reference_note": "直波导基准始终为0.7um脊+7.2um平台，包含在交叉器中的平台过渡损耗不扣除",
            "frequency_settings": frequency_settings,
            "core_top_polygon_um": list(core.exterior.coords),
            "slab_top_polygon_um": list(slab.exterior.coords),
            "port_positions_um": ports, "tensor_index_xyz": [float(no[0]), float(ne[0]), float(no[0])],
            "crystal_axes": {"simulation_x": "crystal_Y", "simulation_y": "crystal_Z", "simulation_z": "crystal_X"},
            "material_note": "目标波长Zelmon无损LN；SiO2使用Palik；非代工实测材料",
            "domain_note": "下方SiN按SiO2填满；顶部1um SiO2后为空气；硅衬底在远场氧化层外，未显式纳入",
            "crossing_origin": "自定义余弦渐扩交叉；参考Zhang形态/双晶向验证方法，不是论文尺寸复现"}


def serial(value):
    if isinstance(value, np.ndarray):
        return serial(value.tolist())
    if isinstance(value, complex):
        return {"real": value.real, "imag": value.imag}
    if isinstance(value, dict):
        return {k: serial(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [serial(v) for v in value]
    return value


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--width-x-um", type=float, default=3.0)
    p.add_argument("--width-y-um", type=float, default=3.0)
    p.add_argument("--half-length-um", type=float, default=12.)
    p.add_argument("--port-width-um", type=float, default=.7)
    p.add_argument("--slab-center-width-um", type=float, default=7.2,
                   help="LN2单臂中心顶宽；在±L平滑接回固定7.2um平台，不是整个十字结构宽度")
    p.add_argument("--slab-flat-half-length-um", type=float, default=6.4,
                   help="LN2中央等宽段半长；初扫固定6.4um，避免窄平台交汇形成锐角刻蚀缝")
    p.add_argument("--mesh-nm", type=float, default=100.)
    p.add_argument("--mesh-z-nm", type=float, default=50.)
    p.add_argument("--slices", type=int, default=4)
    p.add_argument("--xy-padding-um", type=float, default=4.)
    p.add_argument("--z-min-um", type=float, default=-1.3)
    p.add_argument("--z-max-um", type=float, default=3.2)
    p.add_argument("--mask-grid-nm", type=float, default=0., help="0为连续轮廓；1为PDK绘制格点量化检查")
    p.add_argument("--wavelength-nm", type=float, default=1550.)
    p.add_argument("--time-fs", type=float, default=1500.)
    p.add_argument("--shutoff", type=float, default=1e-5)
    p.add_argument("--source", choices=["west", "east", "north", "south"], default="west")
    p.add_argument("--reference-axis", choices=["x", "y"], default=None)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--build-only", action="store_true")
    return p.parse_args()


def main():
    args = parse_args()
    if args.xy_padding_um < 4 or args.z_min_um >= -1.0 or args.z_max_um <= 2.9:
        raise ValueError("边界必须覆盖固定端口区域并保留缓冲距离")
    if args.mask_grid_nm < 0 or not math.isfinite(args.mask_grid_nm):
        raise ValueError("掩膜格点必须为有限非负数")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    summary = args.output_dir / "结果.json"
    project = args.output_dir / "crossing.fsp"
    if summary.exists() or project.exists():
        raise FileExistsError("已有结果或工程，使用新输出目录以保留历史")
    cfg = load_config(DEFAULT_CONFIG)
    core, slab = geometry(args, args.half_length_um+args.xy_padding_um+2.)
    pdk_check = validate_nominal(args, cfg, core, slab)
    api = load_lumapi(cfg)
    print("FDTD session opening", flush=True)
    start = time.time()
    with api.FDTD(hide=True) as f:
        version = f.version()
        info = build(f, args, cfg)
        info["pdk_parameter_precheck"] = pdk_check
        # 只使用一项本机资源，6个物理核心；临时设置在退出前恢复。
        old_resources = {}
        try:
            for key, value in (("processes", "1"), ("threads", "6")):
                old_resources[key] = f.getresource("FDTD", 1, key)
                f.setresource("FDTD", 1, key, value)
            f.save(str(project.resolve()))
            (args.output_dir / "参数与几何.json").write_text(
                json.dumps(serial({"args": {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
                                   "version": version, **info}), ensure_ascii=False, indent=2), encoding="utf-8")
            if args.build_only:
                print("Built "+str(project), flush=True)
                return
            print("FDTD running: "+str(project), flush=True)
            f.run()
            print("FDTD returned; collecting modal S parameters", flush=True)
            results = {}
            for name in info["port_positions_um"]:
                object_name = "FDTD::ports::"+name
                sdata = f.getresult(object_name, "S")
                expand = f.getresult(object_name, "expansion for port monitor")
                actual_nm = assert_wavelength(sdata, args.wavelength_nm, name+" S")
                assert_wavelength(expand, args.wavelength_nm, name+" expansion")
                sval = complex(np.asarray(sdata["S"]).ravel()[0])
                neff = f.getresult(object_name, "neff")
                assert_wavelength(neff, args.wavelength_nm, name+" neff")
                results[name] = {"S": sval, "power": abs(sval)**2, "actual_wavelength_nm": actual_nm,
                                 "power_db": 10*math.log10(max(abs(sval)**2, 1e-30)),
                                 "expansion": expand,
                                 "neff": neff}
                profile = f.getresult(object_name, "mode profiles")
                assert_wavelength(profile, args.wavelength_nm, name+" mode profiles")
                np.savez_compressed(args.output_dir / f"端口模式_{name}.npz",
                                    **{k: v for k, v in profile.items() if isinstance(v, np.ndarray)})
            field = f.getresult("field_xy", "E")
            assert_wavelength(field, args.wavelength_nm, "field_xy")
            np.savez_compressed(args.output_dir / "场分布.npz",
                                **{k: v for k, v in field.items() if isinstance(v, np.ndarray)})
            status = f.getresult("FDTD", "status")
            f.save(str(project.resolve()))
            result = {"status": "本机3D FDTD单点求解；是否收敛以独立比较报告为准", "fdtd_status": status,
                      "model_version": MODEL_VERSION, "spectral_validation_pass": True,
                      "source": args.source, "wavelength_nm": args.wavelength_nm,
                      "ports": results, "elapsed_s": time.time()-start, "project": str(project.resolve()),
                      "reference_axis": args.reference_axis, "geometry": info}
            summary.write_text(json.dumps(serial(result), ensure_ascii=False, indent=2), encoding="utf-8")
            print(json.dumps({k: {"power": v["power"], "power_db": v["power_db"]} for k, v in results.items()}, indent=2), flush=True)
            print("status="+str(status)+" elapsed="+str(time.time()-start), flush=True)
        finally:
            for key, value in old_resources.items():
                f.setresource("FDTD", 1, key, value)


if __name__ == "__main__":
    main()
