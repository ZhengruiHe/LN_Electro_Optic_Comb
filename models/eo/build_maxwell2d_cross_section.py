"""建立10 GHz候选T形电极的二维静电截面并导出电场。

本模型不替代HFSS传输线结果。它只计算单位线电压下，T形帽覆盖区与
主干空隙区在LN波导附近的横向准静电场，供后续光电重叠积分使用。
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG = ROOT / "models" / "rf" / "hfss_cpw_config.json"
RESULTS = ROOT / "results" / "eo"


def installed_aedt_roots() -> dict[str, str]:
    roots = {
        key: value
        for key, value in os.environ.items()
        if key.startswith("ANSYSEM_ROOT") and value and Path(value).is_dir()
    }
    program_files = Path(os.environ.get("ProgramFiles", r"C:\Program Files"))
    for win64_root in (program_files / "AnsysEM").glob("v[0-9][0-9][0-9]/Win64"):
        suffix = win64_root.parent.name.removeprefix("v")
        key = f"ANSYSEM_ROOT{suffix}"
        roots.setdefault(key, str(win64_root))
        os.environ.setdefault(key, str(win64_root))
    return roots


def load_config(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as stream:
        cfg = json.load(stream)
    segmented = cfg["segmented_t_um"]
    stack = cfg["stack_um"]
    rails = cfg["optical_rails_um"]
    if not math.isclose(
        segmented["cap_longitudinal_length_r"] + segmented["unit_gap_c"],
        segmented["period"],
    ):
        raise ValueError("T形周期必须等于帽长与周期空隙之和")
    if not math.isclose(stack["ln1_etch"] + stack["ln_residual_slab"], stack["ln"]):
        raise ValueError("LN1刻蚀深度与残余薄板厚度之和必须等于LN总厚度")
    if not rails["active_region_sin_fully_removed"]:
        raise ValueError("当前有源截面要求下方SiN完全移除并以SiO2回填")
    return cfg


def add_materials(app: Any, cfg: dict[str, Any]) -> dict[str, str]:
    values = cfg["materials"]
    silica = app.materials.add_material("EO_SiO2")
    silica.permittivity = values["silica_relative_permittivity"]
    silicon = app.materials.add_material("EO_HighR_Si")
    silicon.permittivity = values["silicon_relative_permittivity"]
    ln = app.materials.add_material("EO_XCut_LN")
    # Maxwell 2D的X/Y分别为横向/竖直方向；横向对应晶体z轴。
    ln.permittivity = [
        values["ln_relative_permittivity_global_xyz"][0],
        values["ln_relative_permittivity_global_xyz"][2],
        values["ln_relative_permittivity_global_xyz"][1],
    ]
    metal = app.materials.add_material("EO_M1")
    metal.conductivity = values["metal_conductivity_s_per_m"]
    return {
        "silica": silica.name,
        "silicon": silicon.name,
        "ln": ln.name,
        "metal": metal.name,
    }


def create_polygon(app: Any, points: list[list[float]], name: str, material: str) -> Any:
    return app.modeler.create_polyline(
        [[x, y, 0.0] for x, y in points],
        cover_surface=True,
        close_surface=True,
        name=name,
        material=material,
    )


def create_ln_rail(
    app: Any,
    cfg: dict[str, Any],
    center_um: float,
    index: int,
    material: str,
) -> list[Any]:
    stack = cfg["stack_um"]
    rails = cfg["optical_rails_um"]
    ln_bottom = stack["sin"] + stack["ln_sin_interlayer_oxide"]
    slab_height = stack["ln_residual_slab"]
    rib_height = stack["ln1_etch"]
    angle = math.radians(rails["ln_sidewall_angle_deg_from_horizontal"])
    slab_top_width = rails["ln_isolation_width"]
    slab_bottom_width = slab_top_width + 2.0 * slab_height / math.tan(angle)
    rib_top_width = rails["ln_waveguide_width"]
    rib_bottom_width = rib_top_width + 2.0 * rib_height / math.tan(angle)
    slab = create_polygon(
        app,
        [
            [center_um - slab_bottom_width / 2.0, ln_bottom],
            [center_um + slab_bottom_width / 2.0, ln_bottom],
            [center_um + slab_top_width / 2.0, ln_bottom + slab_height],
            [center_um - slab_top_width / 2.0, ln_bottom + slab_height],
        ],
        f"LN_Slab_{index}",
        material,
    )
    rib = create_polygon(
        app,
        [
            [center_um - rib_bottom_width / 2.0, ln_bottom + slab_height],
            [center_um + rib_bottom_width / 2.0, ln_bottom + slab_height],
            [center_um + rib_top_width / 2.0, ln_bottom + stack["ln"]],
            [center_um - rib_top_width / 2.0, ln_bottom + stack["ln"]],
        ],
        f"LN_Rib_{index}",
        material,
    )
    # PyAEDT 1.4.0配合AEDT 2023.1时，会把Maxwell 2D中的各向异性
    # 多边形误判成“不求解内部”；若不显式修正，LN内导出的电场为NaN。
    slab.solve_inside = True
    rib.solve_inside = True
    return [slab, rib]


def add_conductor(
    app: Any,
    x0: float,
    width: float,
    y0: float,
    height: float,
    name: str,
    material: str,
) -> Any:
    return app.modeler.create_rectangle(
        [x0, y0, 0.0], [width, height], name=name, material=material
    )


def build_cross_section(app: Any, cfg: dict[str, Any], section: str) -> dict[str, Any]:
    stack = cfg["stack_um"]
    cpw = cfg["cpw_um"]
    segmented = cfg["segmented_t_um"]
    rails = cfg["optical_rails_um"]
    materials = add_materials(app, cfg)
    app.modeler.model_units = "um"

    model_half_width = 120.0
    substrate_bottom = -60.0
    ln_bottom = stack["sin"] + stack["ln_sin_interlayer_oxide"]
    metal_bottom = ln_bottom + stack["ln"] + stack["top_cladding"]
    metal_top = metal_bottom + stack["metal"]

    substrate = app.modeler.create_rectangle(
        [-model_half_width, substrate_bottom, 0.0],
        [2.0 * model_half_width, substrate_bottom * -1.0 - stack["bottom_oxide"]],
        name="HighR_Si_Handle",
        material=materials["silicon"],
    )
    bottom_oxide = app.modeler.create_rectangle(
        [-model_half_width, -stack["bottom_oxide"], 0.0],
        [2.0 * model_half_width, stack["bottom_oxide"]],
        name="Bottom_Oxide",
        material=materials["silica"],
    )
    lower_fill = app.modeler.create_rectangle(
        [-model_half_width, 0.0, 0.0],
        [2.0 * model_half_width, ln_bottom],
        name="SiN_Removed_And_Interlayer_Oxide",
        material=materials["silica"],
    )

    signal_width = segmented["signal_channel_width"]
    ground_width = segmented["ground_channel_width"]
    h = segmented["neck_lateral_length_h"]
    s = segmented["cap_lateral_width_s"]
    gap = segmented["inner_modulation_gap"]
    base_gap = gap + 2.0 * (h + s)
    signal_left = -0.5 * signal_width
    signal_right = 0.5 * signal_width
    ground_right_inner = signal_right + base_gap
    ground_left_inner = signal_left - base_gap
    rail_offset = signal_right + h + s + 0.5 * gap
    rail_centers = [-rail_offset, rail_offset]

    ln_objects: list[Any] = []
    for index, center in enumerate(rail_centers, start=1):
        ln_objects.extend(create_ln_rail(app, cfg, center, index, materials["ln"]))

    top_oxide = app.modeler.create_rectangle(
        [-model_half_width, ln_bottom, 0.0],
        [2.0 * model_half_width, stack["ln"] + stack["top_cladding"]],
        name="LN_Etch_Fill_And_Top_Cladding",
        material=materials["silica"],
    )
    top_oxide.subtract(ln_objects, keep_originals=True)

    signal_objects = [
        add_conductor(
            app,
            signal_left,
            signal_width,
            metal_bottom,
            stack["metal"],
            "M1_Signal_Trunk",
            materials["metal"],
        )
    ]
    ground_left_x = ground_left_inner - ground_width
    ground_objects = [
        add_conductor(
            app,
            ground_left_x,
            ground_width,
            metal_bottom,
            stack["metal"],
            "M1_Ground_Left_Trunk",
            materials["metal"],
        ),
        add_conductor(
            app,
            ground_right_inner,
            ground_width,
            metal_bottom,
            stack["metal"],
            "M1_Ground_Right_Trunk",
            materials["metal"],
        ),
    ]
    if section == "cap":
        signal_objects.extend(
            [
                add_conductor(
                    app,
                    signal_right + h,
                    s,
                    metal_bottom,
                    stack["metal"],
                    "T_Signal_Right_Cap",
                    materials["metal"],
                ),
                add_conductor(
                    app,
                    signal_left - h - s,
                    s,
                    metal_bottom,
                    stack["metal"],
                    "T_Signal_Left_Cap",
                    materials["metal"],
                ),
            ]
        )
        ground_objects.extend(
            [
                add_conductor(
                    app,
                    ground_right_inner - h - s,
                    s,
                    metal_bottom,
                    stack["metal"],
                    "T_Ground_Right_Cap",
                    materials["metal"],
                ),
                add_conductor(
                    app,
                    ground_left_inner + h,
                    s,
                    metal_bottom,
                    stack["metal"],
                    "T_Ground_Left_Cap",
                    materials["metal"],
                ),
            ]
        )
    elif section != "trunk":
        raise ValueError("截面类型只能是cap或trunk")

    region = app.modeler.create_region(
        [80.0, 80.0, 35.0, 35.0], pad_type="Absolute Offset", name="Air_Region"
    )
    app.assign_voltage(signal_objects, amplitude=1000, name="Signal_1V")
    app.assign_voltage(ground_objects, amplitude=0, name="Ground_0V")
    app.assign_voltage(region.edges, amplitude=0, name="Outer_0V")

    app.mesh.assign_length_mesh(
        [obj.name for obj in ln_objects] + [top_oxide.name],
        inside_selection=True,
        maximum_length="0.08um",
        maximum_elements=None,
        name="Fine_LN_And_Cladding",
    )
    app.mesh.assign_length_mesh(
        [obj.name for obj in signal_objects + ground_objects],
        inside_selection=True,
        maximum_length="0.2um",
        maximum_elements=None,
        name="Fine_Metal",
    )
    setup = app.create_setup(
        name="Setup_EO",
        MaximumPasses=15,
        MinimumPasses=2,
        MinimumConvergedPasses=2,
        PercentError=0.3,
    )
    return {
        "setup": setup,
        "region": region,
        "rail_centers_um": rail_centers,
        "metal_bottom_um": metal_bottom,
        "metal_top_um": metal_top,
        "ln_objects": ln_objects,
        "all_dielectrics": [substrate, bottom_oxide, lower_fill, top_oxide, *ln_objects],
    }


def parse_field_file(path: Path) -> dict[str, np.ndarray]:
    rows: list[list[float]] = []
    number = re.compile(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[Ee][-+]?\d+)?")
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        values = [float(item) for item in number.findall(line)]
        if len(values) >= 6:
            rows.append(values[:6])
    if not rows:
        raise RuntimeError(f"未能解析电场文件：{path}")
    data = np.asarray(rows, dtype=float)
    return {
        "x_m": data[:, 0],
        "y_m": data[:, 1],
        "z_m": data[:, 2],
        "ex_v_per_m": data[:, 3],
        "ey_v_per_m": data[:, 4],
        "ez_v_per_m": data[:, 5],
    }


def export_field(app: Any, section: str, output_dir: Path) -> tuple[Path, Path]:
    # AEDT 2023.1的场导出接口不接受中文文件名，原始场文件保留ASCII名称。
    field_path = output_dir / f"eo_{section}_field_1V.fld"
    x_um = np.arange(-36.0, 36.0 + 0.025, 0.05)
    y_um = np.arange(0.55, 1.30 + 0.0125, 0.025)
    sample_points = [
        [float(x), float(y), 0.0]
        for x in x_um
        for y in y_um
    ]
    exported = app.post.export_field_file(
        quantity="E",
        solution="Setup_EO : LastAdaptive",
        output_file=str(field_path),
        sample_points=sample_points,
        export_with_sample_points=True,
        export_in_si_system=True,
    )
    if not exported:
        raise RuntimeError(f"{section}截面的电场导出失败")
    values = parse_field_file(field_path)
    csv_path = output_dir / f"T形电极_{section}_单位电压电场.csv"
    with csv_path.open("w", newline="", encoding="utf-8-sig") as stream:
        writer = csv.writer(stream)
        writer.writerow(["x_um", "y_um", "Ex_V_per_m", "Ey_V_per_m", "Eabs_V_per_m"])
        for x, y, ex, ey in zip(
            values["x_m"], values["y_m"], values["ex_v_per_m"], values["ey_v_per_m"]
        ):
            writer.writerow([x * 1e6, y * 1e6, ex, ey, math.hypot(ex, ey)])
    return field_path, csv_path


def plot_fields(cfg: dict[str, Any], section_csvs: dict[str, Path], output: Path) -> None:
    plt.rcParams.update(
        {
            "font.sans-serif": ["Microsoft YaHei", "SimHei", "DejaVu Sans"],
            "axes.unicode_minus": False,
            "font.size": 10,
        }
    )
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.0), sharex=True, sharey=True)
    image = None
    for ax, section in zip(axes, ("cap", "trunk")):
        data = np.genfromtxt(section_csvs[section], delimiter=",", names=True, encoding="utf-8-sig")
        x = np.unique(data["x_um"])
        y = np.unique(data["y_um"])
        field = np.abs(data["Ex_V_per_m"]).reshape(x.size, y.size)
        image = ax.pcolormesh(x, y, field.T / 1e5, shading="auto", cmap="viridis")
        ax.set_title("T形帽覆盖区（占周期90%）" if section == "cap" else "主干空隙区（占周期10%）")
        ax.set_xlabel("横向位置（µm）")
        for center in (-30.0, 30.0):
            ax.axvline(center, color="white", linestyle="--", linewidth=0.9, alpha=0.8)
    axes[0].set_ylabel("竖直位置（µm）")
    if image is not None:
        colorbar = fig.colorbar(image, ax=axes.ravel().tolist(), fraction=0.028, pad=0.025)
        colorbar.set_label("横向电场 |Ex|（10⁵ V/m，每1 V线电压）")
    fig.suptitle("10 GHz候选T形电极的二维准静电场（PDK截面）", fontsize=13)
    fig.subplots_adjust(left=0.08, right=0.91, bottom=0.15, top=0.82, wspace=0.10)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=240, facecolor="white")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument(
        "--project",
        type=Path,
        default=RESULTS / "LN_EO_Comb_EO_PDK_10GHz_2023R1_v3.aedt",
    )
    parser.add_argument("--section", choices=["cap", "trunk", "both"], default="both")
    parser.add_argument("--non-graphical", action="store_true", default=True)
    args = parser.parse_args()
    args.project = args.project.resolve()
    cfg = load_config(args.config)
    roots = installed_aedt_roots()
    if not roots:
        raise RuntimeError("没有找到Ansys Electronics Desktop安装目录")
    from ansys.aedt.core import Maxwell2d

    output_dir = args.project.parent
    output_dir.mkdir(parents=True, exist_ok=True)
    sections = ("cap", "trunk") if args.section == "both" else (args.section,)
    section_csvs: dict[str, Path] = {}
    app = None
    try:
        for index, section in enumerate(sections):
            design_name = "EO_Cap_Section" if section == "cap" else "EO_Trunk_Section"
            if index == 0:
                app = Maxwell2d(
                    project=str(args.project),
                    design=design_name,
                    solution_type="ElectrostaticXY",
                    version="2023.1",
                    non_graphical=args.non_graphical,
                    new_desktop=True,
                    close_on_exit=False,
                    remove_lock=True,
                )
            else:
                app.insert_design(design_name, solution_type="ElectrostaticXY")
            build_cross_section(app, cfg, section)
            # 先保存几何和边界，保证场导出失败时已完成的自适应解仍可复用。
            app.save_project(str(args.project))
            app.analyze_setup("Setup_EO", cores=4)
            app.save_project(str(args.project))
            _, csv_path = export_field(app, section, output_dir)
            section_csvs[section] = csv_path
            app.save_project(str(args.project))
        if set(section_csvs) == {"cap", "trunk"}:
            figure = output_dir / "T形电极_单位电压横向电场_截面对比.png"
            plot_fields(cfg, section_csvs, figure)
            print(f"figure={figure}")
        for section, path in section_csvs.items():
            print(f"section={section} field_csv={path}")
        print(f"project={args.project}")
    finally:
        if app is not None:
            app.release_desktop(close_projects=True, close_desktop=True)


if __name__ == "__main__":
    main()
