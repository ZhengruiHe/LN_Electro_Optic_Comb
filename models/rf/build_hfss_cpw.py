"""Build the parameterized LN travelling-wave CPW models in HFSS.

The script deliberately refuses a final solve while the foundry stack is still
marked unconfirmed. Use ``check`` without AEDT. Use ``build`` or ``solve`` only
after installing Ansys Electronics Desktop/HFSS. Engineering placeholders can
be used to exercise model creation with the explicit ``--allow-placeholders``
flag; their results are not design data.
"""

from __future__ import annotations

import argparse
import json
import math
import os
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG = Path(__file__).with_name("hfss_cpw_config.json")
RESULTS = ROOT / "results" / "hfss"
MU_0 = 4.0e-7 * math.pi


def load_config(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as stream:
        data = json.load(stream)
    validate_config(data)
    return data


def validate_config(data: dict[str, Any]) -> None:
    cpw = data["cpw_um"]
    grid = data["screening_grid_um"]
    if cpw["signal_ground_gap"] < 3.0:
        raise ValueError("signal_ground_gap violates the 3-um PDK lower bound")
    if min(grid["signal_ground_gap"]) < 3.0:
        raise ValueError("screening grid violates the 3-um PDK lower bound")
    if len(data["materials"]["ln_relative_permittivity_global_xyz"]) != 3:
        raise ValueError("LN RF permittivity must contain three tensor entries")
    if len(cpw["short_line_lengths"]) != 2:
        raise ValueError("exactly two short-line lengths are required")
    if cpw["short_line_lengths"][0] >= cpw["short_line_lengths"][1]:
        raise ValueError("short-line lengths must be increasing")


def installed_aedt_roots() -> dict[str, str]:
    return {
        key: value
        for key, value in os.environ.items()
        if key.startswith("ANSYSEM_ROOT") and value
    }


def print_check(data: dict[str, Any], config_path: Path) -> None:
    try:
        pyaedt_version = version("pyaedt")
    except PackageNotFoundError:
        pyaedt_version = "NOT_INSTALLED"
    roots = installed_aedt_roots()
    confirmed = data["required_foundry_inputs"]["confirmed"]
    print(f"config={config_path}")
    print(f"pyaedt={pyaedt_version}")
    print(f"aedt_roots={roots if roots else 'NOT_FOUND'}")
    print(f"foundry_inputs_confirmed={confirmed}")
    if not confirmed:
        print("missing_foundry_inputs:")
        for item in data["required_foundry_inputs"]["items"]:
            print(f"- {item}")


def add_materials(hfss: Any, cfg: dict[str, Any]) -> dict[str, str]:
    values = cfg["materials"]

    silica = hfss.materials.add_material("Project_SiO2_RF")
    silica.permittivity = values["silica_relative_permittivity"]
    silica.dielectric_loss_tangent = values["silica_loss_tangent"]

    silicon = hfss.materials.add_material("Project_Si_RF")
    silicon.permittivity = values["silicon_relative_permittivity"]
    silicon.dielectric_loss_tangent = values["silicon_loss_tangent"]

    sin = hfss.materials.add_material("Project_SiN_RF")
    sin.permittivity = values["sin_relative_permittivity"]
    sin.dielectric_loss_tangent = values["sin_loss_tangent"]

    ln = hfss.materials.add_material("Project_LN_RF_Anisotropic")
    ln.permittivity = values["ln_relative_permittivity_global_xyz"]
    ln.dielectric_loss_tangent = values["ln_loss_tangent_global_xyz"]

    metal = hfss.materials.add_material("Project_M1_RF")
    metal.conductivity = values["metal_conductivity_s_per_m"]
    metal.permeability = values["metal_relative_permeability"]

    return {
        "silica": silica.name,
        "silicon": silicon.name,
        "sin": sin.name,
        "ln": ln.name,
        "metal": metal.name,
    }


def create_stack_and_cpw(
    hfss: Any,
    cfg: dict[str, Any],
    line_length_um: float,
    signal_width_um: float,
    gap_um: float,
    ground_width_um: float,
) -> dict[str, Any]:
    stack = cfg["stack_um"]
    cpw = cfg["cpw_um"]
    materials = add_materials(hfss, cfg)
    hfss.modeler.model_units = "um"

    width = stack["model_width"]
    x0 = -0.5 * width
    y0 = 0.0
    z = -stack["silicon_handle"]

    solids: dict[str, Any] = {}
    solids["handle"] = hfss.modeler.create_box(
        [x0, y0, z],
        [width, line_length_um, stack["silicon_handle"]],
        name="Si_Handle",
        material=materials["silicon"],
    )
    z = 0.0
    solids["bottom_oxide"] = hfss.modeler.create_box(
        [x0, y0, z],
        [width, line_length_um, stack["bottom_oxide"]],
        name="Bottom_Oxide",
        material=materials["silica"],
    )
    z += stack["bottom_oxide"]
    solids["sin"] = hfss.modeler.create_box(
        [x0, y0, z],
        [width, line_length_um, stack["sin"]],
        name="SiN_Layer",
        material=materials["sin"],
    )
    z += stack["sin"]
    solids["interlayer_oxide"] = hfss.modeler.create_box(
        [x0, y0, z],
        [width, line_length_um, stack["ln_sin_interlayer_oxide"]],
        name="LN_SiN_Interlayer_Oxide",
        material=materials["silica"],
    )
    z += stack["ln_sin_interlayer_oxide"]
    solids["ln"] = hfss.modeler.create_box(
        [x0, y0, z],
        [width, line_length_um, stack["ln"]],
        name="XCut_LN",
        material=materials["ln"],
    )
    z += stack["ln"]
    solids["top_cladding"] = hfss.modeler.create_box(
        [x0, y0, z],
        [width, line_length_um, stack["top_cladding"]],
        name="Top_Cladding",
        material=materials["silica"],
    )
    z += stack["top_cladding"]
    metal_z = z

    signal_x = -0.5 * signal_width_um
    left_ground_x = signal_x - gap_um - ground_width_um
    right_ground_x = 0.5 * signal_width_um + gap_um
    solids["signal"] = hfss.modeler.create_box(
        [signal_x, y0, metal_z],
        [signal_width_um, line_length_um, stack["metal"]],
        name="M1_Signal",
        material=materials["metal"],
    )
    solids["ground_left"] = hfss.modeler.create_box(
        [left_ground_x, y0, metal_z],
        [ground_width_um, line_length_um, stack["metal"]],
        name="M1_Ground_Left",
        material=materials["metal"],
    )
    solids["ground_right"] = hfss.modeler.create_box(
        [right_ground_x, y0, metal_z],
        [ground_width_um, line_length_um, stack["metal"]],
        name="M1_Ground_Right",
        material=materials["metal"],
    )

    region = hfss.modeler.create_region(
        [
            f"{cpw['air_side_padding']}um",
            f"{cpw['air_side_padding']}um",
            "0um",
            "0um",
            f"{cpw['air_top_padding']}um",
            "0um",
        ],
        pad_type="Absolute Offset",
        name="Air_Region",
    )
    solids["region"] = region

    # Follow the official CPW terminal example: radiation only on the top and
    # bottom thickness faces, Perfect H on the lateral/end faces behind ports.
    radiation_faces: list[int] = []
    perfect_h_faces: list[int] = []
    for face in region.faces:
        normal = face.normal
        if normal and abs(normal[2]) > 0.9:
            radiation_faces.append(face.id)
        else:
            perfect_h_faces.append(face.id)
    hfss.assign_radiation_boundary_to_faces(radiation_faces, name="Open_Z")
    hfss.assign_perfecth_to_sheets(perfect_h_faces, name="Open_XY_PerfectH")

    cpw_outer_left = left_ground_x
    cpw_total_width = (
        2.0 * ground_width_um + 2.0 * gap_um + signal_width_um
    )
    port_z0 = -stack["silicon_handle"]
    port_height = (
        stack["silicon_handle"]
        + metal_z
        + stack["metal"]
        + cpw["port_air_height"]
    )
    p1_sheet = hfss.modeler.create_rectangle(
        "XZ",
        [cpw_outer_left, 0.0, port_z0],
        [cpw_total_width, port_height],
        name="PortSheet_1",
        is_covered=True,
    )
    p2_sheet = hfss.modeler.create_rectangle(
        "XZ",
        [cpw_outer_left, line_length_um, port_z0],
        [cpw_total_width, port_height],
        name="PortSheet_2",
        is_covered=True,
    )
    z_line = metal_z + 0.5 * stack["metal"]
    integration_1 = [
        [0.0, 0.0, z_line],
        [right_ground_x + 0.5 * ground_width_um, 0.0, z_line],
    ]
    integration_2 = [
        [0.0, line_length_um, z_line],
        [right_ground_x + 0.5 * ground_width_um, line_length_um, z_line],
    ]
    refs = [solids["ground_left"], solids["ground_right"]]
    hfss.wave_port(
        p1_sheet,
        reference=refs,
        integration_line=integration_1,
        modes=1,
        impedance=50,
        name="P1",
        renormalize=False,
        characteristic_impedance="Zpi",
    )
    hfss.wave_port(
        p2_sheet,
        reference=refs,
        integration_line=integration_2,
        modes=1,
        impedance=50,
        name="P2",
        renormalize=False,
        characteristic_impedance="Zpi",
    )

    f_mesh_hz = cfg["aedt"]["frequency_stop_ghz"] * 1e9
    sigma = cfg["materials"]["metal_conductivity_s_per_m"]
    skin_depth_um = math.sqrt(2.0 / (2.0 * math.pi * f_mesh_hz * MU_0 * sigma)) * 1e6
    hfss.mesh.assign_skin_depth(
        [
            solids["signal"].name,
            solids["ground_left"].name,
            solids["ground_right"].name,
        ],
        skin_depth=f"{skin_depth_um:.4f}um",
        triangulation_max_length=f"{min(gap_um / 2.0, 5.0):.3f}um",
        layers_number="3",
        name="M1_SkinDepth",
    )

    return solids


def add_solution(hfss: Any, cfg: dict[str, Any]) -> str:
    settings = cfg["aedt"]
    setup = hfss.create_setup(name="Setup_RF")
    setup.enable_adaptive_setup_multifrequency(
        [f"{value}GHz" for value in settings["adaptive_frequencies_ghz"]],
        max_delta_s=settings["maximum_delta_s"],
    )
    setup.props["MaximumPasses"] = settings["maximum_passes"]
    setup.props["MinimumConvergedPasses"] = settings[
        "minimum_converged_passes"
    ]
    setup.props["SaveAnyFields"] = True
    setup.update()

    hfss.create_linear_count_sweep(
        setup=setup.name,
        unit="GHz",
        start_frequency=settings["frequency_start_ghz"],
        stop_frequency=settings["frequency_stop_ghz"],
        num_of_freq_points=246,
        name="Sweep_RF",
        save_fields=False,
        sweep_type="Interpolating",
    )
    hfss.create_single_point_sweep(
        setup=setup.name,
        unit="GHz",
        freq=settings["field_frequency_ghz"],
        name="Field_25GHz",
        save_single_field=True,
        save_fields=True,
    )
    return setup.name


def build_or_solve(args: argparse.Namespace, cfg: dict[str, Any]) -> None:
    if not cfg["required_foundry_inputs"]["confirmed"] and not args.allow_placeholders:
        raise RuntimeError(
            "Foundry inputs are not confirmed. Use --allow-placeholders only "
            "to exercise model creation, never for a final result."
        )
    if not installed_aedt_roots():
        raise RuntimeError(
            "No ANSYSEM_ROOTxxx environment variable was found. Install Ansys "
            "Electronics Desktop/HFSS before building the project."
        )

    from ansys.aedt.core import Hfss

    RESULTS.mkdir(parents=True, exist_ok=True)
    cpw = cfg["cpw_um"]
    line_length_um = args.length_um
    design_name = args.design_name or f"CPW_L{line_length_um:g}um"
    project_file = RESULTS / f"{cfg['aedt']['project_name']}.aedt"
    hfss = Hfss(
        project=str(project_file),
        design=design_name,
        solution_type=cfg["aedt"]["solution_type"],
        version=args.aedt_version or cfg["aedt"]["version"],
        non_graphical=args.non_graphical,
        new_desktop=True,
        close_on_exit=True,
    )
    try:
        create_stack_and_cpw(
            hfss,
            cfg,
            line_length_um=line_length_um,
            signal_width_um=args.signal_width_um or cpw["signal_width"],
            gap_um=args.gap_um or cpw["signal_ground_gap"],
            ground_width_um=args.ground_width_um or cpw["ground_width"],
        )
        setup_name = add_solution(hfss, cfg)
        validation_log = RESULTS / f"{design_name}_validation.log"
        validation_code = hfss.validate_simple(validation_log)
        hfss.save_project(str(project_file))
        print(f"project={project_file}")
        print(f"design={design_name}")
        print(f"validation_code={validation_code}")
        print(f"validation_log={validation_log}")
        if args.mode == "solve":
            if validation_code != 1:
                raise RuntimeError("HFSS validation did not pass; solve aborted")
            if not hfss.analyze(setup=setup_name, cores=args.cores):
                raise RuntimeError("HFSS solve failed")
            touchstone = RESULTS / f"{design_name}.s2p"
            hfss.export_touchstone(
                setup=setup_name,
                sweep="Sweep_RF",
                output_file=str(touchstone),
                renormalization=False,
                gamma_impedance_comments=True,
            )
            print(f"touchstone={touchstone}")
    finally:
        hfss.release_desktop(close_projects=True, close_desktop=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=["check", "build", "solve"], nargs="?", default="check")
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--aedt-version")
    parser.add_argument("--design-name")
    parser.add_argument("--length-um", type=float, default=500.0)
    parser.add_argument("--signal-width-um", type=float)
    parser.add_argument("--gap-um", type=float)
    parser.add_argument("--ground-width-um", type=float)
    parser.add_argument("--cores", type=int, default=4)
    parser.add_argument("--non-graphical", action="store_true")
    parser.add_argument("--allow-placeholders", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    cfg = load_config(args.config)
    if args.mode == "check":
        print_check(cfg, args.config)
        return
    build_or_solve(args, cfg)


if __name__ == "__main__":
    main()
