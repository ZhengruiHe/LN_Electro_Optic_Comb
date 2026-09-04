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
    segmented = data["segmented_t_um"]
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
    if cpw["default_electrode_style"] not in {"regular", "segmented_t"}:
        raise ValueError("default_electrode_style must be regular or segmented_t")
    if not math.isclose(
        segmented["period"],
        segmented["cap_longitudinal_length_r"] + segmented["unit_gap_c"],
    ):
        raise ValueError("segmented period must equal cap length r plus unit gap c")
    for key in ("cap_lateral_width_s", "neck_longitudinal_width_t"):
        if segmented[key] < 2.0:
            raise ValueError(f"{key} violates the 2-um M1 minimum linewidth")
    if segmented["inner_modulation_gap"] < 3.0:
        raise ValueError("segmented inner gap violates the 3-um M1 minimum space")
    if segmented["unit_gap_c"] < 3.0:
        raise ValueError("segmented unit gap violates the 3-um M1 minimum space")
    rules = data["pdk_layout_rules_um"]
    rails = data["optical_rails_um"]
    if rails["ln_waveguide_width"] < rules["ln_min_line"]:
        raise ValueError("LN waveguide violates the PDK minimum linewidth")
    if rails["ln_isolation_width"] < rails["ln_waveguide_width"]:
        raise ValueError("LN isolation width must cover the LN waveguide")
    if not math.isclose(
        data["stack_um"]["ln1_etch"] + data["stack_um"]["ln_residual_slab"],
        data["stack_um"]["ln"],
    ):
        raise ValueError("LN1 etch depth plus residual slab must equal total LN thickness")
    if not 0.0 < rails["ln_sidewall_angle_deg_from_horizontal"] <= 90.0:
        raise ValueError("LN sidewall angle must be in (0, 90] degrees from horizontal")
    if not rails.get("active_region_sin_fully_removed", False):
        raise ValueError("active modulator cross-section requires full SiN removal")


def validate_effective_geometry(
    cfg: dict[str, Any],
    line_length_um: float,
    signal_width_um: float,
    gap_um: float,
    ground_width_um: float,
    electrode_style: str,
) -> None:
    """Reject dimensions that cannot be drawn under the supplied PDK rules."""
    rules = cfg["pdk_layout_rules_um"]
    segmented = cfg["segmented_t_um"]
    rails = cfg["optical_rails_um"]
    minimum_line = rules["m1_min_line"]
    minimum_space = rules["m1_min_space"]

    if line_length_um <= 0 or line_length_um > rules["effective_design_length"]:
        raise ValueError("RF line length exceeds the 21.8-mm effective design length")
    for name, value in (
        ("signal width", signal_width_um),
        ("ground width", ground_width_um),
    ):
        if value < minimum_line:
            raise ValueError(f"{name} violates the {minimum_line:g}-um M1 minimum linewidth")
    if gap_um < minimum_space:
        raise ValueError(f"RF gap violates the {minimum_space:g}-um M1 minimum spacing")

    base_gap_um = gap_um
    if electrode_style == "segmented_t":
        for name in (
            "cap_lateral_width_s",
            "neck_lateral_length_h",
            "neck_longitudinal_width_t",
            "cap_longitudinal_length_r",
        ):
            if segmented[name] < minimum_line:
                raise ValueError(
                    f"{name} violates the {minimum_line:g}-um M1 minimum linewidth"
                )
        if segmented["unit_gap_c"] < minimum_space:
            raise ValueError(
                f"T-unit gap violates the {minimum_space:g}-um M1 minimum spacing"
            )
        base_gap_um = gap_um + 2.0 * (
            segmented["cap_lateral_width_s"]
            + segmented["neck_lateral_length_h"]
        )

    cpw_width_um = signal_width_um + 2.0 * (
        base_gap_um + ground_width_um
    )
    if cpw_width_um > rules["effective_design_width"]:
        raise ValueError("CPW width exceeds the 3.8-mm effective design width")
    if rails["ln_waveguide_width"] >= gap_um:
        raise ValueError("LN waveguide overlaps the M1 electrodes in the modulation gap")

    optical_metal_clearance_um = 0.5 * (
        gap_um - rails["ln_waveguide_width"]
    )
    print(
        "pdk_geometry_check=pass "
        f"cpw_width_um={cpw_width_um:g} "
        f"optical_metal_clearance_um={optical_metal_clearance_um:g} "
        "active_region_sin_fully_removed=true"
    )


def installed_aedt_roots() -> dict[str, str]:
    roots = {
        key: value
        for key, value in os.environ.items()
        if key.startswith("ANSYSEM_ROOT") and value and Path(value).is_dir()
    }
    # Some AEDT installers leave only a stale environment variable after a
    # side-by-side version change.  Discover local Win64 installations before
    # importing PyAEDT, then expose them using the variable names it expects.
    program_files = Path(os.environ.get("ProgramFiles", r"C:\Program Files"))
    for win64_root in (program_files / "AnsysEM").glob("v[0-9][0-9][0-9]/Win64"):
        suffix = win64_root.parent.name.removeprefix("v")
        key = f"ANSYSEM_ROOT{suffix}"
        roots.setdefault(key, str(win64_root))
        os.environ.setdefault(key, str(win64_root))
    return roots


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
    silicon.conductivity = values["silicon_conductivity_s_per_m"]

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
    electrode_style: str,
    simulation_profile: dict[str, Any],
) -> dict[str, Any]:
    stack = cfg["stack_um"]
    cpw = cfg["cpw_um"]
    segmented = cfg["segmented_t_um"]
    rails = cfg["optical_rails_um"]
    validate_effective_geometry(
        cfg,
        line_length_um,
        signal_width_um,
        gap_um,
        ground_width_um,
        electrode_style,
    )
    materials = add_materials(hfss, cfg)
    hfss.modeler.model_units = "um"

    if electrode_style == "segmented_t":
        cap_width = segmented["cap_lateral_width_s"]
        neck_length = segmented["neck_lateral_length_h"]
        base_gap_um = gap_um + 2.0 * (cap_width + neck_length)
        rail_offset_um = (
            0.5 * signal_width_um
            + neck_length
            + cap_width
            + 0.5 * gap_um
        )
    else:
        base_gap_um = gap_um
        rail_offset_um = 0.5 * (signal_width_um + gap_um)
    rail_centers_um = (-rail_offset_um, rail_offset_um)

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
    # The active modulator region contains no residual SiN.  Its original
    # 300-nm vertical interval is refilled with SiO2 across the RF model width.
    solids["sin_removal_oxide_fill"] = hfss.modeler.create_box(
        [x0, y0, z],
        [width, line_length_um, stack["sin"]],
        name="SiN_Removed_Oxide_Fill",
        material=materials["silica"],
    )
    z += stack["sin"]
    solids["interlayer_oxide"] = hfss.modeler.create_box(
        [x0, y0, z],
        [width, line_length_um, stack["ln_sin_interlayer_oxide"]],
        name="LN_SiN_Interlayer_Oxide",
        material=materials["silica"],
    )
    z += stack["ln_sin_interlayer_oxide"]
    # The modulator cross-section is not a continuous 400-nm LN sheet.  LN1
    # removes the upper 200 nm around each rib, leaving a 200-nm residual slab;
    # LN2 removes the remaining slab outside the local 7.2-um LN region.  Build
    # the resulting 400/200/0-nm pattern explicitly.  Both LN1 and LN2 sidewall
    # angles are measured from the horizontal, matching the optical MODE model.
    ln_bottom_z = z
    slab_height = stack["ln_residual_slab"]
    rib_height = stack["ln1_etch"]
    sidewall_angle_rad = math.radians(
        rails["ln_sidewall_angle_deg_from_horizontal"]
    )
    rib_top_width = rails["ln_waveguide_width"]
    rib_bottom_width = rib_top_width + 2.0 * rib_height / math.tan(
        sidewall_angle_rad
    )
    slab_top_width = rails["ln_isolation_width"]
    slab_bottom_width = slab_top_width + 2.0 * slab_height / math.tan(
        sidewall_angle_rad
    )
    ln_rails: list[Any] = []
    for rail_index, rail_center_um in enumerate(rail_centers_um, start=1):
        slab = hfss.modeler.create_polyline(
            [
                [rail_center_um - 0.5 * slab_bottom_width, y0, ln_bottom_z],
                [rail_center_um + 0.5 * slab_bottom_width, y0, ln_bottom_z],
                [rail_center_um + 0.5 * slab_top_width, y0, ln_bottom_z + slab_height],
                [rail_center_um - 0.5 * slab_top_width, y0, ln_bottom_z + slab_height],
            ],
            cover_surface=True,
            close_surface=True,
            name=f"XCut_LN_Slab_{rail_index}",
            material=materials["ln"],
        )
        slab.sweep_along_vector([0.0, line_length_um, 0.0])
        rib = hfss.modeler.create_polyline(
            [
                [rail_center_um - 0.5 * rib_bottom_width, y0, ln_bottom_z + slab_height],
                [rail_center_um + 0.5 * rib_bottom_width, y0, ln_bottom_z + slab_height],
                [rail_center_um + 0.5 * rib_top_width, y0, ln_bottom_z + stack["ln"]],
                [rail_center_um - 0.5 * rib_top_width, y0, ln_bottom_z + stack["ln"]],
            ],
            cover_surface=True,
            close_surface=True,
            name=f"XCut_LN_Rib_{rail_index}",
            material=materials["ln"],
        )
        rib.sweep_along_vector([0.0, line_length_um, 0.0])
        ln_rails.append(slab.unite([rib]))
    solids["ln"] = ln_rails[0].unite(ln_rails[1:])

    # Fill both LN etch levels and the upper cladding with SiO2, then subtract
    # the patterned LN while retaining it as the dielectric waveguide solid.
    solids["top_cladding"] = hfss.modeler.create_box(
        [x0, y0, ln_bottom_z],
        [width, line_length_um, stack["ln"] + stack["top_cladding"]],
        name="LN_Etch_Fill_And_Top_Cladding",
        material=materials["silica"],
    )
    solids["top_cladding"].subtract([solids["ln"]], keep_originals=True)
    z += stack["ln"]
    z += stack["top_cladding"]
    metal_z = z

    signal_x = -0.5 * signal_width_um
    signal_right_x = 0.5 * signal_width_um
    left_ground_x = signal_x - base_gap_um - ground_width_um
    right_ground_x = signal_right_x + base_gap_um
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

    if electrode_style == "segmented_t":
        cap_width = segmented["cap_lateral_width_s"]
        neck_length = segmented["neck_lateral_length_h"]
        neck_y_width = segmented["neck_longitudinal_width_t"]
        cap_y_length = segmented["cap_longitudinal_length_r"]
        period = segmented["period"]
        left_ground_inner_x = left_ground_x + ground_width_um

        signal_segments: list[Any] = []
        left_ground_segments: list[Any] = []
        right_ground_segments: list[Any] = []
        segment_index = 0
        y_center = 0.5 * period
        while y_center < line_length_um:
            segment_index += 1
            y_cap = y_center - 0.5 * cap_y_length
            y_neck = y_center - 0.5 * neck_y_width

            signal_segments.extend(
                [
                    hfss.modeler.create_box(
                        [signal_right_x, y_neck, metal_z],
                        [neck_length, neck_y_width, stack["metal"]],
                        name=f"T_Signal_Right_Neck_{segment_index}",
                        material=materials["metal"],
                    ),
                    hfss.modeler.create_box(
                        [signal_right_x + neck_length, y_cap, metal_z],
                        [cap_width, cap_y_length, stack["metal"]],
                        name=f"T_Signal_Right_Cap_{segment_index}",
                        material=materials["metal"],
                    ),
                    hfss.modeler.create_box(
                        [signal_x - neck_length, y_neck, metal_z],
                        [neck_length, neck_y_width, stack["metal"]],
                        name=f"T_Signal_Left_Neck_{segment_index}",
                        material=materials["metal"],
                    ),
                    hfss.modeler.create_box(
                        [signal_x - neck_length - cap_width, y_cap, metal_z],
                        [cap_width, cap_y_length, stack["metal"]],
                        name=f"T_Signal_Left_Cap_{segment_index}",
                        material=materials["metal"],
                    ),
                ]
            )
            right_ground_segments.extend(
                [
                    hfss.modeler.create_box(
                        [right_ground_x - neck_length, y_neck, metal_z],
                        [neck_length, neck_y_width, stack["metal"]],
                        name=f"T_Ground_Right_Neck_{segment_index}",
                        material=materials["metal"],
                    ),
                    hfss.modeler.create_box(
                        [right_ground_x - neck_length - cap_width, y_cap, metal_z],
                        [cap_width, cap_y_length, stack["metal"]],
                        name=f"T_Ground_Right_Cap_{segment_index}",
                        material=materials["metal"],
                    ),
                ]
            )
            left_ground_segments.extend(
                [
                    hfss.modeler.create_box(
                        [left_ground_inner_x, y_neck, metal_z],
                        [neck_length, neck_y_width, stack["metal"]],
                        name=f"T_Ground_Left_Neck_{segment_index}",
                        material=materials["metal"],
                    ),
                    hfss.modeler.create_box(
                        [left_ground_inner_x + neck_length, y_cap, metal_z],
                        [cap_width, cap_y_length, stack["metal"]],
                        name=f"T_Ground_Left_Cap_{segment_index}",
                        material=materials["metal"],
                    ),
                ]
            )
            y_center += period

        solids["signal"] = solids["signal"].unite(signal_segments)
        solids["ground_left"] = solids["ground_left"].unite(
            left_ground_segments
        )
        solids["ground_right"] = solids["ground_right"].unite(
            right_ground_segments
        )

    region = hfss.modeler.create_region(
        [
            f"{cpw['air_side_padding']}um",
            f"{cpw['air_side_padding']}um",
            "0um",
            "0um",
            f"{cpw['air_top_padding']}um",
            f"{cpw['air_bottom_padding']}um",
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

    # Keep every edge of each port sheet away from the air-region edges and
    # keep the CPW conductors inside, rather than exactly on, the port outline.
    port_x0 = x0
    port_width = width
    port_bottom_margin = 0.5 * cpw["air_bottom_padding"]
    port_z0 = -stack["silicon_handle"] - port_bottom_margin
    port_height = (
        port_bottom_margin + stack["silicon_handle"]
        + metal_z
        + stack["metal"]
        + cpw["port_air_height"]
    )
    p1_sheet = hfss.modeler.create_rectangle(
        "XZ",
        [port_x0, 0.0, port_z0],
        # PyAEDT's XZ rectangle dimensions are ordered as Z-size, X-size.
        [port_height, port_width],
        name="PortSheet_1",
        is_covered=True,
    )
    p2_sheet = hfss.modeler.create_rectangle(
        "XZ",
        [port_x0, line_length_um, port_z0],
        [port_height, port_width],
        name="PortSheet_2",
        is_covered=True,
    )
    # In a three-conductor CPW, the two separated outer grounds cannot both be
    # passed to AutoIdentify as one reference conductor.  Follow the official
    # HFSS Driven Terminal CPW construction instead: use the centre strip as
    # the reference, identify the two outer strips as terminals, and combine
    # those terminals into a mixed-mode pair.  The pair's common mode is the
    # physical CPW mode (both outer grounds at equal potential).
    p1 = hfss.wave_port(
        p1_sheet,
        reference=solids["signal"],
        modes=1,
        impedance=50,
        name="P1",
        renormalize=False,
        characteristic_impedance="Zpi",
    )
    p2 = hfss.wave_port(
        p2_sheet,
        reference=solids["signal"],
        modes=1,
        impedance=50,
        name="P2",
        renormalize=False,
        characteristic_impedance="Zpi",
    )
    if not p1 or not p2:
        raise RuntimeError("HFSS failed to create one or both CPW wave ports")
    terminal_excitations = list(
        hfss.oboundary.GetExcitationsOfType("Terminal")
    )
    print(f"terminal_excitations={terminal_excitations}")
    for port_name in ("P1", "P2"):
        if not hfss.set_differential_pair(
            assignment=f"{port_name}_T1",
            reference=f"{port_name}_T2",
            common_mode=f"{port_name}_CPW",
            differential_mode=f"{port_name}_Odd",
            common_reference=50.0,
            differential_reference=100.0,
        ):
            raise RuntimeError(
                f"HFSS failed to define the mixed-mode pair for {port_name}; "
                f"available terminals: {terminal_excitations}"
            )

    if simulation_profile["assign_skin_depth_mesh"]:
        f_mesh_hz = cfg["aedt"]["frequency_stop_ghz"] * 1e9
        sigma = cfg["materials"]["metal_conductivity_s_per_m"]
        skin_depth_um = (
            math.sqrt(2.0 / (2.0 * math.pi * f_mesh_hz * MU_0 * sigma))
            * 1e6
        )
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
    else:
        print("metal_mesh=screening_without_forced_skin_depth_seeding")

    return solids


def add_solution(
    hfss: Any, cfg: dict[str, Any], simulation_profile: dict[str, Any]
) -> str:
    settings = cfg["aedt"]
    setup = hfss.create_setup(name="Setup_RF")
    setup.enable_adaptive_setup_multifrequency(
        [
            f"{value}GHz"
            for value in simulation_profile["adaptive_frequencies_ghz"]
        ],
        max_delta_s=simulation_profile["maximum_delta_s"],
    )
    setup.props["MaximumPasses"] = simulation_profile["maximum_passes"]
    setup.props["MinimumConvergedPasses"] = simulation_profile[
        "minimum_converged_passes"
    ]
    setup.props["SaveAnyFields"] = simulation_profile["save_field_sweep"]
    setup.update()
    if simulation_profile.get("matrix_convergence_common_mode", False):
        delta_magnitude = simulation_profile["maximum_delta_s"]
        delta_phase = simulation_profile.get("maximum_delta_phase_deg", 5.0)
        common_mode_entries = [
            ["P1_CPW", "P1_CPW", delta_magnitude, delta_phase],
            ["P1_CPW", "P2_CPW", delta_magnitude, delta_phase],
            ["P2_CPW", "P1_CPW", delta_magnitude, delta_phase],
            ["P2_CPW", "P2_CPW", delta_magnitude, delta_phase],
        ]
        if not setup.use_matrix_convergence(
            entry_selection=2,
            ignore_phase_when_mag_is_less_than=0.01,
            custom_entries=common_mode_entries,
        ):
            raise RuntimeError("HFSS物理CPW共模自定义收敛条件设置失败")
        print(
            "matrix_convergence=physical_CPW_common_mode_only "
            f"delta_mag={delta_magnitude} delta_phase_deg={delta_phase}"
        )

    hfss.create_linear_count_sweep(
        setup=setup.name,
        unit="GHz",
        start_frequency=settings["frequency_start_ghz"],
        stop_frequency=simulation_profile.get(
            "frequency_stop_ghz", settings["frequency_stop_ghz"]
        ),
        num_of_freq_points=simulation_profile["sweep_points"],
        name="Sweep_RF",
        save_fields=False,
        sweep_type="Interpolating",
    )
    if simulation_profile["save_field_sweep"]:
        field_frequency_ghz = float(settings["field_frequency_ghz"])
        field_tag = f"{field_frequency_ghz:g}".replace(".", "p")
        hfss.create_single_point_sweep(
            setup=setup.name,
            unit="GHz",
            freq=field_frequency_ghz,
            name=f"Field_{field_tag}GHz",
            save_single_field=True,
            save_fields=True,
        )
    return setup.name


def build_or_solve(args: argparse.Namespace, cfg: dict[str, Any]) -> None:
    if args.project_name:
        cfg["aedt"]["project_name"] = args.project_name
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

    output_dir = args.results_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    cpw = cfg["cpw_um"]
    segmented = cfg["segmented_t_um"]
    if args.t_neck_length_um is not None:
        segmented["neck_lateral_length_h"] = args.t_neck_length_um
    if args.t_cap_width_um is not None:
        segmented["cap_lateral_width_s"] = args.t_cap_width_um
    if args.t_neck_width_um is not None:
        segmented["neck_longitudinal_width_t"] = args.t_neck_width_um
    if args.t_cap_length_um is not None:
        segmented["cap_longitudinal_length_r"] = args.t_cap_length_um
    if args.t_unit_gap_um is not None:
        segmented["unit_gap_c"] = args.t_unit_gap_um
    segmented["period"] = (
        segmented["cap_longitudinal_length_r"] + segmented["unit_gap_c"]
    )
    if segmented["unit_gap_c"] < 3.0:
        raise ValueError("T 形单元纵向间隔违反 PDK 的 3 µm 最小间距")
    profile_name = args.simulation_profile or cfg["simulation_profiles"]["default"]
    simulation_profile = cfg["simulation_profiles"][profile_name]
    output_suffix = f"_{args.output_tag}" if args.output_tag else ""
    electrode_style = args.electrode_style or cpw["default_electrode_style"]
    line_length_um = args.length_um
    design_name = args.design_name or (
        f"CPW_{electrode_style}_L{line_length_um:g}um"
    )
    project_file = output_dir / f"{cfg['aedt']['project_name']}.aedt"
    hfss = Hfss(
        project=str(project_file),
        design=design_name,
        solution_type=cfg["aedt"]["solution_type"],
        version=args.aedt_version or cfg["aedt"]["version"],
        non_graphical=args.non_graphical,
        new_desktop=True,
        close_on_exit=True,
        remove_lock=args.remove_stale_lock,
    )
    try:
        if args.mode == "refine":
            setup_name = "Setup_RF"
            setup = hfss.get_setup(setup_name)
            setup.enable_adaptive_setup_multifrequency(
                [
                    f"{value}GHz"
                    for value in simulation_profile["adaptive_frequencies_ghz"]
                ],
                max_delta_s=simulation_profile["maximum_delta_s"],
            )
            setup.props["MaximumPasses"] = simulation_profile["maximum_passes"]
            setup.props["MinimumConvergedPasses"] = simulation_profile[
                "minimum_converged_passes"
            ]
            setup.props["SaveAnyFields"] = simulation_profile["save_field_sweep"]
            setup.update()
            if simulation_profile.get("matrix_convergence_common_mode", False):
                delta_magnitude = simulation_profile["maximum_delta_s"]
                delta_phase = simulation_profile.get(
                    "maximum_delta_phase_deg", 5.0
                )
                common_mode_entries = [
                    ["P1_CPW", "P1_CPW", delta_magnitude, delta_phase],
                    ["P1_CPW", "P2_CPW", delta_magnitude, delta_phase],
                    ["P2_CPW", "P1_CPW", delta_magnitude, delta_phase],
                    ["P2_CPW", "P2_CPW", delta_magnitude, delta_phase],
                ]
                if not setup.use_matrix_convergence(
                    entry_selection=2,
                    ignore_phase_when_mag_is_less_than=0.01,
                    custom_entries=common_mode_entries,
                ):
                    raise RuntimeError("HFSS物理CPW共模自定义收敛条件设置失败")
                print(
                    "matrix_convergence=physical_CPW_common_mode_only "
                    f"delta_mag={delta_magnitude} delta_phase_deg={delta_phase}"
                )
            print(
                f"refine_existing_design={design_name} "
                f"simulation_profile={profile_name} "
                f"maximum_passes={simulation_profile['maximum_passes']} "
                f"minimum_converged_passes={simulation_profile['minimum_converged_passes']}"
            )
            if not hfss.analyze(setup=setup_name, cores=args.cores):
                raise RuntimeError("HFSS refinement failed")
            touchstone = output_dir / f"{design_name}{output_suffix}.s4p"
            hfss.export_touchstone(
                setup=setup_name,
                sweep="Sweep_RF",
                output_file=str(touchstone),
                renormalization=False,
                gamma_impedance_comments=True,
            )
            hfss.save_project(str(project_file))
            print(f"touchstone={touchstone}")
            return
        create_stack_and_cpw(
            hfss,
            cfg,
            line_length_um=line_length_um,
            signal_width_um=args.signal_width_um
            or (
                segmented["signal_channel_width"]
                if electrode_style == "segmented_t"
                else cpw["signal_width"]
            ),
            gap_um=args.gap_um
            or (
                segmented["inner_modulation_gap"]
                if electrode_style == "segmented_t"
                else cpw["signal_ground_gap"]
            ),
            ground_width_um=args.ground_width_um
            or (
                segmented["ground_channel_width"]
                if electrode_style == "segmented_t"
                else cpw["ground_width"]
            ),
            electrode_style=electrode_style,
            simulation_profile=simulation_profile,
        )
        setup_name = add_solution(hfss, cfg, simulation_profile)
        validation_log = output_dir / f"{design_name}_validation.log"
        validation_code = hfss.validate_simple(validation_log)
        hfss.save_project(str(project_file))
        print(f"project={project_file}")
        print(f"design={design_name}")
        print(f"simulation_profile={profile_name}")
        print(f"validation_code={validation_code}")
        print(f"validation_log={validation_log}")
        if args.mode == "solve":
            if validation_code != 1:
                raise RuntimeError("HFSS validation did not pass; solve aborted")
            if not hfss.analyze(setup=setup_name, cores=args.cores):
                raise RuntimeError("HFSS solve failed")
            # Each physical GSG port is exported as odd/common mixed modes,
            # therefore the Touchstone network contains four ports.
            touchstone = output_dir / f"{design_name}{output_suffix}.s4p"
            hfss.export_touchstone(
                setup=setup_name,
                sweep="Sweep_RF",
                output_file=str(touchstone),
                renormalization=False,
                gamma_impedance_comments=True,
            )
            # 将自适应网格、收敛历史和扫频结果写回工程，供后续诊断导出。
            hfss.save_project(str(project_file))
            print(f"touchstone={touchstone}")
    finally:
        hfss.release_desktop(close_projects=True, close_desktop=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "mode",
        choices=["check", "build", "solve", "refine"],
        nargs="?",
        default="check",
    )
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=RESULTS,
        help="工程、Touchstone和验证日志的保存目录；扫参应指定results/hfss/扫描结果下的子目录",
    )
    parser.add_argument("--aedt-version")
    parser.add_argument(
        "--project-name",
        help="为局部扫参使用独立AEDT工程名，避免占用当前已打开工程",
    )
    parser.add_argument("--design-name")
    parser.add_argument(
        "--output-tag",
        help="附加到Touchstone文件名的短标签，用于保留不同收敛级别的结果",
    )
    parser.add_argument(
        "--electrode-style", choices=["regular", "segmented_t"]
    )
    parser.add_argument(
        "--simulation-profile",
        choices=[
            "screening",
            "common_mode_phase_verification",
            "common_mode_loss_verification",
            "final",
        ],
    )
    parser.add_argument("--length-um", type=float, default=500.0)
    parser.add_argument("--signal-width-um", type=float)
    parser.add_argument("--gap-um", type=float)
    parser.add_argument("--ground-width-um", type=float)
    parser.add_argument("--t-neck-length-um", type=float)
    parser.add_argument(
        "--t-cap-width-um",
        type=float,
        help="T形帽的横向线宽s；不得小于PDK的M1最小线宽",
    )
    parser.add_argument(
        "--t-neck-width-um",
        type=float,
        help="T形颈沿传播方向的线宽t；用于当前PDK局部扫参",
    )
    parser.add_argument("--t-cap-length-um", type=float)
    parser.add_argument("--t-unit-gap-um", type=float)
    parser.add_argument("--cores", type=int, default=4)
    parser.add_argument("--non-graphical", action="store_true")
    parser.add_argument(
        "--remove-stale-lock",
        action="store_true",
        help="让 PyAEDT 移除目标工程的陈旧锁；仅在确认没有 AEDT 进程时使用",
    )
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
