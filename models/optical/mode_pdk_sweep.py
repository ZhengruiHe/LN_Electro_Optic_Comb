"""Run a preliminary 1550-nm MODE FDE sweep for the beta TFLN-on-SiN PDK."""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG = Path(__file__).with_name("mode_pdk_config.json")
RESULTS = ROOT / "results" / "optical"
C_M_PER_S = 299_792_458.0


def load_config(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as stream:
        cfg = json.load(stream)
    stack = cfg["stack_um"]
    waveguide = cfg["waveguide_um"]
    if not math.isclose(
        stack["sin"] + stack["ln_sin_interlayer_oxide"], 0.7
    ):
        raise ValueError("unexpected SiN plus interlayer thickness")
    if max(waveguide["etch_depths"]) >= stack["ln"]:
        raise ValueError("first sweep requires a nonzero residual LN slab")
    sidewall_angle_deg = float(
        waveguide.get("sidewall_angle_deg_from_horizontal", 90.0)
    )
    if not 0.0 < sidewall_angle_deg <= 90.0:
        raise ValueError(
            "LN sidewall angle from horizontal must be in the interval (0, 90]"
        )
    if min(waveguide["widths"]) < 0.3:
        raise ValueError("waveguide width violates the 0.3-um LN PDK rule")
    if waveguide["ln_isolation_width"] < max(waveguide["widths"]):
        raise ValueError("LN isolation width must cover every swept waveguide")
    if not waveguide.get("active_region_sin_fully_removed", False):
        raise ValueError("active modulator cross-section requires full SiN removal")
    return cfg


def zelmon_ln_indices(wavelength_um: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return ordinary and extraordinary congruent-LN refractive indices."""
    wavelength_sq = np.asarray(wavelength_um, dtype=float) ** 2
    no_sq = (
        1.0
        + 2.6734 * wavelength_sq / (wavelength_sq - 0.01764)
        + 1.2290 * wavelength_sq / (wavelength_sq - 0.05914)
        + 12.614 * wavelength_sq / (wavelength_sq - 474.60)
    )
    ne_sq = (
        1.0
        + 2.9804 * wavelength_sq / (wavelength_sq - 0.02047)
        + 0.5981 * wavelength_sq / (wavelength_sq - 0.06660)
        + 8.9543 * wavelength_sq / (wavelength_sq - 416.08)
    )
    return np.sqrt(no_sq), np.sqrt(ne_sq)


def load_lumapi(cfg: dict[str, Any]) -> Any:
    api_path = Path(cfg["lumerical"]["api_path"])
    if not (api_path / "lumapi.py").is_file():
        raise FileNotFoundError(f"lumapi.py was not found under {api_path}")
    sys.path.insert(0, str(api_path))
    import lumapi  # type: ignore

    return lumapi


def ln_index_string(wavelength_um: float) -> str:
    """Return MODE x;y;z indices for X-cut LN at one wavelength."""
    no, ne = zelmon_ln_indices(np.asarray([wavelength_um]))
    # The Ansys x-cut modulator convention aligns the extraordinary axis and
    # r33 with MODE's transverse x direction.  Propagation is along z.
    return f"{ne[0]:.12g};{no[0]:.12g};{no[0]:.12g}"


def add_rect(
    mode: Any,
    name: str,
    x_min_um: float,
    x_max_um: float,
    y_min_um: float,
    y_max_um: float,
    material: str,
    mesh_order: int,
    index: str | None = None,
) -> None:
    mode.addrect()
    mode.set("name", name)
    mode.set("x min", x_min_um * 1e-6)
    mode.set("x max", x_max_um * 1e-6)
    mode.set("y min", y_min_um * 1e-6)
    mode.set("y max", y_max_um * 1e-6)
    mode.set("z span", 1.0e-6)
    mode.set("material", material)
    if index is not None:
        mode.set("index", index)
    mode.set("override mesh order from material database", True)
    mode.set("mesh order", mesh_order)


def add_polygon(
    mode: Any,
    name: str,
    vertices_um: np.ndarray,
    material: str,
    mesh_order: int,
    index: str | None = None,
) -> None:
    """Add one extruded x-y polygon to the 2D Z-normal MODE model."""
    mode.addpoly()
    mode.set("name", name)
    mode.set("vertices", np.asarray(vertices_um, dtype=float) * 1e-6)
    mode.set("z span", 1.0e-6)
    mode.set("material", material)
    if index is not None:
        mode.set("index", index)
    mode.set("override mesh order from material database", True)
    mode.set("mesh order", mesh_order)


def rib_widths_um(
    top_width_um: float,
    etch_depth_um: float,
    sidewall_angle_deg_from_horizontal: float,
) -> tuple[float, float]:
    """Return top and bottom rib widths for a symmetric trapezoid."""
    angle_rad = math.radians(sidewall_angle_deg_from_horizontal)
    lateral_run_um = etch_depth_um / math.tan(angle_rad)
    return top_width_um, top_width_um + 2.0 * lateral_run_um


def build_case(
    mode: Any,
    cfg: dict[str, Any],
    width_um: float,
    etch_depth_um: float,
    include_metal: bool,
    wavelength_um: float,
) -> None:
    mode.switchtolayout()
    mode.deleteall()
    stack = cfg["stack_um"]
    waveguide = cfg["waveguide_um"]
    electrode = cfg["local_electrode_um"]
    materials = cfg["material_models"]
    fde = cfg["fde"]
    half_span = 0.5 * fde["x_span_um"]

    sin_top = stack["sin"]
    ln_bottom = sin_top + stack["ln_sin_interlayer_oxide"]
    ln_top = ln_bottom + stack["ln"]
    top_oxide = ln_top + stack["top_cladding"]
    residual_slab = stack["ln"] - etch_depth_um
    sidewall_angle_deg = float(
        waveguide.get("sidewall_angle_deg_from_horizontal", 90.0)
    )
    rib_top_width_um, rib_bottom_width_um = rib_widths_um(
        width_um,
        etch_depth_um,
        sidewall_angle_deg,
    )
    if rib_bottom_width_um > waveguide["ln_isolation_width"]:
        raise ValueError(
            "trapezoidal LN rib bottom width exceeds the LN isolation region"
        )

    add_rect(
        mode,
        "Bottom_Oxide",
        -half_span,
        half_span,
        fde["y_min_um"],
        0.0,
        materials["sio2"],
        4,
    )
    # The active modulator region contains no residual SiN.  Its original
    # 300-nm vertical interval is refilled with SiO2 across the MODE domain.
    add_rect(
        mode,
        "SiN_Removed_Oxide_Fill",
        -half_span,
        half_span,
        0.0,
        sin_top,
        materials["sio2"],
        3,
    )
    add_rect(
        mode,
        "LN_SiN_Interlayer_Oxide",
        -half_span,
        half_span,
        sin_top,
        ln_bottom,
        materials["sio2"],
        4,
    )
    ln_isolation_top_width_um = waveguide["ln_isolation_width"]
    ln_isolation_bottom_width_um = ln_isolation_top_width_um + (
        2.0 * residual_slab / math.tan(math.radians(sidewall_angle_deg))
    )
    add_polygon(
        mode,
        "LN_Residual_Slab",
        np.asarray(
            [
                [-0.5 * ln_isolation_bottom_width_um, ln_bottom],
                [0.5 * ln_isolation_bottom_width_um, ln_bottom],
                [0.5 * ln_isolation_top_width_um, ln_bottom + residual_slab],
                [-0.5 * ln_isolation_top_width_um, ln_bottom + residual_slab],
            ]
        ),
        "<Object defined dielectric>",
        1,
        ln_index_string(wavelength_um),
    )
    slab_top = ln_bottom + residual_slab
    add_polygon(
        mode,
        "LN_Rib",
        np.asarray(
            [
                [-0.5 * rib_bottom_width_um, slab_top],
                [0.5 * rib_bottom_width_um, slab_top],
                [0.5 * rib_top_width_um, ln_top],
                [-0.5 * rib_top_width_um, ln_top],
            ]
        ),
        "<Object defined dielectric>",
        1,
        ln_index_string(wavelength_um),
    )
    add_rect(
        mode,
        "Top_Oxide",
        -half_span,
        half_span,
        ln_bottom,
        top_oxide,
        materials["sio2"],
        4,
    )
    if include_metal:
        half_gap = 0.5 * electrode["gap"]
        cap_width = electrode["t_cap_width"]
        add_rect(
            mode,
            "M1_Left_T_Cap",
            -half_gap - cap_width,
            -half_gap,
            top_oxide,
            top_oxide + stack["metal"],
            materials["gold"],
            1,
        )
        add_rect(
            mode,
            "M1_Right_T_Cap",
            half_gap,
            half_gap + cap_width,
            top_oxide,
            top_oxide + stack["metal"],
            materials["gold"],
            1,
        )

    mode.addfde()
    mode.set("solver type", "2D Z normal")
    mode.set("x min", -half_span * 1e-6)
    mode.set("x max", half_span * 1e-6)
    mode.set("y min", fde["y_min_um"] * 1e-6)
    mode.set("y max", fde["y_max_um"] * 1e-6)
    mode.set("x min bc", "PML")
    mode.set("x max bc", "PML")
    mode.set("y min bc", "PML")
    mode.set("y max bc", "PML")
    mode.set("mesh cells x", fde["mesh_cells_x"])
    mode.set("mesh cells y", fde["mesh_cells_y"])
    mode.set("wavelength", wavelength_um * 1e-6)
    mode.set("number of trial modes", fde["trial_modes"])
    if "search_neff" in fde:
        mode.set("use max index", False)
        mode.set("n", fde["search_neff"])
    mode.set("calculate group index", False)
    mode.set("detailed dispersion calculation", False)
    mode.set("simulation temperature", cfg["temperature_c"])


def scalar(value: Any) -> complex:
    return complex(np.asarray(value).reshape(-1)[0])


def mode_localization_metrics(
    mode: Any,
    path: str,
    central_half_width_um: float,
) -> dict[str, float]:
    """Return field-based metrics used to distinguish rib and platform modes."""
    x_um = np.asarray(mode.getdata(path, "x"), dtype=float).reshape(-1) * 1e6
    y_um = np.asarray(mode.getdata(path, "y"), dtype=float).reshape(-1) * 1e6
    fields = [
        np.asarray(mode.getdata(path, component)).squeeze()
        for component in ("Ex", "Ey", "Ez")
    ]
    intensity = sum(np.abs(field) ** 2 for field in fields)
    if intensity.shape != (x_um.size, y_um.size):
        intensity = np.reshape(intensity, (x_um.size, y_um.size))
    weights = np.abs(np.gradient(x_um))[:, None] * np.abs(np.gradient(y_um))[None, :]
    lateral_energy = np.sum(intensity * weights, axis=1)
    total = float(np.sum(lateral_energy))
    central = np.abs(x_um) <= central_half_width_um
    centroid_x_um = float(np.sum(x_um * lateral_energy) / total)
    mirrored_overlap = sum(
        np.sum(field * np.conj(field[::-1, :])) for field in fields
    )
    field_norm = sum(np.sum(np.abs(field) ** 2) for field in fields)
    return {
        "central_energy_fraction": float(np.sum(lateral_energy[central]) / total),
        "peak_x_um": float(x_um[int(np.argmax(lateral_energy))]),
        "rms_x_um": float(
            np.sqrt(np.sum((x_um - centroid_x_um) ** 2 * lateral_energy) / total)
        ),
        "parity_correlation": float(np.real(mirrored_overlap / field_norm)),
    }


def solve_wavelength(
    mode: Any,
    cfg: dict[str, Any],
    width_um: float,
    etch_depth_um: float,
    include_metal: bool,
    wavelength_um: float,
) -> list[dict[str, Any]]:
    print(
        f"mode_step=build width_um={width_um:g} etch_um={etch_depth_um:g} "
        f"sidewall_deg={cfg['waveguide_um'].get('sidewall_angle_deg_from_horizontal', 90.0):g} "
        f"wavelength_nm={wavelength_um * 1000:g}",
        flush=True,
    )
    build_case(
        mode,
        cfg,
        width_um,
        etch_depth_um,
        include_metal,
        wavelength_um,
    )
    print("mode_step=findmodes", flush=True)
    found = int(mode.findmodes())
    print(f"mode_step=findmodes_done found={found}", flush=True)
    candidates: list[dict[str, Any]] = []
    for mode_index in range(1, found + 1):
        path = f"FDE::data::mode{mode_index}"
        te_fraction = float(scalar(mode.getdata(path, "TE polarization fraction")).real)
        if te_fraction < cfg["fde"]["te_fraction_min"]:
            continue
        neff = scalar(mode.getdata(path, "neff"))
        loss_db_per_m = float(scalar(mode.getdata(path, "loss")).real)
        localization = mode_localization_metrics(
            mode,
            path,
            float(cfg["fde"].get("central_half_width_um", 1.5)),
        )
        candidates.append(
            {
                "solver_mode": mode_index,
                "neff_real": neff.real,
                "neff_imag": neff.imag,
                "loss_db_per_cm": loss_db_per_m / 100.0,
                "te_fraction": te_fraction,
                **localization,
            }
        )
    candidates.sort(key=lambda item: item["neff_real"], reverse=True)
    return candidates[: cfg["fde"]["modes_to_keep"]]


def solve_case(
    mode: Any,
    cfg: dict[str, Any],
    width_um: float,
    etch_depth_um: float,
    include_metal: bool,
) -> list[dict[str, Any]]:
    center_um = cfg["wavelength_um"]
    step_um = cfg["group_index_wavelength_step_nm"] * 1e-3
    wavelength_results = [
        solve_wavelength(
            mode,
            cfg,
            width_um,
            etch_depth_um,
            include_metal,
            wavelength_um,
        )
        for wavelength_um in (center_um - step_um, center_um, center_um + step_um)
    ]
    keep_count = min(len(rows) for rows in wavelength_results)
    kept: list[dict[str, Any]] = []
    for rank in range(keep_count):
        lower, center, upper = (
            wavelength_results[index][rank] for index in range(3)
        )
        dneff_dlambda = (
            upper["neff_real"] - lower["neff_real"]
        ) / (2.0 * step_um)
        row = dict(center)
        row.update(
            {
                "mode_label": f"TE{rank}",
                "width_um": width_um,
                "etch_depth_um": etch_depth_um,
                "residual_slab_um": cfg["stack_um"]["ln"] - etch_depth_um,
                "sidewall_angle_deg_from_horizontal": float(
                    cfg["waveguide_um"].get(
                        "sidewall_angle_deg_from_horizontal", 90.0
                    )
                ),
                "rib_top_width_um": rib_widths_um(
                    width_um,
                    etch_depth_um,
                    float(
                        cfg["waveguide_um"].get(
                            "sidewall_angle_deg_from_horizontal", 90.0
                        )
                    ),
                )[0],
                "rib_bottom_width_um": rib_widths_um(
                    width_um,
                    etch_depth_um,
                    float(
                        cfg["waveguide_um"].get(
                            "sidewall_angle_deg_from_horizontal", 90.0
                        )
                    ),
                )[1],
                "ln2_top_width_um": cfg["waveguide_um"]["ln_isolation_width"],
                "ln2_bottom_width_um": cfg["waveguide_um"]["ln_isolation_width"]
                + 2.0
                * (cfg["stack_um"]["ln"] - etch_depth_um)
                / math.tan(
                    math.radians(
                        float(
                            cfg["waveguide_um"].get(
                                "sidewall_angle_deg_from_horizontal", 90.0
                            )
                        )
                    )
                ),
                "metal_included": include_metal,
                "ng": center["neff_real"] - center_um * dneff_dlambda,
            }
        )
        kept.append(row)
    return kept


def write_results(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = [
        "width_um",
        "etch_depth_um",
        "residual_slab_um",
        "sidewall_angle_deg_from_horizontal",
        "rib_top_width_um",
        "rib_bottom_width_um",
        "ln2_top_width_um",
        "ln2_bottom_width_um",
        "metal_included",
        "mode_label",
        "solver_mode",
        "neff_real",
        "neff_imag",
        "ng",
        "loss_db_per_cm",
        "te_fraction",
        "central_energy_fraction",
        "peak_x_um",
        "rms_x_um",
        "parity_correlation",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output", type=Path, default=RESULTS / "pdk_1550_mode_sweep.csv")
    parser.add_argument("--include-metal", action="store_true")
    parser.add_argument("--save-projects", action="store_true")
    parser.add_argument("--max-cases", type=int)
    parser.add_argument("--width-um", type=float)
    parser.add_argument("--etch-depth-um", type=float)
    parser.add_argument("--sidewall-angle-deg", type=float)
    parser.add_argument("--ln-isolation-width-um", type=float)
    parser.add_argument("--wavelength-nm", type=float)
    parser.add_argument("--group-index-wavelength-step-nm", type=float)
    parser.add_argument("--mesh-scale", type=float, default=1.0)
    parser.add_argument("--trial-modes", type=int)
    parser.add_argument("--te-fraction-min", type=float)
    parser.add_argument("--modes-to-keep", type=int)
    parser.add_argument("--search-neff", type=float)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    cfg = load_config(args.config)
    if args.width_um is not None:
        cfg["waveguide_um"]["widths"] = [args.width_um]
    if args.etch_depth_um is not None:
        cfg["waveguide_um"]["etch_depths"] = [args.etch_depth_um]
    if args.sidewall_angle_deg is not None:
        if not 0.0 < args.sidewall_angle_deg <= 90.0:
            raise ValueError("侧壁相对水平面的角度必须在0到90度之间")
        cfg["waveguide_um"]["sidewall_angle_deg_from_horizontal"] = (
            args.sidewall_angle_deg
        )
    if args.ln_isolation_width_um is not None:
        if args.ln_isolation_width_um < max(cfg["waveguide_um"]["widths"]):
            raise ValueError("LN2保留区顶宽必须大于波导顶宽")
        cfg["waveguide_um"]["ln_isolation_width"] = args.ln_isolation_width_um
    if args.wavelength_nm is not None:
        if args.wavelength_nm <= 0:
            raise ValueError("中心波长必须大于零")
        cfg["wavelength_um"] = args.wavelength_nm * 1e-3
    if args.group_index_wavelength_step_nm is not None:
        if args.group_index_wavelength_step_nm <= 0:
            raise ValueError("群折射率波长步长必须大于零")
        cfg["group_index_wavelength_step_nm"] = (
            args.group_index_wavelength_step_nm
        )
    if args.mesh_scale <= 0:
        raise ValueError("网格缩放系数必须大于零")
    cfg["fde"]["mesh_cells_x"] = round(
        cfg["fde"]["mesh_cells_x"] * args.mesh_scale
    )
    cfg["fde"]["mesh_cells_y"] = round(
        cfg["fde"]["mesh_cells_y"] * args.mesh_scale
    )
    if args.trial_modes is not None:
        if args.trial_modes < 1:
            raise ValueError("搜索模数必须至少为1")
        cfg["fde"]["trial_modes"] = args.trial_modes
    if args.te_fraction_min is not None:
        if not 0.0 <= args.te_fraction_min <= 1.0:
            raise ValueError("TE偏振分数阈值必须在0到1之间")
        cfg["fde"]["te_fraction_min"] = args.te_fraction_min
    if args.modes_to_keep is not None:
        if args.modes_to_keep < 1:
            raise ValueError("保留模数必须至少为1")
        cfg["fde"]["modes_to_keep"] = args.modes_to_keep
    if args.search_neff is not None:
        if args.search_neff <= 1.0:
            raise ValueError("目标有效折射率必须大于1")
        cfg["fde"]["search_neff"] = args.search_neff
    lumapi = load_lumapi(cfg)
    RESULTS.mkdir(parents=True, exist_ok=True)
    mode = lumapi.MODE(hide=cfg["lumerical"]["hide"])
    rows: list[dict[str, Any]] = []
    try:
        print(f"mode_version={mode.version()}", flush=True)
        cases = [
            (width, etch)
            for etch in cfg["waveguide_um"]["etch_depths"]
            for width in cfg["waveguide_um"]["widths"]
        ]
        if args.max_cases is not None:
            cases = cases[: args.max_cases]
        for width_um, etch_depth_um in cases:
            case_rows = solve_case(
                mode,
                cfg,
                width_um=width_um,
                etch_depth_um=etch_depth_um,
                include_metal=args.include_metal,
            )
            rows.extend(case_rows)
            print(
                f"case=width_{width_um:g}um_etch_{etch_depth_um:g}um "
                f"kept_modes={[(r['mode_label'], round(r['neff_real'], 6), round(r['ng'], 6)) for r in case_rows]}"
            )
            if args.save_projects:
                suffix = "metal" if args.include_metal else "bare"
                sidewall_angle_deg = cfg["waveguide_um"].get(
                    "sidewall_angle_deg_from_horizontal", 90.0
                )
                project = RESULTS / (
                    f"mode_w{width_um:g}_e{etch_depth_um:g}_a{sidewall_angle_deg:g}_"
                    f"ln2w{cfg['waveguide_um']['ln_isolation_width']:g}_"
                    f"wl{cfg['wavelength_um'] * 1000:g}_{suffix}.lms"
                )
                mode.save(str(project))
        write_results(args.output, rows)
        print(f"output={args.output}")
    finally:
        mode.close()


if __name__ == "__main__":
    main()
