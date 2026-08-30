"""Reproduce the compact fully custom EO-comb/SHG area budget.

This script only checks first-order scale estimates. It is not an optical
eigenmode, nonlinear-cavity, RF, or foundry-DRC model.
"""

from __future__ import annotations

import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TARGETS = ROOT / "models" / "system" / "design_targets.json"
C_M_PER_S = 299_792_458.0


def wavelength_span_nm(center_nm: float, span_hz: float) -> float:
    """Small-bandwidth conversion from frequency span to wavelength span."""
    center_m = center_nm * 1e-9
    return center_m**2 * span_hz / C_M_PER_S * 1e9


def ring_fsr_hz(radius_um: float, group_index: float) -> float:
    """First-order ring FSR in frequency units."""
    circumference_m = 2.0 * math.pi * radius_um * 1e-6
    return C_M_PER_S / (group_index * circumference_m)


def radius_for_fsr_um(fsr_ghz: float, group_index: float) -> float:
    """Ring radius required for a target FSR in the same approximation."""
    return C_M_PER_S / (2.0 * math.pi * group_index * fsr_ghz * 1e9) * 1e6


def main() -> None:
    data = json.loads(TARGETS.read_text(encoding="utf-8"))
    comb = data["custom_eo_comb"]
    shg = data["poling_free_shg"]

    for label, count in (
        ("first_pass", comb["first_pass_target_line_count"]),
        ("advanced", comb["advanced_target_line_count"]),
    ):
        span_hz = (count - 1) * comb["rf_frequency_ghz"] * 1e9
        span_nm = wavelength_span_nm(comb["center_wavelength_nm"], span_hz)
        print(f"{label}_comb_span={span_hz / 1e12:.3f} THz ({span_nm:.2f} nm)")

    print("\nPoling-free SHG ring scale")
    print("radius_um,diameter_um,estimated_fsr_ghz")
    for radius_um in shg["radius_scale_points_um"]:
        fsr_ghz = ring_fsr_hz(radius_um, shg["preliminary_group_index"]) / 1e9
        print(f"{radius_um:.1f},{2.0 * radius_um:.1f},{fsr_ghz:.1f}")

    rf_matched_radius = radius_for_fsr_um(
        comb["rf_frequency_ghz"], shg["preliminary_group_index"]
    )
    print(
        f"\nA ring whose FSR is {comb['rf_frequency_ghz']:.1f} GHz would need "
        f"an estimated radius of {rf_matched_radius:.0f} um."
    )
    print(
        "The selected SHG-before-modulation architecture does not impose this "
        "RF-to-optical-FSR condition."
    )

    reserved_area = 0.0
    print("\nFloorplan reservations")
    print("id,length_mm,width_mm,area_mm2")
    for structure in data["structures"]:
        area = structure["reserved_length_mm"] * structure["reserved_width_mm"]
        reserved_area += area
        print(
            f"{structure['id']},{structure['reserved_length_mm']:.2f},"
            f"{structure['reserved_width_mm']:.2f},{area:.2f}"
        )

    die_area = data["die"]["length_mm"] * data["die"]["width_mm"]
    print(f"reserved_total_mm2={reserved_area:.2f}")
    print(f"effective_die_area_mm2={die_area:.2f}")
    print(f"functional_reservation_fraction={reserved_area / die_area:.1%}")
    print(
        "\nWarning: the PDK is used only for stack/rules/die constraints. All "
        "three devices are custom, and all nonlinear and dual-band values "
        "require simulation and written foundry confirmation."
    )


if __name__ == "__main__":
    main()
