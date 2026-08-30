"""Reproduce the preliminary EO-comb/PPLN bandwidth and area budget.

The model is intentionally simple. It is a planning check, not a replacement
for anisotropic optical modes, full-wave RF models, or depleted-pump nonlinear
propagation.
"""

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TARGETS = ROOT / "models" / "system" / "design_targets.json"
C_M_PER_S = 299_792_458.0


def wavelength_span_nm(center_nm: float, span_hz: float) -> float:
    """Small-bandwidth conversion from frequency span to wavelength span."""
    center_m = center_nm * 1e-9
    return center_m**2 * span_hz / C_M_PER_S * 1e9


def qpm_fwhm(length_mm: float, gvm_fs_per_mm: float, center_nm: float) -> tuple[float, float]:
    """Uniform-QPM intensity FWHM from the first-order GVM approximation."""
    gvm_s_per_m = gvm_fs_per_mm * 1e-15 / 1e-3
    length_m = length_mm * 1e-3
    bandwidth_hz = 0.442 / (gvm_s_per_m * length_m)
    return bandwidth_hz, wavelength_span_nm(center_nm, bandwidth_hz)


def low_depletion_efficiency(
    normalized_percent_per_w_cm2: float, length_mm: float, pump_w: float
) -> float:
    """Return P_2w/P_w from eta_norm * P_w * L^2."""
    normalized_per_w_cm2 = normalized_percent_per_w_cm2 / 100.0
    length_cm = length_mm / 10.0
    return normalized_per_w_cm2 * pump_w * length_cm**2


def main() -> None:
    data = json.loads(TARGETS.read_text(encoding="utf-8"))
    comb = data["comb"]
    ppln = data["ppln"]

    comb_span_hz = (comb["target_line_count"] - 1) * comb["rf_frequency_ghz"] * 1e9
    comb_span_nm = wavelength_span_nm(comb["center_wavelength_nm"], comb_span_hz)

    print("PRELIMINARY DESIGN BUDGET")
    print(f"EO comb span: {comb_span_hz / 1e12:.3f} THz ({comb_span_nm:.2f} nm)")
    print("\nUniform-PPLN estimates")
    print("length_mm,fwhm_thz,fwhm_nm,ideal_eta_at_50mW,ideal_eta_at_100mW")
    for length_mm in ppln["doe_lengths_mm"]:
        bw_hz, bw_nm = qpm_fwhm(
            length_mm,
            ppln["reference_gvm_fs_per_mm"],
            comb["center_wavelength_nm"],
        )
        eta_50 = low_depletion_efficiency(
            ppln["reference_normalized_efficiency_percent_per_w_cm2"],
            length_mm,
            0.050,
        )
        eta_100 = low_depletion_efficiency(
            ppln["reference_normalized_efficiency_percent_per_w_cm2"],
            length_mm,
            0.100,
        )
        print(
            f"{length_mm:.1f},{bw_hz / 1e12:.3f},{bw_nm:.2f},"
            f"{eta_50 * 100:.1f}%,{eta_100 * 100:.1f}%"
        )

    reserved_area = 0.0
    print("\nFloorplan reservations")
    print("id,length_mm,width_mm,area_mm2")
    for structure in data["structures"]:
        area = structure["reserved_length_mm"] * structure["reserved_width_mm"]
        reserved_area += area
        print(
            f"{structure['id']},{structure['reserved_length_mm']:.1f},"
            f"{structure['reserved_width_mm']:.1f},{area:.1f}"
        )

    die_area = data["die"]["length_mm"] * data["die"]["width_mm"]
    print(f"reserved_total_mm2={reserved_area:.1f}")
    print(f"die_area_mm2={die_area:.1f}")
    print(f"functional_reservation_fraction={reserved_area / die_area:.1%}")
    print(
        "\nWarning: efficiencies above roughly 20% require a depleted-pump model; "
        "all values depend on the final stack and measured loss."
    )


if __name__ == "__main__":
    main()
