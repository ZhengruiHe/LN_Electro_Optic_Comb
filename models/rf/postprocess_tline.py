"""Extract CPW propagation data from two HFSS line-length Touchstone files."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import numpy as np
import skrf as rf


C_M_PER_S = 299_792_458.0


def tline_abcd(gamma: complex, impedance: complex, length_m: float) -> np.ndarray:
    gl = gamma * length_m
    return np.array(
        [
            [np.cosh(gl), impedance * np.sinh(gl)],
            [np.sinh(gl) / impedance, np.cosh(gl)],
        ],
        dtype=complex,
    )


def extract_from_abcd(
    frequency_hz: np.ndarray,
    short_abcd: np.ndarray,
    long_abcd: np.ndarray,
    delta_length_m: float,
) -> dict[str, np.ndarray]:
    if delta_length_m <= 0:
        raise ValueError("delta length must be positive")
    selected: list[complex] = []
    z0_values: list[complex] = []
    for short_matrix, long_matrix in zip(short_abcd, long_abcd):
        delta_matrix = long_matrix @ np.linalg.inv(short_matrix)
        eigenvalues = np.linalg.eigvals(delta_matrix)
        negative_phase = [value for value in eigenvalues if np.angle(value) <= 0]
        if negative_phase:
            value = min(negative_phase, key=lambda item: abs(abs(item) - 1.0))
        else:
            value = min(eigenvalues, key=abs)
        selected.append(value)

        z0 = np.sqrt(delta_matrix[0, 1] / delta_matrix[1, 0])
        if np.real(z0) < 0:
            z0 = -z0
        z0_values.append(z0)

    eigenvalue = np.asarray(selected)
    alpha_np_per_m = -np.log(np.abs(eigenvalue)) / delta_length_m
    beta_rad_per_m = np.unwrap(-np.angle(eigenvalue)) / delta_length_m
    n_rf = C_M_PER_S * beta_rad_per_m / (2.0 * math.pi * frequency_hz)
    loss_db_per_cm = 8.685889638 * alpha_np_per_m * 0.01
    return {
        "frequency_hz": frequency_hz,
        "alpha_np_per_m": alpha_np_per_m,
        "beta_rad_per_m": beta_rad_per_m,
        "loss_db_per_cm": loss_db_per_cm,
        "n_rf": n_rf,
        "z0": np.asarray(z0_values),
    }


def extract_touchstone(
    short_path: Path,
    long_path: Path,
    delta_length_mm: float,
) -> dict[str, np.ndarray]:
    short = rf.Network(str(short_path))
    long = rf.Network(str(long_path))
    if short.nports != 2 or long.nports != 2:
        raise ValueError("both Touchstone files must describe two-port networks")
    if not np.array_equal(short.f, long.f):
        raise ValueError("short and long files must use identical frequency grids")
    return extract_from_abcd(
        short.f,
        short.a,
        long.a,
        delta_length_mm * 1e-3,
    )


def write_csv(path: Path, data: dict[str, np.ndarray]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "frequency_ghz",
                "z0_real_ohm",
                "z0_imag_ohm",
                "n_rf",
                "alpha_np_per_m",
                "loss_db_per_cm",
            ]
        )
        for index, frequency_hz in enumerate(data["frequency_hz"]):
            writer.writerow(
                [
                    frequency_hz / 1e9,
                    np.real(data["z0"][index]),
                    np.imag(data["z0"][index]),
                    data["n_rf"][index],
                    data["alpha_np_per_m"][index],
                    data["loss_db_per_cm"][index],
                ]
            )


def self_test() -> None:
    frequency_hz = np.linspace(1e9, 50e9, 101)
    expected_n_rf = 2.30
    expected_z0 = 50.0
    expected_alpha = 12.0
    short_length_m = 0.5e-3
    long_length_m = 1.0e-3
    short_abcd = []
    long_abcd = []
    for frequency in frequency_hz:
        beta = 2.0 * math.pi * frequency * expected_n_rf / C_M_PER_S
        gamma = expected_alpha + 1j * beta
        short_abcd.append(tline_abcd(gamma, expected_z0, short_length_m))
        long_abcd.append(tline_abcd(gamma, expected_z0, long_length_m))
    result = extract_from_abcd(
        frequency_hz,
        np.asarray(short_abcd),
        np.asarray(long_abcd),
        long_length_m - short_length_m,
    )
    assert np.max(np.abs(result["n_rf"] - expected_n_rf)) < 1e-10
    assert np.max(np.abs(np.real(result["z0"]) - expected_z0)) < 1e-9
    assert np.max(np.abs(result["alpha_np_per_m"] - expected_alpha)) < 1e-9
    print("two_length_extraction_self_test=PASS")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--short", type=Path)
    parser.add_argument("--long", type=Path)
    parser.add_argument("--delta-length-mm", type=float, default=0.5)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--self-test", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.self_test:
        self_test()
        return
    if not args.short or not args.long or not args.output:
        raise SystemExit("--short, --long and --output are required")
    data = extract_touchstone(args.short, args.long, args.delta_length_mm)
    write_csv(args.output, data)
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
