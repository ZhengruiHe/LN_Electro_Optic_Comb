"""分析 MODE 本征模是否集中在中央 LN 脊波导。

该脚本不按求解器的 mode 编号直接判断单模性，而是读取电场，计算横向
能量分布、峰值位置和中央窗口能量占比，用于区分中央脊模与宽 LN 平台模。
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
LUMAPI = Path(r"C:\Program Files\Lumerical\v232\api\python")


def scalar(value: object) -> complex:
    return complex(np.asarray(value).reshape(-1)[0])


def integration_weights(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float).reshape(-1)
    return np.abs(np.gradient(values))


def analyze(
    project: Path,
    output: Path,
    mode_count: int,
    central_half_width_um: float,
) -> None:
    sys.path.insert(0, str(LUMAPI))
    import lumapi  # type: ignore

    rows: list[dict[str, float | int]] = []
    mode = lumapi.MODE(filename=str(project), hide=True)
    try:
        for mode_index in range(1, mode_count + 1):
            path = f"FDE::data::mode{mode_index}"
            try:
                x_um = np.asarray(mode.getdata(path, "x"), dtype=float).reshape(-1) * 1e6
                y_um = np.asarray(mode.getdata(path, "y"), dtype=float).reshape(-1) * 1e6
                fields = [
                    np.asarray(mode.getdata(path, component)).squeeze()
                    for component in ("Ex", "Ey", "Ez")
                ]
            except Exception:
                break
            intensity = sum(np.abs(field) ** 2 for field in fields)
            if intensity.shape != (x_um.size, y_um.size):
                intensity = np.reshape(intensity, (x_um.size, y_um.size))
            weights = integration_weights(x_um)[:, None] * integration_weights(y_um)[None, :]
            energy = intensity * weights
            lateral_energy = np.sum(energy, axis=1)
            total = float(np.sum(lateral_energy))
            if total <= 0.0:
                continue
            central = np.abs(x_um) <= central_half_width_um
            central_fraction = float(np.sum(lateral_energy[central]) / total)
            centroid_x_um = float(np.sum(x_um * lateral_energy) / total)
            rms_x_um = float(
                np.sqrt(np.sum((x_um - centroid_x_um) ** 2 * lateral_energy) / total)
            )
            peak_x_um = float(x_um[int(np.argmax(lateral_energy))])
            mirrored_overlap = sum(
                np.sum(field * np.conj(field[::-1, :])) for field in fields
            )
            field_norm = sum(np.sum(np.abs(field) ** 2) for field in fields)
            parity_correlation = float(np.real(mirrored_overlap / field_norm))
            rows.append(
                {
                    "solver_mode": mode_index,
                    "neff_real": scalar(mode.getdata(path, "neff")).real,
                    "te_fraction": scalar(
                        mode.getdata(path, "TE polarization fraction")
                    ).real,
                    "loss_db_per_cm": scalar(mode.getdata(path, "loss")).real / 100.0,
                    "central_half_width_um": central_half_width_um,
                    "central_energy_fraction": central_fraction,
                    "peak_x_um": peak_x_um,
                    "centroid_x_um": centroid_x_um,
                    "rms_x_um": rms_x_um,
                    "parity_correlation": parity_correlation,
                }
            )
    finally:
        mode.close()

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"output={output}")
    for row in rows:
        print(
            f"mode={row['solver_mode']} neff={row['neff_real']:.6f} "
            f"TE={row['te_fraction']:.4f} central={row['central_energy_fraction']:.4f} "
            f"peak_x_um={row['peak_x_um']:.3f} rms_x_um={row['rms_x_um']:.3f} "
            f"parity={row['parity_correlation']:.4f}"
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--mode-count", type=int, default=10)
    parser.add_argument("--central-half-width-um", type=float, default=1.5)
    args = parser.parse_args()
    analyze(
        args.project.resolve(),
        args.output,
        args.mode_count,
        args.central_half_width_um,
    )


if __name__ == "__main__":
    main()
