"""用两条不同长度的 HFSS 混合模网络去嵌入传播常数。"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
import skrf as rf


C_M_PER_S = 299_792_458.0


def load_mixed_mode_network(path: Path) -> tuple[rf.Network, rf.Network]:
    network = rf.Network(str(path))
    if network.nports != 4:
        raise ValueError(f"预期 4 端口混合模文件，实际为 {network.nports} 端口: {path}")
    # HFSS 导出顺序：P1_Odd, P1_CPW, P2_Odd, P2_CPW。
    return network, network.subnetwork([1, 3])


def cross_mode_max_db(network: rf.Network, index: int) -> float:
    """返回CPW共模输入耦合到任一Odd端口的最大幅度。"""
    odd_ports = (0, 2)
    common_ports = (1, 3)
    maximum = max(
        abs(network.s[index, odd_output, common_input])
        for odd_output in odd_ports
        for common_input in common_ports
    )
    return 20.0 * np.log10(max(maximum, 1e-30))


def extract(
    short: rf.Network,
    long: rf.Network,
    short_length_um: float,
    long_length_um: float,
    short_full: rf.Network | None = None,
    long_full: rf.Network | None = None,
) -> list[dict[str, float | int]]:
    if not np.allclose(short.f, long.f):
        raise ValueError("两条线的扫频点不一致，不能直接做两长度去嵌入")
    delta_length_m = (long_length_um - short_length_um) * 1e-6
    if delta_length_m <= 0:
        raise ValueError("长线长度必须大于短线长度")

    forward_eigenvalues: list[complex] = []
    for index in range(len(short.f)):
        # 若两条线拥有相同端口夹具，A_long @ inv(A_short) 与额外线长
        # 的传输矩阵相似；相似变换不改变本征值 exp(±gamma*delta_L)。
        delta_abcd = long.a[index] @ np.linalg.inv(short.a[index])
        eigenvalues = np.linalg.eigvals(delta_abcd)
        forward_eigenvalues.append(eigenvalues[np.argmin(np.abs(eigenvalues))])

    eigenvalues = np.asarray(forward_eigenvalues)
    phase_rad = np.unwrap(np.angle(eigenvalues))
    beta_rad_per_m = -phase_rad / delta_length_m
    alpha_np_per_m = -np.log(np.abs(eigenvalues)) / delta_length_m
    n_rf = beta_rad_per_m * C_M_PER_S / (2.0 * np.pi * short.f)
    loss_db_per_cm = 8.685889638 * alpha_np_per_m / 100.0

    rows: list[dict[str, float | int]] = []
    for index, frequency_hz in enumerate(short.f):
        s11 = long.s[index, 0, 0]
        s21 = long.s[index, 1, 0]
        zc = np.lib.scimath.sqrt(long.a[index, 0, 1] / long.a[index, 1, 0])
        if zc.real < 0:
            zc = -zc
        rows.append(
            {
                "frequency_ghz": frequency_hz / 1e9,
                "n_rf": n_rf[index],
                "alpha_np_per_m": alpha_np_per_m[index],
                "loss_db_per_cm": loss_db_per_cm[index],
                "long_s11_db": 20.0 * np.log10(abs(s11)),
                "long_s21_db": 20.0 * np.log10(abs(s21)),
                "long_abcd_zc_real_ohm": zc.real,
                "long_abcd_zc_imag_ohm": zc.imag,
                "propagation_eigenvalue_abs": abs(eigenvalues[index]),
                "propagation_eigenvalue_phase_deg": np.angle(
                    eigenvalues[index], deg=True
                ),
                "short_cross_mode_max_db": (
                    cross_mode_max_db(short_full, index)
                    if short_full is not None
                    else float("nan")
                ),
                "long_cross_mode_max_db": (
                    cross_mode_max_db(long_full, index)
                    if long_full is not None
                    else float("nan")
                ),
                "screening_valid": int(
                    frequency_hz <= 30e9
                    and 1.0 <= n_rf[index] <= 5.0
                    and alpha_np_per_m[index] >= 0.0
                ),
            }
        )
    return rows


def write_csv(path: Path, rows: list[dict[str, float | int]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--short", type=Path, required=True)
    parser.add_argument("--long", type=Path, required=True)
    parser.add_argument("--short-length-um", type=float, required=True)
    parser.add_argument("--long-length-um", type=float, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--target-frequency-ghz", type=float, default=10.0)
    args = parser.parse_args()

    short_full, short = load_mixed_mode_network(args.short)
    long_full, long = load_mixed_mode_network(args.long)
    rows = extract(
        short,
        long,
        args.short_length_um,
        args.long_length_um,
        short_full,
        long_full,
    )
    write_csv(args.output, rows)

    frequency = np.asarray([row["frequency_ghz"] for row in rows])
    target_index = int(np.argmin(np.abs(frequency - args.target_frequency_ghz)))
    at_target = rows[target_index]
    print(f"output={args.output}")
    print(
        f"at_target_{args.target_frequency_ghz:g}GHz="
        f"f={at_target['frequency_ghz']:.3f}GHz, "
        f"nRF={at_target['n_rf']:.6f}, "
        f"loss={at_target['loss_db_per_cm']:.4f}dB/cm, "
        f"S11={at_target['long_s11_db']:.3f}dB, "
        f"S21={at_target['long_s21_db']:.3f}dB, "
        f"cross_mode={at_target['long_cross_mode_max_db']:.1f}dB"
    )


if __name__ == "__main__":
    main()
