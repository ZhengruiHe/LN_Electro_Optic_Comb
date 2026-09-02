"""复用EME工程几何，更新物理端口模式并重新求解、导出S矩阵。"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import numpy as np

from mode_mux_eme import extract_s_matrix, parse_mode_numbers, write_s_matrix
from mode_pdk_sweep import DEFAULT_CONFIG, load_config, load_lumapi


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--project", type=Path, required=True)
    parser.add_argument("--port1-modes", type=parse_mode_numbers, required=True)
    parser.add_argument("--port2-modes", type=parse_mode_numbers, required=True)
    parser.add_argument("--s-output", type=Path, required=True)
    parser.add_argument("--save-project", type=Path)
    return parser.parse_args()


def select_port_modes(mode: Any, port_name: str, selected_modes: list[int]) -> None:
    mode.select(f"EME::Ports::{port_name}")
    mode.set("use full simulation span", True)
    mode.set("mode selection", "user select")
    status = int(mode.updateportmodes(np.asarray(selected_modes, dtype=float)))
    if status != 1:
        raise RuntimeError(f"{port_name}端口模式更新失败")


def main() -> None:
    args = parse_args()
    if not args.project.is_file():
        raise FileNotFoundError(f"找不到EME工程：{args.project}")

    cfg = load_config(args.config)
    lumapi = load_lumapi(cfg)
    mode = lumapi.MODE(filename=str(args.project.resolve()), hide=bool(cfg["lumerical"]["hide"]))
    try:
        print(f"mode_version={mode.version()}", flush=True)
        # MODE 2023 R2不允许在分析态修改端口；回到布局态会清除旧结果，
        # 因此端口更新后必须重新计算EME各单元模式，不能只调用emepropagate。
        mode.switchtolayout()
        select_port_modes(mode, "port_1", args.port1_modes)
        select_port_modes(mode, "port_2", args.port2_modes)
        print("eme_step=calculate_modes", flush=True)
        mode.run()
        print("eme_step=propagate", flush=True)
        mode.emepropagate()

        result: Any | None = None
        for result_name in ("power normalized user s matrix", "user s matrix"):
            try:
                result = mode.getresult("EME", result_name)
                print(f"eme_s_result={result_name}", flush=True)
                break
            except Exception:
                continue
        if result is None:
            raise RuntimeError("EME没有返回用户S矩阵，请先完成一次EME模式求解")

        matrix = extract_s_matrix(result)
        write_s_matrix(args.s_output, matrix)
        if args.save_project is not None:
            args.save_project.parent.mkdir(parents=True, exist_ok=True)
            mode.save(str(args.save_project.resolve()))
            print(f"project={args.save_project}", flush=True)
        print(f"s_output={args.s_output}", flush=True)
        if matrix.shape[0] >= 6 and matrix.shape[1] >= 6:
            print(
                "selected_powers="
                + str(
                    {
                        "主TE0直通": float(abs(matrix[3, 0]) ** 2),
                        "主TE1转辅助TE0": float(abs(matrix[4, 1]) ** 2),
                        "辅助TE0转主TE1": float(abs(matrix[5, 2]) ** 2),
                    }
                ),
                flush=True,
            )
    finally:
        mode.close()


if __name__ == "__main__":
    main()
