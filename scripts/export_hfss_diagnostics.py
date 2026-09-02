"""从已求解HFSS设计中导出收敛、网格和求解概要。"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_PROJECT = ROOT / "results" / "hfss" / "LN_EO_Comb_RF_2023R1.aedt"
DEFAULT_OUTPUT = ROOT / "results" / "hfss" / "diagnostics"


def discover_aedt_roots() -> None:
    program_files = Path(os.environ.get("ProgramFiles", r"C:\Program Files"))
    for win64_root in (program_files / "AnsysEM").glob("v[0-9][0-9][0-9]/Win64"):
        suffix = win64_root.parent.name.removeprefix("v")
        os.environ.setdefault(f"ANSYSEM_ROOT{suffix}", str(win64_root))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", type=Path, default=DEFAULT_PROJECT)
    parser.add_argument("--design", action="append", required=True)
    parser.add_argument("--setup", default="Setup_RF")
    parser.add_argument("--aedt-version", default="2023.1")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--remove-stale-lock", action="store_true")
    args = parser.parse_args()

    discover_aedt_roots()
    from ansys.aedt.core import Hfss

    args.output_dir.mkdir(parents=True, exist_ok=True)
    for design_name in args.design:
        hfss = Hfss(
            project=str(args.project.resolve()),
            design=design_name,
            version=args.aedt_version,
            non_graphical=True,
            new_desktop=True,
            close_on_exit=True,
            remove_lock=args.remove_stale_lock,
        )
        try:
            # AEDT 2023 R1的底层导出接口会截断含非ASCII字符的扩展名，
            # 因此诊断数据文件名保持ASCII；中文解释写在上层报告中。
            convergence_path = args.output_dir / f"{design_name}_convergence.txt"
            profile_path = args.output_dir / f"{design_name}_profile.prof"
            mesh_path = args.output_dir / f"{design_name}_mesh.txt"
            setup_properties_path = (
                args.output_dir / f"{design_name}_setup_properties.json"
            )
            setup_object = hfss.get_setup(args.setup)
            with setup_properties_path.open("w", encoding="utf-8") as stream:
                json.dump(
                    setup_object.props,
                    stream,
                    ensure_ascii=False,
                    indent=2,
                    default=str,
                )
            hfss.export_convergence(
                setup=args.setup, output_file=str(convergence_path.resolve())
            )
            hfss.export_profile(
                setup=args.setup, output_file=str(profile_path.resolve())
            )
            hfss.export_mesh_stats(
                setup=args.setup, output_file=str(mesh_path.resolve())
            )
            print(f"design={design_name}")
            print(f"setup_is_solved={setup_object.is_solved}")
            print(f"existing_analysis_sweeps={hfss.existing_analysis_sweeps}")
            print(f"nominal_adaptive={hfss.nominal_adaptive}")
            print(f"setup_properties={setup_properties_path}")
            print(f"convergence={convergence_path}")
            print(f"profile={profile_path}")
            print(f"mesh={mesh_path}")
        finally:
            hfss.release_desktop(close_projects=True, close_desktop=True)


if __name__ == "__main__":
    main()
