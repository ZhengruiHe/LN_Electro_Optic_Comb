"""通过KLink对当前四程核心工艺GDS运行名义功能DRC。"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from klink import KLinkClient


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_GDS = (
    ROOT
    / "results"
    / "layout"
    / "四程电光梳_10GHz_15mm_3T_3.5T_3T_名义_核心工艺层_未DRC.gds"
)
DEFAULT_RULES = Path(__file__).with_name("pdk_manual_function_drc.lydrc")
DEFAULT_RDB = ROOT / "results" / "layout" / "四程电光梳_10GHz_名义功能DRC.rdb"
DEFAULT_JSON = ROOT / "results" / "layout" / "四程电光梳_10GHz_名义功能DRC.json"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gds", type=Path, default=DEFAULT_GDS)
    parser.add_argument("--rules", type=Path, default=DEFAULT_RULES)
    parser.add_argument("--rdb", type=Path, default=DEFAULT_RDB)
    parser.add_argument("--json", type=Path, default=DEFAULT_JSON)
    parser.add_argument("--top-cell", default="EO4P_10G_15MM_3T_3P5T_3T_NOMINAL")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8765)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if not args.gds.is_file():
        raise FileNotFoundError(f"找不到待检查GDS：{args.gds}")
    if not args.rules.is_file():
        raise FileNotFoundError(f"找不到DRC规则：{args.rules}")
    args.rdb.parent.mkdir(parents=True, exist_ok=True)
    args.json.parent.mkdir(parents=True, exist_ok=True)

    client = KLinkClient(host=args.host, port=args.port)
    client.connect()
    try:
        result = client.drc_run(
            args.rules.read_text(encoding="utf-8"),
            input_layout=str(args.gds.resolve()),
            output_rdb=str(args.rdb.resolve()),
            top_cell=args.top_cell,
            result_mode="full",
        )
    finally:
        client.close()

    args.json.write_text(
        json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(f"rdb={args.rdb}")
    print(f"json={args.json}")
    print(json.dumps(result, ensure_ascii=False, indent=2))
    if result.get("exception") is not None:
        raise RuntimeError("名义功能DRC脚本执行失败")
    summary = result.get("rdb_summary")
    if not args.rdb.is_file() or not isinstance(summary, dict) or "total_items" not in summary:
        raise RuntimeError("DRC未产生有效报告，不能按0项通过处理")
    if summary["total_items"] != 0:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
