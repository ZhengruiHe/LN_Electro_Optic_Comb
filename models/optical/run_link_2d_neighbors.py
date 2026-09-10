"""串行求解实际GDS四组近邻的两路输入；不启动MUX本体、三维或容差扫描。

status=1保留为时间窗初筛并继续收集，不冒充收敛。任何异常退出停止队列。
状态、工程、原始结果及分析均写新目录；源GDS/旧工程不写入。
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time
import psutil

ROOT = Path(__file__).resolve().parents[2]


def select_time_extensions(cases, previous, time_ps):
    """只选择已结束批次的status=1；能量停止案例不得被重复送算。"""
    if previous['status'] not in ('neighbors_screening_completed', 'neighbors_energy_stopped'):
        raise ValueError('旧批次尚未正常结束，不开启时间窗补算')
    screened = {e['case'] for e in previous['completed'] if e['fdtd_status'] == 1}
    old_times = {e['name']: e['time_ps'] for e in previous['cases']}
    selected = []
    for case in cases:
        if case['name'] not in screened:
            continue
        if not time_ps > old_times[case['name']]:
            raise ValueError('补算时间窗必须长于旧时间窗')
        selected.append({**case, 'time_ps': time_ps, 'previous_time_ps': old_times[case['name']]})
    if len(selected) != len(screened) or not selected:
        raise ValueError('没有待补算案例，或案例名与旧批次不一致')
    return selected


def read(path):
    return json.loads(path.read_text(encoding="utf-8"))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preparation", type=Path, required=True)
    parser.add_argument("--model", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--screening-batch", type=Path, help="只补算已结束批次中status=1的案例")
    parser.add_argument("--extension-time-ps", type=float, default=16.)
    args = parser.parse_args()
    prep = read(args.preparation / "几何准备核对.json")
    if not prep["all_geometry_pass"] or len(prep["ready"]) != 4:
        raise ValueError("四组实际GDS几何未准备好，不启动")
    model_sha = hashlib.sha256(args.model.read_bytes()).hexdigest()
    ref_state, ref_result = read(args.reference / "状态.json"), read(args.reference / "结果.json")
    if ref_state["model_sha256"] != model_sha or ref_result["fdtd_status"] != 2:
        raise ValueError("直参考材料或能量停止状态不符")
    args.output_dir.mkdir(parents=True, exist_ok=False)
    cases = []
    for pair in prep["ready"]:
        config = read(Path(pair["config"]))
        gds = Path(config["source_gds"])
        if hashlib.sha256(gds.read_bytes()).hexdigest() != config["source_gds_sha256"]:
            raise ValueError("几何提取后源GDS已变化")
        for source in ("far_a", "far_b"):
            cases.append({"name": pair["pair"] + "_" + source, "pair_config": pair["config"],
                          "source": source, "time_ps": 10 if pair["pair"] in ("left_upper", "right_lower") else 8})
    previous = None
    if args.screening_batch:
        previous = read(args.screening_batch / '批次状态.json')
        if previous['model_sha256'] != model_sha or previous['source_gds_sha256'] != prep['gds_sha256']:
            raise ValueError('补算材料或版图与旧批次不同，不能称仅延长时间窗')
        cases = select_time_extensions(cases, previous, args.extension_time_ps)
        for case in cases:
            old = read(args.screening_batch / case['name'] / '参数与几何.json')
            current = read(Path(case['pair_config']))
            for key in ('core_polygons_um', 'slab_polygons_um', 'ports', 'bounds_um'):
                if old['pair'][key] != current[key]:
                    raise ValueError('补算几何/端口不一致：' + key)
    state = {"status": "starting", "controller_pid": psutil.Process().pid,
             "controller_create_time": psutil.Process().create_time(), "child_pid": None,
             "cases": cases, "total": len(cases), "completed": [], "model_sha256": model_sha,
             "source_gds_sha256": prep["gds_sha256"], "dimension": "2D", "wavelength_nm": 1550,
             "scope": "四组近邻两路输入的等效二维验证；不包括MUX本体、三维和黑盒内部",
             "full_chain_pass": False, "manufacturing_tolerance_sweep": False,
             "screening_batch": str(args.screening_batch.resolve()) if args.screening_batch else None}
    def save(**values):
        state.update(values, updated_at=time.strftime("%Y-%m-%d %H:%M:%S"))
        temporary = args.output_dir / "批次状态.tmp"
        temporary.write_text(json.dumps(state, ensure_ascii=False, indent=2), encoding="utf-8")
        temporary.replace(args.output_dir / "批次状态.json")
    def invoke(command, log_path):
        with log_path.open("x", encoding="utf-8") as log:
            child = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT, cwd=ROOT)
            save(child_pid=child.pid, child_create_time=psutil.Process(child.pid).create_time())
            code = child.wait()
        save(child_pid=None)
        if code:
            raise RuntimeError(f"子程序返回{code}；输出保存在{log_path.name}")
    save()
    try:
        for index, case in enumerate(cases, 1):
            directory = args.output_dir / case["name"]
            command = [sys.executable, "-X", "utf8", str(ROOT / "models/optical/link_2d_component.py"),
                       "--model", str(args.model.resolve()), "--output-dir", str(directory.resolve()),
                       "--kind", "near_pair", "--pair-config", case["pair_config"],
                       "--source-port", case["source"], "--mesh-nm", "100", "--time-ps", str(case["time_ps"])]
            save(status="running", current_index=index, current_case=case["name"])
            invoke(command, args.output_dir / (case["name"] + "_求解控制器.log"))
            result = read(directory / "结果.json")
            save(status="analyzing")
            analysis = args.output_dir / (case["name"] + "_分析")
            invoke([sys.executable, "-X", "utf8", str(ROOT / "models/optical/analyze_link_2d_pair.py"),
                    "--case", str(directory.resolve()), "--reference", str(args.reference.resolve()),
                    "--output-dir", str(analysis.resolve())], args.output_dir / (case["name"] + "_分析控制器.log"))
            metrics = read(analysis / "近邻局域模式投影.json")
            if not metrics["basis_checks_pass"]:
                raise RuntimeError("局域模式投影未通过基底核对，不继续传播不可靠的模式标签")
            comparison = None
            if previous:
                old_entry = next(e for e in previous['completed'] if e['case'] == case['name'])
                old_metrics = read(Path(old_entry['analysis']))
                delta_il = metrics['reference_normalized_target_insertion_loss_db'] - old_metrics['reference_normalized_target_insertion_loss_db']
                delta_xt = metrics['local_crosstalk_relative_db'] - old_metrics['local_crosstalk_relative_db']
                comparison = {
                    'previous_time_ps': case['previous_time_ps'], 'new_max_time_ps': case['time_ps'],
                    'previous_fdtd_status': old_metrics['fdtd_status'], 'new_fdtd_status': metrics['fdtd_status'],
                    'delta_normalized_insertion_loss_db': delta_il, 'delta_crosstalk_db': delta_xt,
                    'time_window_metrics_stable': abs(delta_il) <= .01 and abs(delta_xt) <= .5,
                    'criteria': '本轮工程数值检查：插损变化不超过0.01dB、串扰变化不超过0.5dB；不是PDK指标或网格收敛',
                    'geometry_and_material_unchanged': True,
                }
                with (analysis / '时间窗比较.json').open('x', encoding='utf-8') as stream:
                    json.dump(comparison, stream, ensure_ascii=False, indent=2)
            state["completed"].append({"case": case["name"], "fdtd_status": result["fdtd_status"],
                                       "analysis": str((analysis / "近邻局域模式投影.json").resolve()),
                                       "target_insertion_loss_db": metrics["local_target_insertion_loss_db"],
                                       "crosstalk_relative_db": metrics["local_crosstalk_relative_db"],
                                       "basis_checks_pass": metrics["basis_checks_pass"],
                                       "time_window_comparison": comparison})
            save()
        statuses = [v["fdtd_status"] for v in state["completed"]]
        save(status="neighbors_energy_stopped" if all(v == 2 for v in statuses) else "neighbors_screening_completed",
             current_case=None, pending=["MUX实际外接续", "必要的同网格参考及数值一致性复核"],
             note="求解结束/能量停止不是全链路通过；未执行工艺容差")
    except Exception as exc:
        save(status="failed", error=str(exc))
        raise


if __name__ == "__main__":
    main()
