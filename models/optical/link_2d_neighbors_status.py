"""只读核对二维近邻批次及真实进程树；不停止/重启/启动求解器。"""
import argparse
import json
from pathlib import Path
import re
import psutil


def read(path):
    return json.loads(path.read_text(encoding="utf-8"))


def process_info(pid, created=None):
    if not pid:
        return None
    try:
        proc = psutil.Process(pid)
        matching = created is None or abs(proc.create_time() - created) < .01
        children = proc.children(recursive=True) if matching else []
        return {"pid": pid, "alive_matching_identity": matching, "name": proc.name(),
                "cpu_s": sum(proc.cpu_times()[:2]),
                "children": [{"pid": child.pid, "name": child.name(),
                              "cpu_s": sum(child.cpu_times()[:2])} for child in children]}
    except (psutil.NoSuchProcess, psutil.AccessDenied):
        return {"pid": pid, "alive_matching_identity": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batch", type=Path, required=True)
    args = parser.parse_args()
    state = read(args.batch / "批次状态.json")
    cancellation = args.batch / '用户取消补算.json'
    result = {"batch_status": state["status"], "updated_at": state["updated_at"],
              "completed_count": len(state["completed"]), "total": state["total"],
              "current_case": state.get("current_case"),
              "controller": process_info(state["controller_pid"], state.get("controller_create_time")),
              "child": process_info(state.get("child_pid"), state.get("child_create_time")),
              "full_chain_pass": False}
    if state.get("current_case"):
        directory = args.batch / state["current_case"]
        if (directory / "状态.json").is_file():
            child = read(directory / "状态.json")
            result["case_status"] = {k: child[k] for k in ("status", "updated_at", "version",
                                     "fdtd_status", "elapsed_s", "error") if k in child}
            # 日志头可能含许可证上下文，只输出严格以百分比起始的数值进度行。
            for log in directory.glob("*_p0.log"):
                progress = [line.strip() for line in log.read_text(encoding="utf-8", errors="replace").splitlines()
                            if re.match(r"^\s*\d+(?:\.\d+)?%\s", line)]
                if progress:
                    result["last_numeric_progress"] = progress[-1]
    if state["status"] in ("starting", "running", "analyzing") and not result["controller"]["alive_matching_identity"]:
        result["attention"] = "状态声称运行但控制器异常消失；没有自动重启"
    if state["status"] == "failed":
        result["attention"] = state.get("error", "批次失败")
    if cancellation.is_file():
        result['original_batch_status'] = result['batch_status']
        result['batch_status'] = 'cancelled_by_user_no_restart'
        result['user_decision'] = read(cancellation)
        result['attention'] = '用户已取消补算，保留原始状态文件，不重启队列'
    result["completed"] = state["completed"]
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
