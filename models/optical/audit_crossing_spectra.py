"""只读核查原始端口、模式及场数据的真实波长，另写审计报告，不篡改历史数值。"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np


def audit(directory):
    path=directory/"结果.json"
    data=json.loads(path.read_text(encoding="utf-8"))
    target=float(data["wavelength_nm"])
    def close(values):
        return len(values)==1 and abs(values[0]-target)<1e-6
    ports={}
    for name,result in data["ports"].items():
        actual=np.asarray(result["expansion"]["lambda"],dtype=float).ravel()*1e9
        ports[name]=actual.tolist()
    with np.load(directory/"场分布.npz") as field:
        field_nm=(field["lambda"].ravel()*1e9).tolist()
    modes={}
    for name in ports:
        with np.load(directory/f"端口模式_{name}.npz") as mode:
            modes[name]=(mode["lambda"].ravel()*1e9).tolist()
    passed=all(close(v) for v in ports.values()) and close(field_nm) and all(close(v) for v in modes.values())
    return {"directory":str(directory.resolve()),"target_label_nm":target,
            "port_actual_nm":ports,"field_actual_nm":field_nm,"port_mode_actual_nm":modes,
            "actual_spectral_match":passed,"model_version":data.get("model_version","historical_v1"),
            "result_sha256":hashlib.sha256(path.read_bytes()).hexdigest(),
            "conclusion":"波长一致" if passed else "不可作为标注波长的验证通过证据；只保留原始频点的历史初筛"}


def main():
    p=argparse.ArgumentParser()
    p.add_argument("directory",type=Path)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():
        raise FileExistsError("审计报告已存在，请使用新文件名")
    rows=[audit(path.parent) for path in a.directory.rglob("结果.json")]
    result={"audited":len(rows),"mismatched":sum(not r["actual_spectral_match"] for r in rows),
            "rows":rows,"note":"此审计优先于历史扫描状态中的nominal_crossing_pass；不修改历史原始结果"}
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({"audited":result["audited"],"mismatched":result["mismatched"],"report":str(a.output)},ensure_ascii=False))


if __name__=="__main__": main()
