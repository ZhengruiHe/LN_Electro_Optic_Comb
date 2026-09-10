"""联合核对扫描记录、Windows进程和实际输出；只读，不把陈旧running当作仍在计算。"""
from __future__ import annotations
import argparse
import ctypes
from ctypes import wintypes
import json
from pathlib import Path
import time


def process_info(pid):
    if not pid:
        return {"alive":False}
    k=ctypes.WinDLL("kernel32",use_last_error=True)
    k.OpenProcess.argtypes=[wintypes.DWORD,wintypes.BOOL,wintypes.DWORD]
    k.OpenProcess.restype=wintypes.HANDLE
    k.GetExitCodeProcess.argtypes=[wintypes.HANDLE,ctypes.POINTER(wintypes.DWORD)]
    k.CloseHandle.argtypes=[wintypes.HANDLE]
    handle=k.OpenProcess(0x1000,False,int(pid))
    if not handle:
        error=ctypes.get_last_error()
        return {"alive":False if error==87 else None,"win_error":error}
    try:
        code=wintypes.DWORD()
        if not k.GetExitCodeProcess(handle,ctypes.byref(code)):
            return {"alive":None,"win_error":ctypes.get_last_error()}
        return {"alive":code.value==259,"exit_code":code.value}
    finally:
        k.CloseHandle(handle)


def main():
    p=argparse.ArgumentParser()
    p.add_argument("status",type=Path)
    a=p.parse_args()
    s=json.loads(a.status.read_text(encoding="utf-8"))
    process=process_info(s.get("pid"))
    rows=s.get("completed",[])
    files_ok=all((Path(r["directory"])/"结果.json").is_file() for r in rows)
    effective=s["state"]
    if s["state"]=="running" and process["alive"] is False:
        effective="interrupted_or_stale"
    age=time.time()-a.status.stat().st_mtime
    print(json.dumps({"recorded_state":s["state"],"effective_state":effective,"pid":s.get("pid"),
                      "controller_process":process,"child_process":process_info(s.get("child_pid")),
                      "record_age_seconds":round(age,1),"current":s.get("current"),
                      "completed":len(rows),"planned":len(s.get("plan",[])),
                      "completed_result_files_present":files_ok,
                      "confirmed_nominal_pass":s.get("confirmed_nominal_pass",False),
                      "checks":s.get("checks"),"error":s.get("error")},ensure_ascii=False,indent=2))


if __name__=="__main__": main()
