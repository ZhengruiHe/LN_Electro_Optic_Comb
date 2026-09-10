# 版图脚本导航

当前采用 v13_22 固定窗口版本，具体路径、哈希及交付缺口见 [阶段交付核查](../../docs/阶段交付核查与进展_20260911.md)。不是目录名数字最大的文件就自动成为采用版。

## 当前生成与核查流程

1. `route_horizontal_access.py`：10 GHz整数/半整数/整数周期条件下的横向接入及独立时延预算，默认3/3.5/4周期。
2. `build_15mm_ysj_merge.py`、`klayout_build_upper_15mm.py`：从本地冻结源构建四程部分，保留YSJ原结构记录合并到单个block。
3. `finalize_cycle_layout.py`：针对新生成目录回读、间距、时延和预览。
4. `recheck_release_layout.py`：对现有采用版仅写新检查报告，不改GDS，不启动仿真。
5. `audit_merged_block.py`、`audit_route_spacing_gds.py`、`audit_route_proximity.py`：几何保留、同层冲突、GSG连通和近邻范围核对。

运行例子：

```powershell
.\.venv\Scripts\python.exe -X utf8 models/layout/recheck_release_layout.py --directory results/layout/YSJ合并15mm_横向接入_v13_22_PDK固定窗口 --output-dir results/layout/自己的只读核查
.\.venv\Scripts\python.exe -m unittest discover -s models/layout -p "test_*.py"
```

必须提供已有生成目录中的元数据和所引用的授权源文件。Git不包括GDS、PDK和结果，因此空仓库克隆后不能凭空产生相同黑盒；先按中文复现说明准备本地输入。

## 仿真几何与预览

`klayout_extract_optical_roi.py`、`klayout_extract_mux_validation.py`直接读取实际芯层，用于二维几何一致性核对。`preview_*`、`explain_mux_location.py`和`inspect_*`是读取或可视化辅助，不代表光学验证。

## 历史过程脚本

`route_compact_*`、`route_five_ports_global.py`、`route_cycle_matched.py`、`replan_mux_connections.py`、早期GSG/窗口/尾端修正等记录了先前方案；其中一部分仍提供当前脚本使用的几何函数。保留这些依赖，不要为“清理”随意移动Python模块导致导入失效。旧方案的限制和否决原因以当前阶段说明为准。

本地 `probe_*.py` 为一次性中间调试，不作为基准发布。脚本里的扫描/可调参数保留；生成的中间版图只做可恢复归档。

## 规则口径

`.lydrc` 文件及原生检查只是项目已实现的规则子集，不是代工官方全套检查。当前保留40条核心线宽、40条刻蚀间距标记。只读核查通过不等于这些标记已豁免，更不等于生产签核。
