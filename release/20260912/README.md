# 2026-09-12 自有轻量快照

这是可纳入Git的学习/复核快照，**不是完整芯片或生产签核包**。完整说明见 [Git提交与复现说明](../../docs/Git提交与复现说明_20260912.md)。

## 内容

- self_cells/：从v7主干逐一导出的5个自定义组件，每个GDS只有一个顶层、没有子实例。未复制PDK交叉器、端面耦合器或YSJ单元。
- resistor_readback.json：有效电阻尺寸和名义50 Ω预算；地桥平面跨光芯的原始检查保留。
- mux_local_geometry_check.json：当前M04单元局部线宽/间距标记为0，非全片DRC。
- O_input_GSG150_8_12GHz.s4p、hfss_input_cpw_summary.json：局部RF输入原始数值与物理CPW通道摘要。
- Al_absorption_simplified.json、TiN_absorption_simplified.json：简化光学吸收对照，保留模型限制。
- archive_index.json：99个目录、425个文件的相对路径和哈希索引；不包含归档文件本体。
- snapshot_manifest.json：上述导出文件的大小和SHA256。

## 使用边界

GDS单元必须结合本地PDK工艺映射才能使用，不能直接把独立组件拼在一起就宣称完成全芯片。光学与RF数值只适用于文件中列明的简化几何/材料；当前用户选择不再深入这些初筛问题，不等于官方DRC或数值收敛通过。

导出过程由 scripts/export_git_snapshot.py 复现，需要现有本地v7和历史结果输入以及KLink。原始GDS、结果、PDK文件均不修改。
