# 光学模型

此处存放本项目的LN波导模式、群折射率、模式复用器、交叉器、欧拉弯及二维近邻验证脚本。

2026-09-11阶段决定：用户接受现有二维仿真，不再补跑。四组近邻共8例有结果，其中6例触发能量停止、2例时间窗初筛；保留原始状态，不宣称全部数值收敛。当前交付口径以 [阶段交付核查](../../docs/阶段交付核查与进展_20260911.md) 为准；下文保留早期有源和MUX基线说明，不表示此刻要执行新的容差任务。

当前只读入口为 `link_2d_neighbors_status.py`、`audit_mux_external_results.py`。二维复现入口是 `link_2d_effective_model.py`、`prepare_link_2d_layout_validation.py`、`link_2d_component.py`、`run_link_2d_neighbors.py`，均需用户主动执行才启动计算。`scan_crossing_ln2.py`等参数扫描能力保留，结果不上传。

本地历史脚本 `mux_external_taper_eme.py` 已查出几何和模式标签问题，不纳入当前基准提交，也不将其结果用于版图通过判断。MUX本体的原始`mode_mux_eme.py`基线仍保留。

- `mode_pdk_config.json`：南智 TFLN-on-SiN 堆栈、待扫描波导尺寸和 MODE 求解设置。
- `mode_pdk_sweep.py`：调用本机 Lumerical MODE 2023 R2，扫描1550 nm下的TE0/TE1有效折射率和群折射率。当前采用LN1刻蚀200 nm并保留200 nm薄膜、LN2继续刻透的两级截面；LN1和LN2两次刻蚀形成的侧壁均按相对水平面70°建模，版图宽度解释为梯形顶宽。有源区下方的SiN全部移除，原300 nm SiN高度以SiO₂填满。
- `analyze_mode_localization.py`：从模式场中计算中心能量占比、横向展宽和奇偶对称性，用于区分波导芯区的TE₀/TE₁与LN平台模式。
- `mode_mux_cross_section_sweep.py`：扫描主波导TE₁与辅助波导TE₀的局部避免交叉。
- `mode_mux_eme.py`：建立两级梯形LN结构的三维EME模式复用器，LN1和LN2均以配置中的“与水平面夹角”生成。
- `../../docs/模式复用器EME仿真思路与建立方法.md`：从避免交叉、各向异性坐标、EME单元与端口设置到S矩阵和传播场判读的完整中文学习说明。
- `repropagate_eme_ports.py`：制造偏差导致端口本征模编号重排时，复用工程几何重新选择物理端口模式、重算EME基底并导出S矩阵；MODE 2023 R2切回布局态会清除旧结果，因此不能只做传播。
- `export_mode_mux_fields.py`：导出输入端、避免交叉中心和输出端的TE0直通分支及TE1→TE0转换分支光强图。
- `export_active_mode_fields.py`：导出当前1.33 µm、70°有源LN脊波导的TE0/TE1光强图。

## 输出目录约定

- `results/optical/`根目录：只保留当前名义结构的可复现工程、几何轨迹、S矩阵和阶段报告；
- `results/optical/扫描结果/波导截面/`：保存宽度、刻蚀深度、波长等FDE扫描输出；
- `results/optical/扫描结果/模式复用器/`：保存局部避免交叉、长度、间隙、波长和工艺角点扫描输出；
- 扫描工程和数值结果由Git忽略，不上传仓库；仓库只跟踪建模脚本、配置和中文使用说明。

运行扫描时若显式指定`--output`、`--project`、`--geometry-output`或`--s-output`，也应遵守上述目录约定，避免覆盖基准工程。

## 当前光学结论

四程有源段不是传统的单模波导，而是有意使用TE₀和TE₁两个芯区模式：接入波导保持TE₀单模，模式复用器在有源段选择性地产生TE₁。用户已把LN1和LN2两次刻蚀侧壁角都更新为相对水平面70°。此前1.28 μm及其群折射率数据来自60°模型，只保留为历史对照，不能再用于设计。

当前70°有源波导名义顶宽选为1.33 µm。1550 nm严格网格结果为TE₀群折射率2.22130、TE₁群折射率2.22061，两者差0.00069，射频慢波目标取平均值2.22095。名义模式复用器、有源直波导以及与Maxwell 2D射频场的Pockels重叠积分已经完成；后续仍需完成宽度/刻蚀/间隙联合容差、弯曲与整链路级联。

模式复用器已完成一个较危险的组合角点：主波导和辅助波导同时增宽50 nm，耦合间隙同时缩小100 nm。端口光场识别表明输出端物理模式编号由名义值重排为`[1,2,4]`，校正端口后，TE₀直通率为99.9886%，主路TE₁转辅助TE₀为95.6672%（0.1924 dB），另一独立转换通道为99.6724%。这说明机制仍成立，但TE₁转换已对组合偏差表现出敏感性；反向线宽/间隙偏差以及LN刻蚀深度、膜厚角点尚未完成。

详细结果见 `results/optical/光学双模候选报告.md`和 `results/optical/70度模式复用器阶段报告.md`。

运行完整扫描（结果默认写入`results/optical/扫描结果/波导截面/`）：

```powershell
.\.venv\Scripts\python.exe -u models\optical\mode_pdk_sweep.py
```
