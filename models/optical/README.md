# 光学模型

此处存放波导模式、色散、跑道形谐振器谐振、波长选择性耦合器和 PPLN 模场重叠的源模型及可复现脚本。

- `mode_pdk_config.json`：南智 TFLN-on-SiN 堆栈、待扫描波导尺寸和 MODE 求解设置。
- `mode_pdk_sweep.py`：调用本机 Lumerical MODE 2023 R2，扫描1550 nm下的TE0/TE1有效折射率和群折射率。当前采用LN1刻蚀200 nm并保留200 nm薄膜、LN2继续刻透的两级截面；LN1和LN2两次刻蚀形成的侧壁均按相对水平面70°建模，版图宽度解释为梯形顶宽。有源区下方的SiN全部移除，原300 nm SiN高度以SiO₂填满。
- `analyze_mode_localization.py`：从模式场中计算中心能量占比、横向展宽和奇偶对称性，用于区分波导芯区的TE₀/TE₁与LN平台模式。
- `mode_mux_cross_section_sweep.py`：扫描主波导TE₁与辅助波导TE₀的局部避免交叉。
- `mode_mux_eme.py`：建立两级梯形LN结构的三维EME模式复用器，LN1和LN2均以配置中的“与水平面夹角”生成。
- `export_mode_mux_fields.py`：导出输入端、避免交叉中心和输出端的TE0直通分支及TE1→TE0转换分支光强图。
- `export_active_mode_fields.py`：导出当前1.33 µm、70°有源LN脊波导的TE0/TE1光强图。

## 当前光学结论

四程有源段不是传统的单模波导，而是有意使用TE₀和TE₁两个芯区模式：接入波导保持TE₀单模，模式复用器在有源段选择性地产生TE₁。用户已把LN1和LN2两次刻蚀侧壁角都更新为相对水平面70°。此前1.28 μm及其群折射率数据来自60°模型，只保留为历史对照，不能再用于设计。

当前70°有源波导名义顶宽选为1.33 µm。1550 nm严格网格结果为TE₀群折射率2.22130、TE₁群折射率2.22061，两者差0.00069，射频慢波目标取平均值2.22095。名义模式复用器、有源直波导以及与Maxwell 2D射频场的Pockels重叠积分已经完成；后续仍需完成宽度/刻蚀/间隙联合容差、弯曲与整链路级联。

详细结果见 `results/optical/光学双模候选报告.md`和 `results/optical/70度模式复用器阶段报告.md`。

运行完整扫描：

```powershell
.\.venv\Scripts\python.exe -u models\optical\mode_pdk_sweep.py
```
