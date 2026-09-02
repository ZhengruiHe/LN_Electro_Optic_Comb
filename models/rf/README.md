# 射频电极模型

本目录保存行波电极的 HFSS 建模、求解与传输线参数后处理脚本。

- `hfss_cpw_config.json`：PDK 堆栈、材料参数、普通 CPW 与 T 形分段电极参数、筛选范围和验收指标。
- `build_hfss_cpw.py`：通过 PyAEDT 调用 HFSS 2023 R1，建立 Driven Terminal 模型、端口、混合模端子、网格与扫频设置。
- `postprocess_tline.py`：利用两种长度的未归一化 Touchstone 数据，提取特性阻抗、传播损耗和微波有效折射率。
- `../../scripts/analyze_hfss_two_length.py`：对 HFSS 导出的四端口混合模数据进行两长度 ABCD 本征值去嵌入，避免强反射时直接相减 S21 造成错误。

## 当前10 GHz T形电极候选结构

当前默认采用 T 形周期微结构电极。Kharel 等人的 Optica 论文只用于确定“减轻窄间隙电流拥挤并利用周期电容调节微波速度”的机制；当前尺寸是面向南智 TFLN-on-SiN PDK 的候选起点，不是论文复现：

- 信号主干宽度：43 µm；
- 两侧地电极主干宽度：100 µm；
- 电光调制内间隙：5 µm；
- T 形帽横向宽度：2 µm；
- T形颈横向长度：4 µm；
- T形颈沿传播方向宽度：10 µm；
- T形帽沿传播方向长度：45 µm；
- 单元间隔：5 µm；
- 周期：50 µm。

T形帽占空比固定为90%，没有通过降低占空比实现速度匹配。当前候选通过缩短横向颈长并加宽沿传播方向的颈部来降低串联电感。严格网格在9.96 GHz得到微波有效折射率2.20009、阻抗48.12 Ω、筛选损耗1.66 dB/cm；与70°、1.33 µm有源波导的TE₀/TE₁平均群折射率2.22095相差-0.02086。当前有源调制区已按用户确认改为完全移除SiN，并在原300 nm SiN高度内全部填入SiO₂；旧25 GHz结果及68%占空比10 GHz种子均保留为历史对照。

同一截面的Maxwell 2D单位线电压场已完成：波导中心处T形帽覆盖区横向场约7.85×10⁴ V/m/V，主干空隙区约2.28×10⁴ V/m/V；与TE0/TE1复数光场重叠后，周期加权半波电压长度积分别为8.15和9.36 V·cm。电光系数采用Ansys官方示例的r33=30.9 pm/V、r13=9.6 pm/V，仍需代工或实测确认。

旧60°侧壁光学模型下，TE0群折射率为 `2.2252`，TE1为 `2.3011`；旧 `nRF` 在25.07 GHz扫频点处位于两者之间。这些数据全部保留，但不作为10 GHz定版依据。25.07 GHz只是旧线性扫频中离25.00 GHz最近的网格点。

当前10 GHz建模配置显式包含高阻硅、8.2 µm底氧、原SiN层位置的300 nm SiO₂填充、400 nm层间氧化层以及400/200/0 nm两级LN图形。LN1和LN2侧壁均按相对水平面70°建模。有源区下方没有残余SiN；射频材料损耗、金属工艺电导率与粗糙度及背面测试边界仍待确认，所以设计状态仍是“PDK兼容候选”，不是“可投片定版”。

参考文献：P. Kharel 等，《利用微结构电极突破集成铌酸锂调制器的电压—带宽限制》，Optica 8, 357–363 (2021)，[论文预印本](https://arxiv.org/abs/2011.13422)，[DOI](https://doi.org/10.1364/OPTICA.416155)。

## 已验证的建模命令

仅检查环境和配置，不启动 HFSS：

```powershell
.\.venv\Scripts\python.exe models\rf\build_hfss_cpw.py check
.\.venv\Scripts\python.exe models\rf\postprocess_tline.py --self-test
```

建立 500 µm 的 T 形电极短线模型：

```powershell
.\.venv\Scripts\python.exe models\rf\build_hfss_cpw.py build --allow-placeholders --non-graphical --aedt-version 2023.1 --electrode-style segmented_t --length-um 500
```

`--allow-placeholders` 只允许检查建模链路，不能把结果当作最终器件数据。当前 PDK 中仍未确认 LN1/LN2 刻蚀深度、射频材料损耗、金属工艺电导率与粗糙度、背面载台边界和探针焊盘，因此正式扫频前必须补齐或做敏感性分析。

生成的 AEDT 工程、网格、Touchstone 数据、场文件和预览图统一保存在 `results/hfss/`，并由 Git 忽略。
