# 光电重叠仿真

本目录用于把光学模式与电极电场合并，而不是重复计算HFSS传输线参数。

当前流程：

1. 用HFSS两长度模型提取传播常数、阻抗和射频损耗；
2. 用Maxwell 2D按同一PDK截面求单位线电压下的横向准静电场；
3. 对T形帽覆盖区和主干空隙区分别求解，再按45/50与5/50的周期占比加权；
4. 与MODE导出的TE0、TE1复电场做Pockels张量重叠积分；
5. 输出单程半波电压长度积，并送入四程系统模型。

当前电极为10 GHz候选：信号主干43 µm、T形帽内间隙5 µm、横向颈长4 µm、
纵向颈宽10 µm、帽长45 µm、周期空隙5 µm、周期50 µm。LN波导顶宽
1.33 µm，LN1/LN2侧壁均按与水平面70°处理，有源区下方SiN完全去除并
用SiO2填充。

电光系数不属于当前PDK已确认参数。后续计算会把所用的文献起始值及其
不确定性单独写入结果，不能把估算的半波电压当作代工保证值。

运行顺序：

```powershell
.\.venv\Scripts\python.exe models\eo\build_maxwell2d_cross_section.py --section both --project results\eo\LN_EO_Comb_EO_PDK_10GHz_2023R1_v3.aedt --non-graphical
.\.venv\Scripts\python.exe models\eo\export_optical_fields_for_overlap.py
.\.venv\Scripts\python.exe models\eo\calculate_eo_overlap.py
```

AEDT 2023.1不能把场直接导出到中文文件名，所以原始场文件使用英文文件名，
之后生成的CSV、图片和中文报告仍使用中文名称。脚本还显式设置LN多边形内部参与
求解，以规避PyAEDT 1.4.0在Maxwell 2D各向异性多边形上的默认属性问题。
