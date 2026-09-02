# 系统模型

此处存放基频与二次谐波频率梳动力学的耦合模方程和参数文件。

- `simulate_four_pass_comb.py`：合成MODE光学群折射率、HFSS射频传播参数、Maxwell 2D/MODE光电重叠结果、2T/2.5T/2T回路延迟和理想贝塞尔梳谱，输出四程频率响应、半波电压工程预测和带门限定义的梳齿数。

当前主设计频率已改为10 GHz。四程回路当前采用2T/2.5T/2T，中间的0.5T用于补偿GSG上下调制区的电场极性反转。第三回路由1T改为2T，以增加可实现的被动回环长度。

在70°光学有源波导群折射率和10 GHz HFSS新数据生成后，运行当前10 GHz四程候选：

```powershell
.\.venv\Scripts\python.exe models\system\simulate_four_pass_comb.py --rf-adaptive-converged --rf-port-extra-mode-cleared --output-dir results\system\four_pass_pdk_10GHz
```

当前光电重叠积分得到TE0和TE1的单程半波电压长度积分别为8.15和9.36 V·cm。1 cm四程在完全同相、不计射频损耗时的理想半波电压为2.18 V；接入HFSS沿程损耗与回波后，10 GHz工程预测为2.39 V。28 dBm、50 Ω匹配负载的理想纯相位调制模型给出27根高于最强梳齿−20 dB的谱线。该梳齿数尚未计入模式复用器、回环、交叉和各程幅度不均衡造成的光学损耗。

已有的60 µm信号主干、5 µm内间隙和34/16 µm加载/空隙尺寸来自25 GHz历史筛选，现仅作为10 GHz重新扫参的初值。原25 GHz工程、CSV、报告和梳谱预测全部保留，不删除、不覆盖；新10 GHz结果必须使用带`10GHz`的独立文件名。
