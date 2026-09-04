# 系统模型

此处存放基频与二次谐波频率梳动力学的耦合模方程和参数文件。

- `simulate_four_pass_comb.py`：合成MODE光学群折射率、HFSS射频传播参数、Maxwell 2D/MODE光电重叠结果、2T/2.5T/2T回路延迟和理想贝塞尔梳谱，输出四程频率响应、半波电压工程预测和带门限定义的梳齿数。

当前主设计频率已改为10 GHz。四程回路当前采用2T/2.5T/2T，中间的0.5T用于补偿GSG上下调制区的电场极性反转。第三回路由1T改为2T，以增加可实现的被动回环长度。

在70°光学有源波导群折射率和10 GHz HFSS双长度数据生成后，运行当前15 mm T形电极候选：

```powershell
.\.venv\Scripts\python.exe models\system\simulate_four_pass_comb.py `
  --electrode-length-cm 1.5 `
  --eo-overlap-json results\eo\电光重叠_半波电压摘要.json `
  --output-dir results\system\扫描结果\电极对照\15mm\T形电极
```

普通CPW对照必须同时换用普通CPW的HFSS传播数据和Maxwell/MODE重叠结果：

```powershell
.\.venv\Scripts\python.exe models\system\simulate_four_pass_comb.py `
  --rf-csv results\hfss\pdk_10GHz_a70_regular_w43_g5_phasecheck_500_1000.csv `
  --electrode-length-cm 1.5 `
  --eo-overlap-json results\eo\regular_cpw\普通CPW_电光重叠_半波电压摘要.json `
  --output-dir results\system\扫描结果\电极对照\15mm\普通CPW
```

当前光电重叠积分得到T形电极TE0/TE1单程半波电压长度积8.15/9.36 V·cm，普通CPW为7.38/8.46 V·cm。普通CPW的静电场重叠更强，但10 GHz微波有效折射率只有约1.825；T形加载将其提高到约2.200，更接近光学平均群折射率2.221。15 mm四程工程模型中，T形电极与普通CPW的有效半波电压分别约1.671 V和1.684 V，28 dBm理想纯相位模型均给出32根高于最强梳齿−20 dB的谱线。也就是说，当前频点和长度下T形电极的主要价值是速度匹配，不是显著增强局部电场；更严格的优势应由频带、长度和损耗联合比较判断。

所有电极长度扫描和结构对照结果只放在`results/system/扫描结果/`，由Git忽略。基准脚本和本说明上传仓库，避免把大型工程和扫描数据提交到Git。

上述梳齿数尚未计入模式复用器、回环、交叉和各程幅度不均衡造成的光学损耗。

已有的60 µm信号主干、5 µm内间隙和34/16 µm加载/空隙尺寸来自25 GHz历史筛选，现仅作为10 GHz重新扫参的初值。原25 GHz工程、CSV、报告和梳谱预测全部保留，不删除、不覆盖；新10 GHz结果必须使用带`10GHz`的独立文件名。
