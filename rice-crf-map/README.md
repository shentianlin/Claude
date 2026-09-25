# 稻季控释配方地图

全球 44 个稻作区、73 个稻季的“包膜尿素 + 尿素”一次性基施掺混配方模型。

- `zones.py`：稻作分区数据库，包括代表站月均温、种植日历、品种类型、目标产量、土壤供氮、当地常规施氮量。
- `model.py`：包膜尿素温度释放模型（Arrhenius + Weibull）、水稻需氮模型（积温驱动的 S 形曲线）、肥料氮库质量平衡与线性规划优化。
- `build.py`：计算全部稻季，输出 `data.json` 和交互地图 `index.html`。
- `template.html`：地图页面模板。

运行：

```bash
pip install numpy scipy
python3 build.py
```

所有参数都是文献初值，需要用实测数据替换：各档产品在 5 个温度下的释放曲线（Ea、β）、田间埋袋数据（f_soil）、无氮区吸氮量（INS）、¹⁵N 试验（k_loss、η）、逐日气象数据。
