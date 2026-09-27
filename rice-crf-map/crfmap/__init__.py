"""稻季控释配方模型：包膜尿素温度释放模型、水稻需氮模型、肥料氮库与线性规划优化。

模块
  params   读取 params.toml
  data     读取 data/ 下的稻作区气温与作物日历
  climate  月均温插值成逐日气温、日期换算
  release  温度—释放模型（Arrhenius + Weibull）
  demand   水稻需氮模型（积温驱动的 S 形吸氮曲线）
  optimize 肥料氮库质量平衡、线性规划与配方取整
"""
