# 稻季控释配方地图

全球 50 个稻作区、73 个稻季的“包膜尿素 + 尿素”一次性基施掺混配方模型。模型的作用：

1. 用当地气温，把 25 °C 静水控释期换算成稻田中的实际释放曲线；
2. 与该季水稻的吸氮进程逐日比对；
3. 求出“全季不断顿、总氮量最少”的掺混方案。

最后生成一张可交互的地图页面 `index.html`，同时评估各方案的环境效应（氨挥发、N₂O、淋溶径流、碳足迹）。

> 所有参数和数据都是文献初值，尚未用您的产品实测数据和田间试验标定。

## 目录结构

```
rice-crf-map/
├── build.py            主程序：计算全部稻季，生成 data.json、index.html、BUILD_INFO.json
├── params.toml         参数文件：释放模型、氮库、需氮模型、品种类型、优化与情景设置
├── crfmap/             模型代码
│   ├── params.py         读取 params.toml
│   ├── data.py           读取 data/ 下的稻作区数据
│   ├── climate.py        月均温 → 逐日气温，日期换算
│   ├── release.py        温度—释放模型（Arrhenius + Weibull）
│   ├── demand.py         水稻需氮模型（积温驱动的 S 形吸氮曲线）
│   └── optimize.py       肥料氮库质量平衡、线性规划、配方取整、分次施尿素对照
├── data/
│   ├── zones.csv         50 个稻作区：代表站、经纬度、一年几季、1–12 月平均气温（°C）
│   ├── seasons.csv       73 个稻季：种植方式、品种类型、施肥日、生育天数、目标产量、土壤供氮、当地常规施氮量
│   ├── countries.csv     ISO 国家代码 → 中文国名
│   └── countries-50m.json  世界地图边界（TopoJSON）
├── template.html       页面模板（样式、地图、图表的前端代码）
├── tests/test_model.py 基本检查与回归测试
├── requirements.txt    Python 依赖
├── data.json           生成结果：全部稻季的计算结果（页面里内嵌同一份）
├── index.html          生成结果：交互地图页面
└── BUILD_INFO.json     生成记录：源代码 commit、运行环境、输出文件的 sha256
```

环境效应和 S 型情景的计算代码在仓库根目录的 `../common/`（与小麦版共用）：
- `env.py`：计算代码和全部环境参数；
- `env_section.html`、`env.js`：页面里的环境效应部分；
- `scenarios.py`、`scenario.js`：S 型产品情景的计算和页面切换。

`build.py` 会自动读取它们，所以请保持仓库的目录结构不变。

## 依赖

- Python ≥ 3.11（用标准库 `tomllib` 读取参数文件）
- numpy、scipy（见 `requirements.txt`，已验证版本为 numpy 2.4.6、scipy 1.17.1）

查看页面不需要安装任何东西。页面在浏览器中打开时，会从 CDN 加载 d3 7.9.0、topojson-client 3.1.0 和 Google 字体，所以需要联网。

## 运行

```bash
cd rice-crf-map
python -m venv .venv && source .venv/bin/activate      # Windows：.venv\Scripts\activate
pip install -r requirements.txt

python -m unittest discover -s tests    # 约 2 秒
python build.py                         # 4 核约 1.5–2 分钟；--jobs 1 为单进程
```

运行后会更新以下三个文件，并在终端打印全部稻季的配方汇总表：
- `data.json`
- `index.html`：用浏览器直接打开即可；
- `BUILD_INFO.json`。

## 产品方案情景

页面顶部的“产品方案”可以在三种方案间切换，地图提示、详情曲线、配方表和环境效应都会随之更新：

| 方案 | 含义 |
|---|---|
| 线性型（现有产品） | 默认方案：β=1.3，尿素 + 1–2 个控释期（CR 档） |
| S 型单档 | 假设的 S 型包膜尿素（β=`params.toml` 中的 `stype_beta`，默认 2.5，有滞后期）：尿素 + 1 个 S 型档位 |
| S 型两档 | 尿素 + 至多 2 个 S 型档位，取整后若只剩一档，则与单档相同 |

- S 型档位的控释期同样按 25 °C 静水 80% 释放天数标定。
- S 型的滞后期对温度更敏感，所以 S 型方案都要求同一配方在常年、偏冷 1.5 °C、偏暖 1.5 °C 三种年份都不断顿，靠 `crfmap.optimize.optimize(robust_dT=...)` 实现。线性型方案本来就满足这个要求。
- S 型方案的计算代码在 `../common/scenarios.py`（与小麦版共用），前端切换在 `../common/scenario.js`。
- 在链接末尾加 `#s1` 或 `#s2`，打开页面时就直接显示对应方案。

## 修改数据和参数

- **气温、种植日历、产量、土壤供氮**：直接编辑 `data/zones.csv`、`data/seasons.csv`。可以用 Excel 打开，保存时选“CSV UTF-8”。新增稻作区时，`id` 不能重复；`seasons.csv` 的 `zone_id` 必须对应 `zones.csv` 中的 `id`。同时要在 `../common/env.py` 的 `RICE_ZONE` 里补上该区的土壤 pH、气候、化肥产地和 CH₄ 区域。
- **模型参数**：编辑 `params.toml`，每一项都有注释。浮点数和整数的写法请保持原样，比如 `17.0` 和 `7`。
- **环境效应参数**：编辑 `../common/env.py` 顶部，包括减排比例、排放因子、损害成本和碳价。
- 也可以用环境变量指定另一份参数文件做情景对比：`CRFMAP_PARAMS=/path/to/other.toml python build.py`。

改完后重新运行 `python build.py`。

### 需要实测的数据（按重要性排序）

1. 每个控释期档位在 15、20、25、30、35 °C 下的静水释放曲线，用来拟合 `Ea`、`beta`、`f0`；
2. 稻田埋网袋的实际释放数据，用来标定 `f_soil`；
3. 各区无氮区的地上部吸氮量，填入 `seasons.csv` 的 `INS` 列；
4. ¹⁵N 示踪或氮素表观平衡试验，用来标定 `k_loss`、`eta`、`urea_loss`；
5. 各站近 20 年的逐日气象数据（ERA5 或 NASA POWER），替换月均温插值。

## 版本与可复现性

每次运行 `build.py` 都会在 `BUILD_INFO.json` 中记录：
- 源代码 commit，以及有没有未提交的改动；
- Python、numpy、scipy 的版本；
- 输出文件的 sha256。

要核对某份 HTML 是哪个版本生成的，比较它的 sha256 即可。

当前线上发布的页面对应 `index.html`，sha256 为：
```
7f0998644ad11e30db6dbbf2ee16c41d23149c4d18f691ee8b576f41112eb673
```
这个文件由 commit `618c51c`（“Add environmental effect assessment to rice and wheat maps”）的源代码生成。之后把代码整理成现在的包结构，只重组了代码、没有改动任何计算。用整理后的代码重新生成，得到的 `index.html` 和 `data.json` 与原来逐字节一致。

加入 S 型产品情景后，页面和 `data.json` 都多了可选方案的数据，所以 sha256 变了；默认线性型方案的每一个数值都与之前相同。最新的对应关系以 `BUILD_INFO.json` 为准。

## 局限

- 气温用的是各代表站的月均温近似值；种植日历、产量、土壤供氮、当地常规施氮量取自文献和公开统计的代表值。同一分区内部，品种和土壤的差别很大。
- 旱直播（美国中南部、巴西、澳洲）播后约有 30 天是旱地，模型只用较高的氮损失率粗略处理。
- 只计算氮；磷、钾请按土壤测试结果一次基施。

## 第三方数据与库

- `data/countries-50m.json` 来自 [world-atlas](https://github.com/topojson/world-atlas) 2.0.2（ISC 许可），底图数据源自 Natural Earth（公有领域）。
- 页面运行时加载的 d3 7.9.0、topojson-client 3.1.0 均为 ISC 许可。
