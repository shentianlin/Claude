"""控释掺混肥的环境效应：氨挥发、氧化亚氮、淋溶径流、温室气体、足迹与社会损害成本。

方法（仿 Gu et al. 2023 Nature 与 IPCC 2019 Tier 1 / LCA 做法）
  1. 基线损失 = 施氮量 × 排放因子（按作物、气候、土壤 pH 取值）。
  2. 控释效应 = 按掺混中包膜尿素所占比例 c，对损失乘以 (1 − c·RR)，
     RR 为 meta 分析得到的“控释尿素相对普通尿素（同氮量）”的减排比例。
  3. 间接 N2O = EF4 × 挥发的 NH3-N + EF5 × 淋溶径流 N。
  4. 温室气体（LCA，摇篮到农场大门，化肥部分）= 化肥生产运输 + 田间 N2O（直接 + 间接）
     [+ 水稻 CH4 变化，单列]。
  5. 足迹：以每吨籽粒为功能单位（假定三种方案产量相同，偏保守）。
  6. 社会损害成本：各形态氮的单位损害成本（USD/kg N）× 损失量；N2O、CH4 与生产排放按碳价折算。
  7. 不确定性：对排放因子、RR、损害成本做 Monte Carlo（默认 2000 次），给出 5–95% 区间。

三种方案
  FP  农户常规：当地常规施氮量，全部普通尿素
  OPT 同标准分次施尿素：本模型算出的“管饱”最少尿素量
  CRF 推荐控释掺混：本模型推荐的总氮量与包膜比例
"""
import numpy as np

GWP_N2O = 273.0      # IPCC AR6，100 年
GWP_CH4 = 27.0       # IPCC AR6，非化石来源 CH4
N2O_N_TO_N2O = 44.0 / 28.0

# ---------------- 控释尿素相对普通尿素的减排比例（同氮量） ----------------
# 取 (中值, 低, 高)，Monte Carlo 用三角分布
RR = {
    "rice": dict(
        NH3=(0.35, 0.25, 0.45),     # Jiang et al. 2021: −35.9%；Lu et al. 2023 稻季: −34.8%
        N2O=(0.27, 0.15, 0.38),     # Jiang 2021: −25.6%；Lu 2023 稻季: −29.1%
        LEACH=(0.15, 0.05, 0.27),   # Jiang 2021: 淋溶 −14.9%
        CH4=(0.18, 0.00, 0.30),     # Jiang 2021: −18.3%（机理证据较弱，单列）
    ),
    "wheat": dict(
        NH3=(0.20, 0.10, 0.35),     # Lu et al. 2023 麦季: −19.8%；Yang et al. 2021 全作物: −24–46%
        N2O=(0.23, 0.12, 0.35),     # Lu 2023 麦季: −22.8%
        LEACH=(0.27, 0.15, 0.40),   # Zhang et al. 2019（旱地玉米）: −27.1%，作为旱地作物代理
        CH4=(0.0, 0.0, 0.0),
    ),
}

# ---------------- 排放因子 ----------------
# 普通尿素的氨挥发比例（中性土壤），土壤 pH 修正
NH3_BASE = {"rice": (0.20, 0.12, 0.30), "wheat": (0.11, 0.06, 0.18)}
PH_MULT = {"acid": 0.7, "neutral": 1.0, "alkaline": 1.4}
# 直接 N2O（IPCC 2019 Refinement 表 11.1）
EF1 = {"flooded": (0.004, 0.001, 0.008), "wet": (0.016, 0.013, 0.019), "dry": (0.005, 0.001, 0.011)}
EF4 = {"wet": (0.014, 0.011, 0.017), "dry": (0.005, 0.000, 0.011)}
EF5 = (0.011, 0.000, 0.020)
# 淋溶 + 径流比例（IPCC FracLEACH-(H)=0.24 为湿润区旱地通用值；稻田与干旱区按实测文献下调）
LEACH_FRAC = {"paddy": (0.10, 0.05, 0.18), "wet": (0.20, 0.10, 0.30), "dry": (0.05, 0.00, 0.10)}

# 化肥生产与运输（kg CO2e / kg N）：中国煤头约 8（Zhang et al. 2013；Luo et al. 2024），全球平均约 4（Menegat et al. 2022）
FERT_PROD = {"china": 8.0, "other": 4.0}
COATING_EXTRA = 0.43   # 包膜聚合物：约 6% 质量、~3 kg CO2e/kg 聚合物 → 每 kg 包膜尿素 N 增加 ~0.43 kg CO2e

# 水稻 CH4 基线（IPCC 2019 表 5.11，持续淹水、无有机肥，kg CH4/ha/天）
CH4_DAILY = {"east_asia": 1.32, "southeast_asia": 1.22, "south_asia": 0.85, "europe": 1.56,
             "north_america": 0.65, "africa": 1.19, "south_america": 1.27, "west_asia": 1.19, "oceania": 1.19}

# 单位社会损害成本（USD/kg N），(中值, 低, 高)
# NH3：PM2.5 健康 + 生态（van Grinsven et al. 2013 欧洲 €2–20 健康 + €2–10 生态；Keeler et al. 2016 显示因地而异跨数量级）
# 淋溶径流：水体富营养化 + 饮用水（van Grinsven 2013 €5–20）
DAMAGE = {"NH3": (15.0, 5.0, 40.0), "LEACH": (10.0, 3.0, 25.0)}
CARBON_PRICE = (50.0, 25.0, 190.0)   # USD/t CO2e；190 为美国 EPA 2023 社会碳成本

# 各分区的土壤 pH 类别、气候类别、化肥产地、CH4 区域
RICE_ZONE = {
    # id: (pH, 气候 wet/dry, 化肥产地, CH4 区域)
    "cn-hlj": ("neutral", "wet", "china", "east_asia"), "cn-jl": ("neutral", "wet", "china", "east_asia"),
    "cn-hh": ("alkaline", "wet", "china", "east_asia"), "cn-lyz": ("neutral", "wet", "china", "east_asia"),
    "cn-myz": ("acid", "wet", "china", "east_asia"), "cn-hn": ("acid", "wet", "china", "east_asia"),
    "cn-hs": ("acid", "wet", "china", "east_asia"), "cn-sc": ("neutral", "wet", "china", "east_asia"),
    "cn-yn": ("acid", "wet", "china", "east_asia"), "jp": ("acid", "wet", "other", "east_asia"),
    "kr": ("acid", "wet", "other", "east_asia"), "tw": ("acid", "wet", "other", "east_asia"),
    "in-pb": ("alkaline", "dry", "other", "south_asia"), "in-bh": ("neutral", "wet", "other", "south_asia"),
    "in-wb": ("neutral", "wet", "other", "south_asia"), "in-cg": ("acid", "wet", "other", "south_asia"),
    "in-ap": ("neutral", "wet", "other", "south_asia"), "in-tn": ("neutral", "wet", "other", "south_asia"),
    "bd": ("neutral", "wet", "other", "south_asia"), "pk-pb": ("alkaline", "dry", "other", "south_asia"),
    "pk-sd": ("alkaline", "dry", "other", "south_asia"), "np": ("neutral", "wet", "other", "south_asia"),
    "lk": ("neutral", "wet", "other", "south_asia"), "mm": ("acid", "wet", "other", "southeast_asia"),
    "th-c": ("acid", "wet", "other", "southeast_asia"), "th-ne": ("acid", "wet", "other", "southeast_asia"),
    "vn-m": ("acid", "wet", "other", "southeast_asia"), "vn-r": ("neutral", "wet", "other", "southeast_asia"),
    "kh": ("acid", "wet", "other", "southeast_asia"), "ph": ("neutral", "wet", "other", "southeast_asia"),
    "id": ("acid", "wet", "other", "southeast_asia"), "my": ("acid", "wet", "other", "southeast_asia"),
    "ir": ("neutral", "wet", "other", "west_asia"), "eg": ("alkaline", "dry", "other", "africa"),
    "sn": ("neutral", "dry", "other", "africa"), "ml": ("acid", "dry", "other", "africa"),
    "ng": ("acid", "dry", "other", "africa"), "tz": ("acid", "wet", "other", "africa"),
    "mg": ("acid", "wet", "other", "africa"), "it": ("neutral", "wet", "other", "europe"),
    "es": ("alkaline", "dry", "other", "europe"), "ru": ("neutral", "dry", "other", "europe"),
    "us-ar": ("neutral", "wet", "other", "north_america"), "us-la": ("acid", "wet", "other", "north_america"),
    "us-ca": ("neutral", "dry", "other", "north_america"), "br": ("acid", "wet", "other", "south_america"),
    "uy": ("neutral", "wet", "other", "south_america"), "co": ("acid", "wet", "other", "south_america"),
    "pe": ("alkaline", "dry", "other", "south_america"), "au": ("neutral", "dry", "other", "oceania"),
}
# 小麦：气候类别按 IPCC 定义，灌溉区计为湿润
WHEAT_ZONE = {
    "cn-hbn": ("alkaline", "wet", "china"), "cn-hhs": ("alkaline", "wet", "china"),
    "cn-yz": ("neutral", "wet", "china"), "cn-sc": ("neutral", "wet", "china"),
    "cn-lp": ("alkaline", "dry", "china"), "cn-xj": ("alkaline", "wet", "china"),
    "cn-nx": ("alkaline", "wet", "china"), "cn-db": ("neutral", "dry", "china"),
    "in-pb": ("alkaline", "wet", "other"), "in-up": ("alkaline", "wet", "other"),
    "in-mp": ("alkaline", "wet", "other"), "pk": ("alkaline", "wet", "other"),
    "ir": ("alkaline", "dry", "other"), "tr": ("alkaline", "dry", "other"),
    "eg": ("alkaline", "wet", "other"), "ma": ("alkaline", "dry", "other"),
    "et": ("acid", "wet", "other"), "fr": ("neutral", "wet", "other"),
    "uk": ("neutral", "wet", "other"), "de": ("neutral", "wet", "other"),
    "it": ("alkaline", "dry", "other"), "ua": ("neutral", "dry", "other"),
    "ru-s": ("neutral", "dry", "other"), "kz": ("alkaline", "dry", "other"),
    "us-ks": ("neutral", "dry", "other"), "us-nd": ("neutral", "dry", "other"),
    "us-pnw": ("neutral", "dry", "other"), "ca": ("neutral", "dry", "other"),
    "mx": ("alkaline", "wet", "other"), "ar": ("neutral", "wet", "other"),
    "br": ("acid", "wet", "other"), "au-wa": ("acid", "dry", "other"),
    "au-nsw": ("acid", "dry", "other"),
}

N_MC = 2000


def _draw(rng, spec, n):
    mid, lo, hi = spec
    if hi <= lo:
        return np.full(n, mid)
    return rng.triangular(lo, mid, hi, n)


def season_env(crop, zone_id, s, rng):
    """返回一个稻季/麦季三种方案的环境指标（中值与 5–95% 区间）。"""
    if crop == "rice":
        ph, clim, prod, ch4reg = RICE_ZONE[zone_id]
        ef1_key, leach_key = "flooded", "paddy"
    else:
        ph, clim, prod = WHEAT_ZONE[zone_id]
        ch4reg = None
        ef1_key, leach_key = clim, clim
    n = N_MC
    rr = {k: _draw(rng, v, n) for k, v in RR[crop].items()}
    ef_nh3 = _draw(rng, NH3_BASE[crop], n) * PH_MULT[ph]
    ef1 = _draw(rng, EF1[ef1_key], n)
    ef4 = _draw(rng, EF4[clim], n)
    ef5 = _draw(rng, EF5, n)
    fl = _draw(rng, LEACH_FRAC[leach_key], n)
    d_nh3 = _draw(rng, DAMAGE["NH3"], n)
    d_leach = _draw(rng, DAMAGE["LEACH"], n)
    cprice = _draw(rng, CARBON_PRICE, n) / 1000.0     # USD/kg CO2e

    c_crf = sum(p for k, p in s["rec"] if k != "urea") / 100.0
    scen = {"FP": (float(s["FN"]), 0.0), "OPT": (float(s["split_total"]), 0.0), "CRF": (float(s["total"]), c_crf)}
    ch4_base = CH4_DAILY[ch4reg] * s["days"] if crop == "rice" else 0.0

    out = {}
    draws = {}
    for key, (N, c) in scen.items():
        nh3 = N * ef_nh3 * (1 - c * rr["NH3"])
        n2o_d = N * ef1 * (1 - c * rr["N2O"])
        leach = N * fl * (1 - c * rr["LEACH"])
        n2o_i = ef4 * nh3 + ef5 * leach
        n2o = n2o_d + n2o_i
        nr = nh3 + n2o + leach
        ghg_n2o = n2o * N2O_N_TO_N2O * GWP_N2O
        ghg_prod = N * FERT_PROD[prod] + N * c * COATING_EXTRA
        ch4 = ch4_base * (1 - c * rr["CH4"])
        ghg = ghg_n2o + ghg_prod
        ghg_all = ghg + ch4 * GWP_CH4
        dmg = nh3 * d_nh3 + leach * d_leach + ghg * cprice
        Y = s["Y"]
        m = dict(N=N, NH3=nh3, N2O_direct=n2o_d, N2O_indirect=n2o_i, N2O=n2o, LEACH=leach, Nr=nr,
                 GHG_N2O=ghg_n2o, GHG_prod=ghg_prod, GHG=ghg, CH4=ch4, GHG_all=ghg_all,
                 NF=nr / Y, CF=ghg / Y, CF_all=ghg_all / Y, DMG=dmg)
        draws[key] = m
        out[key] = {k: _stat(v) for k, v in m.items()}
    # CRF 相对 FP 的变化，并分解为“减量效应”（同为尿素、氮量 FP→CRF）与“控释效应”（CRF 氮量下 c·RR）
    diff = {}
    for k in ("NH3", "N2O", "LEACH", "Nr", "GHG", "GHG_all", "DMG", "CF", "NF", "CH4"):
        a, b = draws["FP"][k], draws["CRF"][k]
        diff[k] = dict(abs=_stat(b - a), pct=_stat(np.where(a > 0, 100 * (b - a) / np.maximum(a, 1e-9), 0.0)))
    # 分解（中值）
    def at(N, c):
        nh3 = N * ef_nh3 * (1 - c * rr["NH3"])
        leach = N * fl * (1 - c * rr["LEACH"])
        n2o = N * ef1 * (1 - c * rr["N2O"]) + ef4 * nh3 + ef5 * leach
        return float(np.median(nh3 + leach + n2o))
    FPn, CRn = scen["FP"][0], scen["CRF"][0]
    rate_eff = at(CRn, 0.0) - at(FPn, 0.0)
    prod_eff = at(CRn, c_crf) - at(CRn, 0.0)
    return dict(ph=ph, climate=clim, prod=prod, ch4_region=ch4reg, c_crf=round(c_crf, 2),
                scen=out, diff=diff, decomp=dict(rate=round(rate_eff, 2), product=round(prod_eff, 2)))


def _stat(x):
    q = np.percentile(x, [5, 50, 95])
    return [round(float(q[0]), 3), round(float(q[1]), 3), round(float(q[2]), 3)]


def add_env(crop, zones, seed=20260927):
    rng = np.random.default_rng(seed)
    for z in zones:
        for s in z["seasons"]:
            s["env"] = season_env(crop, z["id"], s, rng)
    return dict(
        crop=crop, RR=RR[crop], NH3_BASE=NH3_BASE[crop], PH_MULT=PH_MULT, EF1=EF1, EF4=EF4, EF5=EF5,
        LEACH_FRAC=LEACH_FRAC, FERT_PROD=FERT_PROD, COATING_EXTRA=COATING_EXTRA, DAMAGE=DAMAGE,
        CARBON_PRICE=CARBON_PRICE, GWP_N2O=GWP_N2O, GWP_CH4=GWP_CH4, CH4_DAILY=CH4_DAILY, N_MC=N_MC,
    )


def add_env_alt(crop, zones, key, seed=20260927):
    """为可选情景 s["alt"][key] 计算环境效应。

    随机数种子与遍历顺序和 add_env 相同，因此同一季别在各情景下使用同一组 Monte Carlo 抽样
    （公共随机数），情景之间的差别不受抽样噪声影响。
    """
    rng = np.random.default_rng(seed)
    for z in zones:
        for s in z["seasons"]:
            a = s["alt"][key]
            s["alt"][key]["env"] = season_env(crop, z["id"], dict(s, rec=a["rec"], total=a["total"]), rng)
