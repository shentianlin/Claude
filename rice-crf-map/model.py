"""控释尿素温度释放模型 + 水稻需氮模型 + 掺混配方优化（氮库质量平衡 + 线性规划）。

1. 释放模型（文献依据见 README）
   - “控释期”D = 25 °C 静水中累计释放 80% 的天数（与 GB/T 23348、HG/T 4215 定义一致）。
   - 温度只改变释放速率、不改变曲线形状（Hara 2000）。用 Arrhenius 方程把逐日温度折算成
     “25 °C 等效天数” τ：
         a(T) = exp[ Ea/R · (1/298.15 − 1/(T+273.15)) ]，   τ(t) = f_soil · Σ a(T_d)
     Ea 默认 50 kJ/mol（20–30 °C 区间 Q10≈2，Gandeza et al. 1991；Du et al. 2006 实测
     37–46 kJ/mol；Liu 2011 报道 Q10 1.4–2.8）。
   - 曲线形状：带初期溶出的 Weibull 函数
         F(τ) = f0 + (1−f0)·[1 − exp(−(τ/λ)^β)]，由 F(D)=0.80 解出 λ。
2. 需氮模型
   - 作物累计吸氮比例是“有效积温进度” x = GDD(t)/GDD(成熟) 的 S 形曲线，
     用幼穗分化期、抽穗期两点的累计吸氮比例标定。
   - 作物总吸氮 = 目标产量 × 单位籽粒需氮量；其中土壤基础供氮 INS 按同一曲线分摊，
     其余部分必须由肥料供给。
3. 肥料氮库（每天）
       库(t+1) = 库(t)·(1 − k_loss) + 当日释放 − 当日肥料氮吸收 / η
     k_loss：淹水稻田有效氮库的日损失率（氨挥发、硝化—反硝化、淋溶、固定）。
     η：根系对库中氮的最大截获效率。尿素施入后另有一次性的初始损失 L_urea。
   约束“管饱”：库(t) ≥ 未来 buffer 天的吸氮需求（任何一天都不断顿）。
   目标“不浪费”：总施氮量最小（成熟时仍未释放的氮自然计为浪费）。
   在 30–360 天（10 天一档）的系列中枚举“尿素 + 至多 2 个控释期”，每个组合解一个线性规划。
"""
import math
import numpy as np
from scipy.optimize import linprog, minimize

R_GAS = 8.314
PRODUCTS = list(range(30, 361, 10))  # 控释期（天）

P = dict(
    Ea=50.0,          # kJ/mol，Arrhenius 表观活化能
    f_soil=0.90,      # 淹水土壤中相对静水的释放速率系数
    f0=0.03,          # 初期溶出比例（24 h）
    beta=1.3,         # Weibull 形状参数：>1 表示有短暂滞后、中段近线性
    paddy_dT=1.0,     # 田面水/表层土温 = 气温 + 1 °C
    urea_tau=2.0,     # 尿素水解时间常数（天）
    urea_loss=0.15,   # 基施尿素的初始损失（主要是氨挥发，混施入土的情形）
    k_loss=0.013,     # 有效氮库日损失率
    eta=0.80,         # 根系最大截获效率
    buffer=7,         # 库中至少保有未来 7 天的需求
    early_pool=0.20,  # 分蘖期（施肥后到幼穗分化）库中至少保有“肥料需供吸氮量”的 20%，保证分蘖所需浓度
    n_cru=0.42,       # 包膜尿素含氮量
    n_urea=0.46,
)

CROP = {
    "IND": dict(nreq=17.0, pPI=0.40, pH=0.85, label="常规籼稻"),
    "HYB": dict(nreq=17.0, pPI=0.35, pH=0.78, label="杂交籼稻"),
    "JAP": dict(nreq=19.0, pPI=0.40, pH=0.82, label="温带粳稻"),
    "TJ": dict(nreq=17.0, pPI=0.38, pH=0.80, label="热带粳稻/长粒型"),
    "AROM": dict(nreq=18.0, pPI=0.35, pH=0.80, label="香稻/感光品种"),
}
EST = {
    "TPR": dict(label="移栽", dPI=0.0, dH=0.0, k_mult=1.0),
    "WSR": dict(label="水直播/湿直播", dPI=-0.08, dH=-0.03, k_mult=1.0),
    "DSR": dict(label="旱直播（晚淹水）", dPI=-0.08, dH=-0.03, k_mult=1.15),
}
MIN_FERT_DEMAND = 12.0  # kg N/ha：肥料需供吸氮量的下限


# ---------------- 气温 ----------------
def daily_temps(T12, start_doy, n):
    mid = np.array([15.2 + 30.44 * m for m in range(12)])
    xs = np.concatenate([mid - 365, mid, mid + 365])
    ys = np.concatenate([T12, T12, T12])
    d = (start_doy + np.arange(n)) % 365
    return np.interp(d, xs, ys)


def doy(mmdd):
    m, d = map(int, mmdd.split("-"))
    return [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334][m - 1] + d - 1


# ---------------- 释放模型 ----------------
def arrhenius(T, Ea):
    return np.exp(Ea * 1000 / R_GAS * (1 / 298.15 - 1 / (np.asarray(T) + 273.15)))


def lam(D, f0, beta):
    return D / (math.log((1 - f0) / 0.2)) ** (1 / beta)


def release(tau, D, f0=None, beta=None):
    f0 = P["f0"] if f0 is None else f0
    beta = P["beta"] if beta is None else beta
    tau = np.maximum(tau, 0)
    return f0 + (1 - f0) * (1 - np.exp(-(tau / lam(D, f0, beta)) ** beta))


def cum_tau(Tpaddy, Ea, f_soil):
    """第 t 天结束时累计的 25 °C 等效天数。"""
    return np.cumsum(arrhenius(Tpaddy, Ea) * f_soil)


# ---------------- 需氮模型 ----------------
def logistic_norm(x, x0, k):
    L = lambda z: 1 / (1 + np.exp(-k * (z - x0)))
    return (L(x) - L(0)) / (L(1) - L(0))


def fit_logistic(xPI, pPI, xH, pH):
    def err(p):
        return (logistic_norm(xPI, *p) - pPI) ** 2 + (logistic_norm(xH, *p) - pH) ** 2
    best = None
    for x0 in (0.3, 0.45, 0.6):
        for k in (5, 9, 14):
            r = minimize(err, [x0, k], method="L-BFGS-B", bounds=[(0.05, 1.0), (2.0, 30.0)])
            if best is None or r.fun < best.fun:
                best = r
    return float(best.x[0]), float(best.x[1])


def season_setup(zone, s, prm=P):
    L = s["days"]
    n = L + 1                                   # 第 0 … L 天
    start = doy(s["app"])
    Tair = daily_temps(zone["T"], start, n)
    gdd = np.concatenate([[0], np.cumsum(np.clip(np.minimum(Tair, 30) - 10, 0, None))])[:n]
    x = gdd / gdd[-1]
    Tm_end = Tair[-35:].mean()
    gfill = int(round(np.clip(30 + (24 - Tm_end) * 2.0, 25, 45)))
    tH = L - gfill
    tPI = max(tH - 28, 15)
    c, e = CROP[s["type"]], EST[s["est"]]
    pPI, pH = c["pPI"] + e["dPI"], c["pH"] + e["dH"]
    x0, k = fit_logistic(x[tPI], pPI, x[tH], pH)
    Fcrop = logistic_norm(x, x0, k)
    uptake_total = s["Y"] * c["nreq"]
    fert_demand = max(MIN_FERT_DEMAND, uptake_total - s["INS"])   # 需由肥料提供的吸氮量
    d = np.diff(Fcrop, prepend=0.0) * fert_demand                  # 每日肥料氮吸收需求
    return dict(L=L, n=n, start=start, Tair=Tair, Tpaddy=Tair + prm["paddy_dT"], gdd=gdd,
                tPI=tPI, tH=tH, gfill=gfill, Fcrop=Fcrop, d=d, uptake_total=uptake_total,
                fert_demand=fert_demand, x0=x0, k=k, pPI=pPI, pH=pH,
                k_loss=prm["k_loss"] * e["k_mult"])


# ---------------- 氮库线性模型 ----------------
def decay_conv(inc, k):
    """pool[t] = Σ_{s≤t} inc[s]·(1−k)^(t−s)"""
    out = np.empty_like(inc)
    acc = 0.0
    q = 1 - k
    for i, v in enumerate(inc):
        acc = acc * q + v
        out[i] = acc
    return out


def urea_increment(n, prm=P):
    t = np.arange(n + 1)
    cum = 1 - np.exp(-t / prm["urea_tau"])
    return np.diff(cum) * (1 - prm["urea_loss"])


def cr_increment(tau, D):
    cum = np.concatenate([[0.0], release(tau, D)])
    return np.diff(cum)


def season_basis(ss, Ea=None, f_soil=None, dT=0.0, prm=P):
    """返回每种产品 1 kg N 在第 t 天末对氮库的贡献 G_k(t)，以及需求侧 Dreq(t)。"""
    Ea = prm["Ea"] if Ea is None else Ea
    f_soil = prm["f_soil"] if f_soil is None else f_soil
    k = ss["k_loss"]
    tau = cum_tau(ss["Tpaddy"] + dT, Ea, f_soil)
    G = {"urea": decay_conv(urea_increment(ss["n"], prm), k)}
    Fend = {}
    for D in PRODUCTS:
        G[D] = decay_conv(cr_increment(tau, D), k)
        Fend[D] = float(release(tau[-1:], D)[0])
    draw = ss["d"] / prm["eta"]
    Dcum = decay_conv(draw, k)
    buf = buffer_need(ss, prm)
    return G, Dcum + buf, Fend, tau


def buffer_need(ss, prm=P):
    buf = np.array([ss["d"][t + 1:t + 1 + prm["buffer"]].sum() for t in range(ss["n"])])
    early = np.zeros(ss["n"])
    early[1:ss["tPI"]] = prm["early_pool"] * ss["fert_demand"]
    return np.maximum(buf, early) / prm["eta"]


def solve_lp(G, need, keys, cost):
    A = -np.column_stack([G[k] for k in keys])
    c = np.array([cost[k] for k in keys])
    r = linprog(c, A_ub=A, b_ub=-need, bounds=[(0, None)] * len(keys), method="highs")
    if r.status != 0:
        return None
    return r.x


def optimize(ss, Ea=None, f_soil=None, max_cr=2, prm=P, cost_cr=1.02):
    G, need, Fend, tau = season_basis(ss, Ea, f_soil, prm=prm)
    L = ss["L"]
    cands = [D for D in PRODUCTS if Fend[D] >= 0.30]
    cost = {"urea": 1.0, **{D: cost_cr for D in PRODUCTS}}
    res = []
    for i, D1 in enumerate(cands):
        x = solve_lp(G, need, ["urea", D1], cost)
        if x is not None:
            res.append((float(x @ [cost["urea"], cost[D1]]), ("urea", D1), x))
    best1 = min(res, key=lambda r: r[0])
    if max_cr >= 2:
        for i, D1 in enumerate(cands):
            for D2 in cands[i + 1:]:
                x = solve_lp(G, need, ["urea", D1, D2], cost)
                if x is not None:
                    res.append((float(x @ [cost[k] for k in ("urea", D1, D2)]), ("urea", D1, D2), x))
    best = min(res, key=lambda r: r[0])
    # 理论下限：34 个控释期全部可用
    xall = solve_lp(G, need, ["urea"] + cands, cost)
    ideal = float(xall.sum()) if xall is not None else None
    return best, best1, ideal, G, need, Fend, tau


def round_recipe(keys, x, G, need, step=5):
    """比例取 5% 整数倍、去掉 <8% 的组分；在 LP 最优解附近枚举取整方案，选总量最小者（总量取整到 5 kg N）。"""
    import itertools
    tot = x.sum()
    parts = [(k, v / tot) for k, v in zip(keys, x) if v / tot >= 0.08]
    ks = [k for k, _ in parts]
    base = [round(100 * f / step) * step for _, f in parts]
    best = None
    for delta in itertools.product(range(-2 * step, 2 * step + 1, step), repeat=len(ks) - 1):
        pct = [b + d for b, d in zip(base[:-1], delta)]
        last = 100 - sum(pct)
        pct.append(last)
        if min(pct) < 0 or ("urea" in ks and pct[ks.index("urea")] < 10 and len(ks) > 1 and base[ks.index("urea")] >= 10):
            continue
        shape = sum(p / 100 * G[k] for k, p in zip(ks, pct))
        total = float(np.max(need / np.maximum(shape, 1e-9)))
        if best is None or total < best[1]:
            best = ([(k, p) for k, p in zip(ks, pct) if p > 0], total)
    return best[0], 5 * math.ceil(best[1] / 5)


def simulate(ss, recipe, total, G, need, Fend, prm=P):
    pool_supply = sum(total * p / 100 * G[k] for k, p in recipe)
    margin = pool_supply - need
    Dsum = ss["fert_demand"]
    unrel = sum(total * p / 100 * (1 - Fend[k]) for k, p in recipe if k != "urea")
    released = total - unrel
    return dict(
        pool=pool_supply - (need - buffer_need(ss, prm)),
        margin=margin,
        short_days=int(np.sum(margin < -0.5)),
        RE=float(Dsum / total),
        unreleased_pct=float(100 * unrel / total),
        released=float(released),
    )
