"""小麦版：控释尿素温度释放模型 + 小麦需氮模型 + 掺混配方优化（氮库质量平衡 + 线性规划）。

与水稻版的差别
  - 旱地土壤：释放系数 f_soil 取 0.80（湿润土壤中慢于静水，Du et al. 2006）；
    土温 = 气温 + 1 °C，气温低于 0 °C 时按积雪/土层缓冲取气温的 40%；
    土温 0–3 °C 线性降到 0：冻土中包膜内水分冻结，释放停止。
  - 氮库日损失率随土温变化：k(T) = k20 · 2^((T−20)/10) · wet，
    冬季低温时损失小，冬春多雨地区 wet>1（淋溶、渍害）。
  - 需氮曲线用积温（基温 0 °C）进度 x 的 S 形曲线，由拔节期、开花期的累计吸氮比例标定，
    生育期直接取当地农时（播种、拔节、开花、收获日期）。
  - “管饱”约束：施肥后到拔节，氮库不低于肥料需供吸氮量的 15%（冬前分蘖、返青）；
    之后不低于未来 7 天的需求。
"""
import math
import numpy as np
from scipy.optimize import linprog, minimize

R_GAS = 8.314
PRODUCTS = list(range(30, 361, 10))  # 控释期（天）

P = dict(
    Ea=50.0,          # kJ/mol，Arrhenius 表观活化能
    f_soil=0.80,      # 湿润旱地土壤中相对静水的释放速率系数
    f0=0.03,          # 初期溶出比例（24 h）
    beta=1.3,         # Weibull 形状参数：>1 表示有短暂滞后、中段近线性
    soil_dT=1.0,      # 5–10 cm 土温 = 气温 + 1 °C
    urea_tau=3.0,     # 尿素水解时间常数（天）
    urea_loss=0.10,   # 基施尿素的初始损失（混施入土、条施的情形）
    k_loss=0.010,     # 20 °C 时有效氮库日损失率（硝态氮淋溶、反硝化、固定、挥发）
    eta=0.80,         # 根系最大截获效率
    buffer=7,         # 库中至少保有未来 7 天的需求
    early_pool=0.15,  # 施肥后到拔节，库中至少保有“肥料需供吸氮量”的 15%（冬前分蘖、返青）
    n_cru=0.42,       # 包膜尿素含氮量
    n_urea=0.46,
)

# 单位籽粒需氮量（kg N/t，含秸秆），拔节、开花时的累计吸氮比例（从施肥日起算）
CROP = {
    "WW": dict(nreq=28.0, pPI=0.30, pH=0.85, label="冬小麦（中强筋）"),
    "SOFT": dict(nreq=25.0, pPI=0.32, pH=0.88, label="弱筋/软质小麦"),
    "SW": dict(nreq=29.0, pPI=0.25, pH=0.82, label="春小麦（硬红春）"),
    "DUR": dict(nreq=29.0, pPI=0.28, pH=0.83, label="硬粒小麦"),
}
EST = {
    "SOW": dict(label="秋冬播种时基施", k_mult=1.0),
    "SPR": dict(label="春播时基施", k_mult=1.0),
    "LW": dict(label="返青前一次施用", k_mult=1.0),
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


def freeze_ramp(T):
    return np.clip(np.asarray(T) / 3.0, 0.0, 1.0)


def soil_temp(Tair, dT):
    T = np.asarray(Tair, dtype=float)
    return np.where(T >= 0, T + dT, T * 0.4 + dT * 0.5)


def cum_tau(Tsoil, Ea, f_soil):
    """第 t 天结束时累计的 25 °C 等效天数（冻土不释放）。"""
    return np.cumsum(arrhenius(Tsoil, Ea) * freeze_ramp(Tsoil) * f_soil)


def loss_rate(Tsoil, k20, wet=1.0):
    T = np.maximum(np.asarray(Tsoil, dtype=float), 0.0)
    return k20 * wet * 2.0 ** ((T - 20.0) / 10.0)


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
    start = doy(s["app"])
    off = lambda mmdd: (doy(mmdd) - start) % 365
    L = off(s["harv"])
    n = L + 1
    Tair = daily_temps(zone["T"], start, n)
    Tsoil = soil_temp(Tair, prm["soil_dT"])
    gdd = np.concatenate([[0], np.cumsum(np.clip(np.minimum(Tair, 30), 0, None))])[:n]
    x = gdd / gdd[-1]
    tPI, tH = off(s["joint"]), off(s["anth"])     # 拔节、开花
    c, e = CROP[s["type"]], EST[s["mode"]]
    pPI, pH = c["pPI"], c["pH"]
    if s["mode"] == "LW":                          # 返青前施用：冬前吸氮不计入
        pPI = max(0.15, pPI - 0.12)
    x0, k = fit_logistic(x[tPI], pPI, x[tH], pH)
    Fcrop = logistic_norm(x, x0, k)
    uptake_total = s["Y"] * c["nreq"]
    fert_demand = max(MIN_FERT_DEMAND, uptake_total - s["INS"])
    d = np.diff(Fcrop, prepend=0.0) * fert_demand
    sow_off = (doy(s["sow"]) - start) % 365
    sow_off = sow_off - 365 if sow_off > 0 else sow_off   # 播种日在施肥日之前或同一天
    return dict(L=L, n=n, start=start, Tair=Tair, Tpaddy=Tsoil, gdd=gdd,
                tPI=tPI, tH=tH, gfill=L - tH, Fcrop=Fcrop, d=d, uptake_total=uptake_total,
                fert_demand=fert_demand, x0=x0, k=k, pPI=pPI, pH=pH, sow_off=sow_off,
                wet=s.get("wet", 1.0) * e["k_mult"],
                k_loss=loss_rate(Tsoil, prm["k_loss"], s.get("wet", 1.0) * e["k_mult"]))


# ---------------- 氮库线性模型 ----------------
def decay_conv(inc, k):
    """pool[t] = Σ_{s≤t} inc[s]·(1−k)^(t−s)"""
    out = np.empty_like(inc)
    acc = 0.0
    q = 1 - np.broadcast_to(np.asarray(k, dtype=float), inc.shape)
    for i, v in enumerate(inc):
        acc = acc * q[i] + v
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
    Tsoil = ss["Tpaddy"] + dT
    k = loss_rate(Tsoil, prm["k_loss"], ss["wet"]) if dT else ss["k_loss"]
    tau = cum_tau(Tsoil, Ea, f_soil)
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
