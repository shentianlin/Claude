"""玉米版：田间释放模型 + 玉米需氮模型 + 掺混配方优化（氮库质量平衡 + 线性规划）。

与小麦版的差别
  - 释放直接用田间释放模型（../common/field_release.py）：包膜处土温按种植系统预设（露地 +1 °C；
    地膜前 60 天 +3.5 °C，之后 40 天内降到 +1 °C），对昼夜温差积分，土温 0–3 °C 线性停释；
    土壤水分因子按各区水分状态（雨养湿润、半干旱、灌溉）；f_soil = 0.90，破损率 2% 计入初期溶出。
  - 需氮曲线用积温（基温 10 °C、上限 30 °C）进度的 S 形曲线，由 V10（约 10 叶期）与吐丝期（R1）
    的累计吸氮比例标定（Ciampitti & Vyn 2012、2013；Bender et al. 2013）。
  - “管饱”约束：施肥后到 V10，氮库不低于肥料需供吸氮量的 10%；之后不低于未来 7 天的需求。
  - 氮库日损失率随土温变化（同小麦版），多雨地区 wet > 1。
  - 同一配方须在常年、偏冷 1.5 °C、偏暖 1.5 °C 三种年份都不断顿。
"""
import itertools
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import linprog, minimize

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "common"))
import field_release as FR  # noqa: E402

PRODUCTS = list(range(30, 361, 10))  # 控释期（天）

P = dict(
    Ea=FR.DEFAULTS["Ea"],               # kJ/mol，Arrhenius 表观活化能
    f_soil=FR.DEFAULTS["f_soil"],       # 土壤中相对静水的基础释放系数
    f0=FR.DEFAULTS["f0"],               # 初期溶出比例（24 h）
    damage=FR.DEFAULTS["damage"],       # 搬运掺混造成的破损率（破损颗粒施后即溶）
    moisture_sens=FR.DEFAULTS["moisture_sens"],
    beta=1.3,                           # Weibull 形状参数
    urea_tau=3.0,                       # 尿素水解时间常数（天）
    urea_loss=0.10,                     # 基施尿素的初始损失（混施入土）
    k_loss=0.010,                       # 20 °C 时有效氮库日损失率
    eta=0.80,                           # 根系最大截获效率
    buffer=7,                           # 库中至少保有未来 7 天的需求
    early_pool=0.10,                    # 施肥后到 V10，库中至少保有“肥料需供吸氮量”的 10%
    n_cru=0.42, n_urea=0.46,
)

# 单位籽粒需氮量（kg N/t，含秸秆），V10、吐丝时的累计吸氮比例
CROP = {"M": dict(nreq=19.0, pPI=0.25, pH=0.60, label="玉米（杂交种）")}
MIN_FERT_DEMAND = 15.0
GDD_BASE, GDD_CAP = 10.0, 30.0


# ---------------- 气温与日期 ----------------
def daily_temps(T12, start_doy, n):
    return FR.daily_temps(np.asarray(T12, float), start_doy, n)


def doy(mmdd):
    m, d = map(int, mmdd.split("-"))
    return [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334][m - 1] + d - 1


# ---------------- 释放模型 ----------------
def f0_eff(f0=None):
    f0 = P["f0"] if f0 is None else f0
    return f0 + (1 - f0) * P["damage"]


def release(tau, D, f0=None, beta=None):
    beta = P["beta"] if beta is None else beta
    return FR.release(tau, D, f0_eff(f0), beta)


def coating_temp(ss, dT=0.0):
    return FR.soil_temperature(ss["Tair"] + dT, ss["system"], "incorporated")


def cum_tau(Tsoil, Ea, f_soil, ss):
    """第 t 天结束时累计的 25 °C 等效天数：昼夜温差积分的温度因子 × 水分因子 × f_soil。"""
    return np.cumsum(FR.temp_factor(Tsoil, ss["amp"], Ea) * ss["aW"] * f_soil)


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
    sysd = FR.SYSTEMS[s["system"]]
    psi = FR.matric_potential(s["moist"], n)
    ss = dict(Tair=Tair, system=s["system"], moist=s["moist"], amp=sysd["amp"],
              aW=FR.moisture_factor(psi, prm["moisture_sens"]))
    Tsoil = coating_temp(ss)
    gdd = np.concatenate([[0], np.cumsum(np.clip(np.minimum(Tair, GDD_CAP) - GDD_BASE, 0, None))])[:n]
    x = gdd / gdd[-1]
    tPI, tH = off(s["v10"]), off(s["r1"])
    c = CROP[s["type"]]
    x0, k = fit_logistic(x[tPI], c["pPI"], x[tH], c["pH"])
    Fcrop = logistic_norm(x, x0, k)
    uptake_total = s["Y"] * c["nreq"]
    fert_demand = max(MIN_FERT_DEMAND, uptake_total - s["INS"])
    d = np.diff(Fcrop, prepend=0.0) * fert_demand
    ss.update(L=L, n=n, start=start, Tpaddy=Tsoil, gdd=gdd, tPI=tPI, tH=tH, gfill=L - tH, Fcrop=Fcrop, d=d,
              uptake_total=uptake_total, fert_demand=fert_demand, x0=x0, k=k, pPI=c["pPI"], pH=c["pH"],
              sow_off=0, wet=s["wet"], k_loss=loss_rate(Tsoil, prm["k_loss"], s["wet"]))
    return ss


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


def buffer_need(ss, prm=P):
    buf = np.array([ss["d"][t + 1:t + 1 + prm["buffer"]].sum() for t in range(ss["n"])])
    early = np.zeros(ss["n"])
    early[1:ss["tPI"]] = prm["early_pool"] * ss["fert_demand"]
    return np.maximum(buf, early) / prm["eta"]


def season_basis(ss, Ea=None, f_soil=None, dT=0.0, prm=P):
    """返回每种产品 1 kg N 在第 t 天末对氮库的贡献 G_k(t)，以及需求侧约束 need(t)。"""
    Ea = prm["Ea"] if Ea is None else Ea
    f_soil = prm["f_soil"] if f_soil is None else f_soil
    Tsoil = coating_temp(ss, dT)
    k = loss_rate(Tsoil, prm["k_loss"], ss["wet"]) if dT else ss["k_loss"]
    tau = cum_tau(Tsoil, Ea, f_soil, ss)
    G = {"urea": decay_conv(urea_increment(ss["n"], prm), k)}
    Fend = {}
    for D in PRODUCTS:
        G[D] = decay_conv(cr_increment(tau, D), k)
        Fend[D] = float(release(tau[-1:], D)[0])
    Dcum = decay_conv(ss["d"] / prm["eta"], k)
    return G, Dcum + buffer_need(ss, prm), Fend, tau


def solve_lp(G, need, keys, cost):
    A = -np.column_stack([G[k] for k in keys])
    c = np.array([cost[k] for k in keys])
    r = linprog(c, A_ub=A, b_ub=-need, bounds=[(0, None)] * len(keys), method="highs")
    return None if r.status != 0 else r.x


ROBUST_DT = (-1.5, 0.0, 1.5)   # 冷年、常年、暖年同时满足“管饱”


def robust_basis(ss, Ea=None, f_soil=None, prm=P):
    parts = [season_basis(ss, Ea, f_soil, dT=dT, prm=prm) for dT in ROBUST_DT]
    G = {k: np.concatenate([p[0][k] for p in parts]) for k in parts[0][0]}
    need = np.concatenate([p[1] for p in parts])
    _, _, Fend, tau = parts[ROBUST_DT.index(0.0)]
    return G, need, Fend, tau


def optimize(ss, Ea=None, f_soil=None, max_cr=2, prm=P, cost_cr=1.02):
    G, need, Fend, tau = robust_basis(ss, Ea, f_soil, prm=prm)
    cands = [D for D in PRODUCTS if Fend[D] >= 0.30]
    cost = {"urea": 1.0, **{D: cost_cr for D in PRODUCTS}}
    res = []
    for D1 in cands:
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
    xall = solve_lp(G, need, ["urea"] + cands, cost)
    ideal = float(xall.sum()) if xall is not None else None
    return best, best1, ideal, G, need, Fend, tau


def round_recipe(keys, x, G, need, step=5):
    """比例取 5% 整数倍；在 LP 最优解附近枚举取整方案（含把小组分去掉），选总量最小者（总量取整到 5 kg N）。

    与稻麦版不同，占比 2.5–8% 的组分也保留到枚举里：玉米的最优解常只有约 5% 尿素，
    直接去掉会让苗期供氮不足、总氮量大增（东北春玉米由 175 升到 240 kg N）。"""
    tot = x.sum()
    parts = [(k, v / tot) for k, v in zip(keys, x) if v / tot >= 0.025]
    ks = [k for k, _ in parts]
    base = [max(step, round(100 * f / step) * step) for _, f in parts]
    best = None
    for delta in itertools.product(range(-2 * step, 2 * step + 1, step), repeat=len(ks) - 1):
        pct = [b + d for b, d in zip(base[:-1], delta)]
        pct.append(100 - sum(pct))
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
    unrel = sum(total * p / 100 * (1 - Fend[k]) for k, p in recipe if k != "urea")
    return dict(
        pool=pool_supply - (need - buffer_need(ss, prm)),
        margin=margin,
        short_days=int(np.sum(margin < -0.5)),
        RE=float(ss["fert_demand"] / total),
        unreleased_pct=float(100 * unrel / total),
        released=float(total - unrel),
    )
