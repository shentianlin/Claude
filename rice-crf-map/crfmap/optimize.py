"""肥料氮库质量平衡 + 线性规划优化。

    库(t+1) = 库(t)·(1 − k_loss) + 当日释放 − 当日肥料氮吸收 / η
约束“管饱”：库(t) ≥ max(分蘖期库存下限, 未来 buffer 天的需求)，任何一天都不断顿。
目标“不浪费”：总施氮量最小（成熟时仍未释放的氮、过早释放而损失的氮都计入总量）。
在控释期档位中枚举“尿素 + 至多 max_cr 个控释期”，每个组合解一个线性规划。
"""
import itertools
import math

import numpy as np
from scipy.optimize import linprog

from .params import OPT, P, PRODUCTS, SCEN
from .release import cr_increment, cum_tau, release, urea_increment


def decay_conv(inc, k):
    """pool[t] = Σ_{s≤t} inc[s]·(1−k)^(t−s)"""
    out = np.empty_like(inc)
    acc = 0.0
    q = 1 - k
    for i, v in enumerate(inc):
        acc = acc * q + v
        out[i] = acc
    return out


def buffer_need(ss, prm=P):
    """“管饱”下限：分蘖期库存下限与未来 buffer 天需求中的较大者（按 η 折算）。"""
    buf = np.array([ss["d"][t + 1:t + 1 + prm["buffer"]].sum() for t in range(ss["n"])])
    early = np.zeros(ss["n"])
    early[1:ss["tPI"]] = prm["early_pool"] * ss["fert_demand"]
    return np.maximum(buf, early) / prm["eta"]


def season_basis(ss, Ea=None, f_soil=None, dT=0.0, prm=P):
    """返回每种产品 1 kg N 在第 t 天末对氮库的贡献 G_k(t)、需求侧约束 need(t)、成熟时释放比例、τ。"""
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


def robust_basis(ss, dTs, Ea=None, f_soil=None, prm=P):
    """把多个气温偏移年份的约束纵向拼接：同一配方须在这些年份都不断顿。"""
    parts = [season_basis(ss, Ea, f_soil, dT=dT, prm=prm) for dT in dTs]
    G = {k: np.concatenate([p[0][k] for p in parts]) for k in parts[0][0]}
    need = np.concatenate([p[1] for p in parts])
    base = parts[list(dTs).index(0.0)] if 0.0 in dTs else parts[0]
    return G, need, base[2], base[3]


def solve_lp(G, need, keys, cost):
    A = -np.column_stack([G[k] for k in keys])
    c = np.array([cost[k] for k in keys])
    r = linprog(c, A_ub=A, b_ub=-need, bounds=[(0, None)] * len(keys), method="highs")
    if r.status != 0:
        return None
    return r.x


def optimize(ss, Ea=None, f_soil=None, max_cr=None, prm=P, cost_cr=None, robust_dT=None):
    """返回 (最优组合, 最优单一控释期组合, 全档位理论下限, G, need, Fend, tau)。
    robust_dT 给出时（如 (-1.5, 0.0, 1.5)），要求配方在这些气温偏移年份都不断顿。"""
    max_cr = OPT["max_cr"] if max_cr is None else max_cr
    cost_cr = OPT["cost_cr"] if cost_cr is None else cost_cr
    if robust_dT:
        G, need, Fend, tau = robust_basis(ss, robust_dT, Ea, f_soil, prm=prm)
    else:
        G, need, Fend, tau = season_basis(ss, Ea, f_soil, prm=prm)
    cands = [D for D in PRODUCTS if Fend[D] >= OPT["min_release_at_maturity"]]
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
    xall = solve_lp(G, need, ["urea"] + cands, cost)      # 理论下限：全部档位可用
    ideal = float(xall.sum()) if xall is not None else None
    return best, best1, ideal, G, need, Fend, tau


def round_recipe(keys, x, G, need, step=None):
    """比例取 step% 整数倍、去掉占比过小的组分；在 LP 最优解附近枚举取整方案，选总量最小者。"""
    step = OPT["round_step"] if step is None else step
    tot = x.sum()
    parts = [(k, v / tot) for k, v in zip(keys, x) if v / tot >= OPT["min_component"]]
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
    ts = OPT["total_step"]
    return best[0], ts * math.ceil(best[1] / ts)


def best_recipe(ss, Ea=None, beta=None, robust_dT=None):
    """优化并取整；单一控释期方案多用的氮不超过容差时优先选它（生产上 SKU 越少越好）。"""
    beta0 = P["beta"]
    if beta is not None:
        P["beta"] = beta
    try:
        best, best1, ideal, G, need, Fend, tau = optimize(ss, Ea=Ea, robust_dT=robust_dT)
        r2, t2 = round_recipe(best[1], best[2], G, need)
        r1, t1 = round_recipe(best1[1], best1[2], G, need)
    finally:
        P["beta"] = beta0
    pair, pair_total = r2, t2               # 不做“单档优先”取舍的两档方案（供 S 型两档情景使用）
    if t1 <= t2 * OPT["single_cr_tolerance"]:
        r2, t2 = r1, t1
    return dict(rec=r2, total=t2, simple=r1, simple_total=t1, ideal=ideal,
                pair=pair, pair_total=pair_total), G, need, Fend, tau


def split_urea(ss, need, fr=None):
    """对照：同一“管饱”标准下，三次分施尿素（基肥、分蘖肥、穗肥）所需的最少总氮量。"""
    fr = SCEN["split_fractions"] if fr is None else fr
    k, n = ss["k_loss"], ss["n"]
    days = [0, max(SCEN["split_tiller_min_day"], ss["tPI"] - SCEN["split_tiller_before_pi"]), ss["tPI"]]
    base = urea_increment(n)
    shape = np.zeros(n)
    for f, d0 in zip(fr, days):
        inc = np.zeros(n)
        inc[d0:] = base[: n - d0]
        shape += f * decay_conv(inc, k)
    ts = OPT["total_step"]
    return ts * math.ceil(float(np.max(need / np.maximum(shape, 1e-9))) / ts), days


def simulate(ss, recipe, total, G, need, Fend, prm=P):
    """给定配方，计算氮库轨迹、断顿天数、模拟回收率与成熟时未释放比例。"""
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
