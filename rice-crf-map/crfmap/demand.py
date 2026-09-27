"""水稻需氮模型。

作物累计吸氮比例是“有效积温进度” x = GDD(t)/GDD(成熟) 的 S 形曲线，
用幼穗分化期、抽穗期两点的累计吸氮比例标定。
作物总吸氮 = 目标产量 × 单位籽粒需氮量；土壤基础供氮 INS 按同一曲线分摊，
其余部分必须由肥料供给。
"""
import numpy as np
from scipy.optimize import minimize

from .climate import daily_temps, doy
from .params import CROP, EST, MIN_FERT_DEMAND, P, PHEN


def logistic_norm(x, x0, k):
    """归一化 logistic：x=0 时为 0，x=1 时为 1。"""
    L = lambda z: 1 / (1 + np.exp(-k * (z - x0)))
    return (L(x) - L(0)) / (L(1) - L(0))


def fit_logistic(xPI, pPI, xH, pH):
    """求 (x0, k)，使曲线通过（幼穗分化, pPI）与（抽穗, pH）两点。"""
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
    """一个稻季的逐日气温、生育期、吸氮曲线与每日肥料氮需求。"""
    L = s["days"]
    n = L + 1                                   # 第 0 … L 天
    start = doy(s["app"])
    Tair = daily_temps(zone["T"], start, n)
    gdd = np.concatenate([[0], np.cumsum(np.clip(np.minimum(Tair, PHEN["gdd_cap"]) - PHEN["gdd_base"], 0, None))])[:n]
    x = gdd / gdd[-1]
    Tm_end = Tair[-35:].mean()
    gfill = int(round(np.clip(PHEN["gfill_base"] + (PHEN["gfill_ref_T"] - Tm_end) * PHEN["gfill_slope"],
                              PHEN["gfill_min"], PHEN["gfill_max"])))
    tH = L - gfill
    tPI = max(tH - PHEN["pi_to_heading"], PHEN["pi_min"])
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
