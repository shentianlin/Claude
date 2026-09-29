"""S 型产品情景：在同一“管饱、总氮最少”框架下，改用 S 型（有滞后期）包膜尿素的单档 / 两档方案。

两张地图（水稻、小麦）共用。每个季别在现行“线性型”方案之外，再算两套可选方案：
  s1  尿素 + 1 个 S 型档位
  s2  尿素 + 至多 2 个 S 型档位（取整后若只剩 1 档，则与 s1 相同）
S 型与线性型的控释期同口径（25 °C 静水累计释放 80% 的天数），只是曲线形状参数 β 不同。

各作物模型通过 adapter 接入，adapter 需提供：
  P                       模型参数字典（含可临时修改的 "beta"）
  best_recipe(ss, Ea, beta) -> (dict(rec,total,simple,simple_total,pair,pair_total,ideal), G, need, Fend, tau)
  nominal_basis(ss)       -> (G, need, Fend, tau)   常年气温下的氮库贡献（用于模拟轨迹）
  season_basis(ss, dT)    -> (G, need, Fend, tau)   气温整体偏移 dT 时的基础量
  simulate(ss, rec, total, G, need, Fend) -> dict(pool, RE, unreleased_pct, ...)
  buffer_need(ss), cum_tau(T, Ea, f_soil), release(tau, D, beta=), urea_increment(n)
"""
from contextlib import contextmanager

import numpy as np

SCEN_LABELS = {"lin": "线性型（现有产品）", "s1": "S 型单档", "s2": "S 型两档"}


def rec_str(rec, prefix="CR"):
    return " + ".join(("尿素" if k == "urea" else f"{prefix}{k}") + f" {p}%" for k, p in rec)


@contextmanager
def beta_set(P, beta):
    b0 = P["beta"]
    P["beta"] = beta
    try:
        yield
    finally:
        P["beta"] = b0


def _fields(ad, ss, rec, total, simple, simple_total, ideal, sens, beta, idx, robust_dT, prefix):
    P = ad.P
    with beta_set(P, beta):
        G0, need0, Fend0, _ = ad.nominal_basis(ss)
        sim = ad.simulate(ss, rec, total, G0, need0, Fend0)
        robust = {}
        for dT in robust_dT:
            G2, need2, Fend2, _ = ad.season_basis(ss, dT)
            sup = sum(total * p / 100 * G2[k] for k, p in rec)
            margin = sup - need2
            pool = sup - (need2 - ad.buffer_need(ss))
            robust[str(dT)] = dict(
                short_days=int(np.sum(pool < -0.5)),
                min_buffer_pct=float(100 * np.min(margin / np.maximum(need2, 1e-6))),
                unreleased_pct=float(100 * sum(total * p / 100 * (1 - Fend2[k]) for k, p in rec if k != "urea") / total))
    tau_all = ad.cum_tau(ss["Tpaddy"], P["Ea"], P["f_soil"])
    comps = {}
    for k, p in rec:
        kg = total * p / 100
        if k == "urea":
            cum = np.cumsum(ad.urea_increment(ss["n"])) / (1 - P["urea_loss"])
            name = "尿素"
        else:
            cum = ad.release(tau_all, k, beta=beta)
            name = f"{prefix}{k}"
        comps[name] = [round(float(kg * cum[i]), 1) for i in idx]
    n_cr = sum(p for k, p in rec if k != "urea") / 100
    kg_cr = total * n_cr / P["n_cru"]
    kg_u = total * (1 - n_cr) / P["n_urea"]
    return dict(
        rec=rec, rec_str=rec_str(rec, prefix), total=total,
        simple=simple, simple_str=rec_str(simple, prefix), simple_total=simple_total,
        ideal=round(ideal) if ideal else None,
        kg_cr=round(kg_cr), kg_urea=round(kg_u), kg_product=round(kg_cr + kg_u),
        RE=round(sim["RE"], 3), unreleased=round(sim["unreleased_pct"], 1),
        sens=sens, robust=robust, beta=beta,
        series_comps=comps, series_pool=[round(float(sim["pool"][i]), 1) for i in idx],
    )


def s_type_alternatives(ad, ss, idx, beta_s, ea_list, robust_dT):
    """返回 {"s1": {...}, "s2": {...}}，字段与页面上线性型方案的字段一一对应。"""
    base, *_ = ad.best_recipe(ss, None, beta_s)
    sens1, sens2 = {}, {}
    for Ea in ea_list:
        r, *_ = ad.best_recipe(ss, Ea, beta_s)
        sens1[str(Ea)] = dict(rec_str=rec_str(r["simple"], "S"), total=r["simple_total"])
        sens2[str(Ea)] = dict(rec_str=rec_str(r["pair"], "S"), total=r["pair_total"])
    s1 = _fields(ad, ss, base["simple"], base["simple_total"], base["simple"], base["simple_total"],
                 base["ideal"], sens1, beta_s, idx, robust_dT, "S")
    s2 = _fields(ad, ss, base["pair"], base["pair_total"], base["simple"], base["simple_total"],
                 base["ideal"], sens2, beta_s, idx, robust_dT, "S")
    s2["same_as_s1"] = [list(x) for x in base["pair"]] == [list(x) for x in base["simple"]]
    return {"s1": s1, "s2": s2}
