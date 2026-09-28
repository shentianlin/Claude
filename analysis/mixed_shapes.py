"""线性型 + S 型混合掺混的评估。

对每个稻季、麦季，在同一“管饱、总氮最少”的线性规划框架下比较：
  A  尿素 + 1 个线性型档位            （地图现行方案的结构）
  B  尿素 + 2 个线性型档位
  C  尿素 + 1 个 S 型档位
  F  尿素 + 2 个 S 型档位
  D  尿素 + 1 个线性型 + 1 个 S 型     （线性型管前期，S 型管后期/越冬后）
  E  尿素 + 全部线性型 + 全部 S 型档位  （任意多档组合的理论下限）
线性型 β=1.3（现有产品）；S 型 β 默认 2.5，另对 2.0、3.5 只算理论下限 E。
两类产品的控释期定义相同：25 °C 静水累计释放 80% 的天数。
水稻按常年气温；小麦沿用地图的稳健约束（偏冷 1.5 °C、常年、偏暖 1.5 °C 同时不断顿）。

运行：python analysis/mixed_shapes.py   （输出 analysis/mixed_shapes_results.csv 与汇总）
"""
import csv
import json
import statistics as st
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "rice-crf-map"))
sys.path.insert(0, str(ROOT / "wheat-crf-map"))
sys.path.insert(0, str(ROOT / "common"))

import crfmap.optimize as RO               # noqa: E402  水稻
from crfmap.data import ZONES as RZ        # noqa: E402
from crfmap.demand import season_setup as r_setup  # noqa: E402
from crfmap.params import P as RP          # noqa: E402
import model as WM                         # noqa: E402  小麦
from zones import ZONES as WZ              # noqa: E402

BETA_L = 1.3
BETA_S = 2.5
BETA_S_EXTRA = (2.0, 3.5)
COST = {"urea": 1.0, "cr": 1.02}


def family_basis(crop, ss, beta):
    """一个形状族（给定 β）全部档位的氮库贡献 G、需求约束 need、成熟时释放比例 Fend。"""
    P = RP if crop == "rice" else WM.P
    b0 = P["beta"]
    P["beta"] = beta
    try:
        if crop == "rice":
            G, need, Fend, tau = RO.season_basis(ss)
        else:
            G, need, Fend, tau = WM.robust_basis(ss)
    finally:
        P["beta"] = b0
    return G, need, Fend, tau


def lp(G, need, keys):
    A = -np.column_stack([G[k] for k in keys])
    c = np.array([COST["urea"] if k == "urea" else COST["cr"] for k in keys])
    from scipy.optimize import linprog
    r = linprog(c, A_ub=A, b_ub=-need, bounds=[(0, None)] * len(keys), method="highs")
    return None if r.status != 0 else r.x


def best_of(G, need, combos):
    best = None
    for keys in combos:
        x = lp(G, need, keys)
        if x is None:
            continue
        cost = sum((COST["urea"] if k == "urea" else COST["cr"]) * v for k, v in zip(keys, x))
        if best is None or cost < best[0]:
            best = (cost, keys, x)
    return best


def analyse(args):
    crop, zi, si = args
    if crop == "rice":
        z = RZ[zi]; s = z["seasons"][si]; ss = r_setup(z, s); rnd = RO.round_recipe
        frost = False
    else:
        z = WZ[zi]; s = z["seasons"][si]; ss = WM.season_setup(z, s); rnd = WM.round_recipe
        frost = bool(np.min(ss["Tpaddy"]) < 3.0)        # 有土温 < 3 °C 的冻土/休眠期
    GL, need, FL, tau = family_basis(crop, ss, BETA_L)
    GS, _, FS, _ = family_basis(crop, ss, BETA_S)
    G = {"urea": GL["urea"]}
    L = [("L", D) for D in GL if D != "urea" and FL[D] >= 0.30]
    S = [("S", D) for D in GS if D != "urea" and FS[D] >= 0.30]
    for (_, D) in L:
        G[("L", D)] = GL[D]
    for (_, D) in S:
        G[("S", D)] = GS[D]
    res = {}
    res["A"] = best_of(G, need, [["urea", k] for k in L])
    res["B"] = best_of(G, need, [["urea", a, b] for i, a in enumerate(L) for b in L[i + 1:]])
    res["C"] = best_of(G, need, [["urea", k] for k in S])
    res["F"] = best_of(G, need, [["urea", a, b] for i, a in enumerate(S) for b in S[i + 1:]])
    res["D"] = best_of(G, need, [["urea", a, b] for a in L for b in S])
    xE = lp(G, need, ["urea"] + L + S)
    out = dict(crop=crop, zone=z["id"], zone_name=z["name"], season=s["name"], frost=frost,
               mode=s.get("mode", s.get("est")), days=ss["L"])
    for k in ("A", "B", "C", "F", "D"):
        r = res[k]
        out[k] = round(float(r[2].sum()), 1) if r else None
        if r:
            rec, tot = rnd(r[1], r[2], G, need)
            out[k + "_round"] = tot
            out[k + "_rec"] = " + ".join(("尿素" if kk == "urea" else f"{'线' if kk[0] == 'L' else 'S'}{kk[1]}") + f" {p}%"
                                         for kk, p in rec)
    out["E"] = round(float(xE.sum()), 1) if xE is not None else None
    if xE is not None:
        keysE = ["urea"] + L + S
        out["E_n_used"] = int(sum(1 for k, v in zip(keysE, xE) if k != "urea" and v > 0.5))
    # S 型形状更强或更弱时的理论下限
    for b in BETA_S_EXTRA:
        Gb, _, Fb, _ = family_basis(crop, ss, b)
        Gx = {"urea": GL["urea"], **{("L", D): GL[D] for (_, D) in L},
              **{("S", D): Gb[D] for D in Gb if D != "urea" and Fb[D] >= 0.30}}
        x = lp(Gx, need, list(Gx))
        out[f"E_beta{b}"] = round(float(x.sum()), 1) if x is not None else None
    # 混合方案中线性、S 型两部分在“休眠期前/后”的释放（仅小麦有冻土期时）
    if crop == "wheat" and frost and res["D"]:
        cold = np.where(ss["Tpaddy"] < 3.0)[0]
        t_end = int(cold.max())                            # 冻土/休眠期结束
        keys, x = res["D"][1], res["D"][2]
        kl, ks = keys[1], keys[2]
        tauW = WM.cum_tau(ss["Tpaddy"], WM.P["Ea"], WM.P["f_soil"])
        rl = float(WM.release(tauW[t_end:t_end + 1], kl[1], beta=BETA_L)[0])
        rs = float(WM.release(tauW[t_end:t_end + 1], ks[1], beta=BETA_S)[0])
        out["D_linear_released_before_regreen"] = round(100 * rl)
        out["D_S_released_before_regreen"] = round(100 * rs)
    return out


def main():
    jobs = [("rice", zi, si) for zi, z in enumerate(RZ) for si in range(len(z["seasons"]))]
    jobs += [("wheat", zi, si) for zi, z in enumerate(WZ) for si in range(len(z["seasons"]))]
    with Pool() as p:
        rows = p.map(analyse, jobs)
    keys = sorted({k for r in rows for k in r}, key=lambda k: (k not in ("crop", "zone", "zone_name", "season"), k))
    with open(ROOT / "analysis" / "mixed_shapes_results.csv", "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)

    def summ(sel, label):
        if not sel:
            return
        pct = lambda k, base="A": [100 * (r[base] - r[k]) / r[base] for r in sel if r.get(k) is not None]
        kg = lambda k, base="A": [r[base] - r[k] for r in sel if r.get(k) is not None]
        line = f"{label}（{len(sel)} 季）"
        for k, n in (("B", "2 线性"), ("C", "1 S"), ("F", "2 S"), ("D", "线性+S"), ("E", "任意多档")):
            v = pct(k)
            line += f" | {n}: 中位 {st.median(v):.1f}% 最大 {max(v):.1f}%"
        print(line)
        v = pct("D", "C")
        print(f"    线性+S 相对 1 个 S 型再省：中位 {st.median(v):.1f}%  最大 {max(v):.1f}%；"
              f"线性+S 比现方案省 ≥5% 的季数 {sum(x >= 5 for x in pct('D'))}；省氮 kg 中位 {st.median(kg('D')):.1f} 最大 {max(kg('D')):.1f}")
        rr = [100 * (r["A_round"] - r["D_round"]) / r["A_round"] for r in sel if r.get("D_round")]
        print(f"    按 5% 取整后的实用配方：线性+S 比现方案省 中位 {st.median(rr):.1f}%  最大 {max(rr):.1f}%")
        for b in BETA_S_EXTRA:
            v = pct(f"E_beta{b}")
            print(f"    若 S 型 β={b}：任意多档比现方案省 中位 {st.median(v):.1f}%  最大 {max(v):.1f}%")

    rice = [r for r in rows if r["crop"] == "rice"]
    wheat = [r for r in rows if r["crop"] == "wheat"]
    summ(rice, "水稻 全部")
    summ(wheat, "小麦 全部")
    summ([r for r in wheat if r["frost"] and r["mode"] == "SOW"], "小麦 秋播且有越冬冻土期")
    summ([r for r in wheat if not r["frost"] and r["mode"] == "SOW"], "小麦 秋冬播、无冻土期")
    summ([r for r in wheat if r["mode"] == "LW"], "小麦 返青前施用（西欧）")
    summ([r for r in wheat if r["mode"] == "SPR"], "小麦 春播")
    print("\n线性+S 省氮最多的 10 季：")
    for r in sorted(rows, key=lambda r: -(r["A"] - r["D"]) / r["A"])[:10]:
        extra = ""
        if "D_S_released_before_regreen" in r:
            extra = f"；返青前已释放：线性 {r['D_linear_released_before_regreen']}%、S 型 {r['D_S_released_before_regreen']}%"
        print(f"  {r['crop']} {r['zone_name']} {r['season']}: 现方案 {r['A_round']}（{r['A_rec']}） → "
              f"线性+S {r['D_round']}（{r['D_rec']}），理论下限 {r['E']}{extra}")


if __name__ == "__main__":
    main()
