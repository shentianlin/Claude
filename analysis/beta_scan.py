"""控释肥形状参数 β 扫描：用田间释放模型（common/field_release.py）找对产量、环境最有利的 β。

对每个稻季（73）、麦季，β 取 1.0–5.0，每个 β 都让产品自选最合适的控释期 D 和尿素比例
（尿素 + 1 个控释档，与地图推荐结构相同），然后比较两件事：

一、等氮比较（看增产）
    施氮量固定为“现有线性产品（β=1.3）刚好不断顿所需的氮量”的 100% 和 85%（减氮 15%）。
    每个 β 按常年、偏冷 1.5 °C、偏暖 1.5 °C 三种年份的平均吸氮缺口最小来选 D 和尿素比例；
    再在 5 种年份/田间条件下检验：常年、偏冷 3 °C、偏暖 3 °C、田间释放慢 25%、快 33%
    （后两者代表土壤、水分、施肥位置带来的不确定性）。
    吸氮缺口 = 氮库不够时作物少吸的氮。产量损失按“缺口 ÷ 作物总吸氮”估算
    （缺氮范围内产量大致与吸氮量成正比，是粗略的上限估计）。

二、等产比较（看环境）
    每个 β 求“三种年份都不断顿”所需的最少氮量（与地图相同的线性规划），
    再用 common/env.py 计算氨挥发、N₂O、淋溶径流、温室气体和社会损害成本（同一组 Monte Carlo 抽样），
    并给出收获时未释放的氮。

田间释放：水稻用 rice_flooded（旱直播稻季用 rice_dsr），小麦用 wheat_irrigated，施肥方式为混施入土、壤土。
土温按种植系统预设，昼夜温差积分，冻土期停释（与模拟器一致）。

运行：python analysis/beta_scan.py   （约 2–4 分钟；输出 analysis/beta_scan_results.csv、beta_scan.png 和汇总）
"""
import csv
import statistics as st
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "rice-crf-map"))
sys.path.insert(0, str(ROOT / "wheat-crf-map"))
sys.path.insert(0, str(ROOT / "common"))

import crfmap.optimize as RO               # noqa: E402
from crfmap.data import ZONES as RZ        # noqa: E402
from crfmap.demand import season_setup as r_setup  # noqa: E402
from crfmap.params import P as RP          # noqa: E402
import env as ENV                          # noqa: E402
import field_release as FR                 # noqa: E402
import model as WM                         # noqa: E402
from zones import ZONES as WZ              # noqa: E402

BETAS = [1.0, 1.3, 1.6, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0]
DS = list(range(30, 361, 10))
UREA = np.arange(0, 0.61, 0.05)
DESIGN = (-1.5, 0.0, 1.5)
STRESS = [("常年", 0.0, 1.0), ("偏冷3°C", -3.0, 1.0), ("偏暖3°C", 3.0, 1.0),
          ("释放慢25%", 0.0, 0.75), ("释放快33%", 0.0, 1.33)]
N_LEVELS = (1.0, 0.85)


# ---------------------------------------------------------------- 季别与田间释放
def setup(crop, zi, si):
    if crop == "rice":
        z = RZ[zi]; s = z["seasons"][si]; ss = r_setup(z, s); P = RP
        system = "rice_dsr" if s["est"] == "DSR" else "rice_flooded"
    else:
        z = WZ[zi]; s = z["seasons"][si]; ss = WM.season_setup(z, s); P = WM.P
        system = "wheat_irrigated"
    return z, s, ss, P, system


def k_loss(crop, ss, P, dT):
    if crop == "rice" or dT == 0:
        return np.broadcast_to(np.asarray(ss["k_loss"], float), (ss["n"],)).copy()
    return WM.loss_rate(ss["Tpaddy"] + dT, P["k_loss"], ss["wet"])


def field_tau(ss, system, dT=0.0, mult=1.0):
    """田间释放模型的 25 °C 等效天数（第 t 天结束时）。"""
    n = ss["n"]
    Tair = ss["Tair"] + dT
    Ts = FR.soil_temperature(Tair, system, "incorporated")
    sysd = FR.SYSTEMS[system]
    aT = FR.temp_factor(Ts, sysd["amp"], FR.DEFAULTS["Ea"])
    psi = FR.matric_potential(sysd["moisture"], n)
    if "switch" in sysd:
        d0, reg = sysd["switch"]
        psi[d0:] = FR.matric_potential(reg, n - d0)
    aW = FR.moisture_factor(psi, FR.DEFAULTS["moisture_sens"])
    return np.cumsum(aT * aW * FR.DEFAULTS["f_soil"] * mult)


def cr_inc(tau, D, beta):
    f0e = FR.DEFAULTS["f0"] + (1 - FR.DEFAULTS["f0"]) * FR.DEFAULTS["damage"]
    F = FR.release(tau, D, f0e, beta)
    return np.diff(np.concatenate([[0.0], F])), float(F[-1])


def decay(inc, k):
    out = np.empty_like(inc)
    acc = 0.0
    for i, v in enumerate(inc):
        acc = acc * (1 - k[i]) + v
        out[i] = acc
    return out


# ---------------------------------------------------------------- 一、等产（最少氮量）
def design(crop, ss, P, system, beta, mod):
    """三种年份都不断顿的最少氮量：返回 (总氮, D, 尿素比例, 成熟时释放比例)。"""
    parts = []
    for dT in DESIGN:
        k = k_loss(crop, ss, P, dT)
        tau = field_tau(ss, system, dT)
        G = {"urea": decay(mod.urea_increment(ss["n"], P), k)}
        Fend = {}
        for D in DS:
            inc, Fe = cr_inc(tau, D, beta)
            G[D] = decay(inc, k)
            Fend[D] = Fe
        need = decay(ss["d"] / P["eta"], k) + mod.buffer_need(ss, P)
        parts.append((G, need, Fend))
    G = {key: np.concatenate([p[0][key] for p in parts]) for key in parts[0][0]}
    need = np.concatenate([p[1] for p in parts])
    Fend = parts[1][2]
    best = None
    for D in DS:
        if Fend[D] < 0.30:
            continue
        A = -np.column_stack([G["urea"], G[D]])
        r = linprog([1.0, 1.02], A_ub=A, b_ub=-need, bounds=[(0, None)] * 2, method="highs")
        if r.status != 0:
            continue
        tot = float(r.x.sum())
        if best is None or tot < best[0]:
            best = (tot, D, float(r.x[0] / tot), Fend[D])
    return best


# ---------------------------------------------------------------- 二、等氮（吸氮缺口）
def shortfall(N, u, sup_u, sup_cr, d, k, eta):
    """N kg/ha、尿素比例 u（数组，候选方案）时，整季吸氮缺口（kg N/ha）。sup_cr 形状 (候选, 天)。"""
    supply = N * (u[:, None] * sup_u[None, :] + (1 - u[:, None]) * sup_cr)
    p = np.zeros(len(u))
    short = np.zeros(len(u))
    for t in range(len(d)):
        p = p * (1 - k[t]) + supply[:, t]
        draw = d[t] / eta
        got = np.minimum(p, draw)
        p -= got
        short += d[t] - got * eta
    return short


def candidates(crop, ss, P, system, beta, dT, mult, mod):
    k = k_loss(crop, ss, P, dT)
    tau = field_tau(ss, system, dT, mult)
    sup_u = mod.urea_increment(ss["n"], P)
    rows, Ds, us = [], [], []
    for D in DS:
        inc, _ = cr_inc(tau, D, beta)
        for u in UREA:
            rows.append(inc); Ds.append(D); us.append(u)
    return np.array(rows), np.array(Ds), np.array(us), sup_u, k


def equal_n(crop, ss, P, system, beta, N, mod, cache):
    """在设计年份（常年、±1.5 °C）里选平均缺口最小的 (D, 尿素比例)，再检验 5 种条件。"""
    tot = np.zeros(len(DS) * len(UREA))
    for dT in DESIGN:
        key = (beta, dT, 1.0)
        if key not in cache:
            cache[key] = candidates(crop, ss, P, system, beta, dT, 1.0, mod)
        sup, Ds, us, sup_u, k = cache[key]
        tot += shortfall(N, us, sup_u, sup, ss["d"], k, P["eta"])
    # 缺口相同（都为 0）时，选控释期最短的方案（收获时残留最少）
    order = np.lexsort((Ds, np.round(tot, 2)))
    i = order[0]
    D, u = int(Ds[i]), float(us[i])
    res = {}
    for name, dT, mult in STRESS:
        key = (beta, dT, mult)
        if key not in cache:
            cache[key] = candidates(crop, ss, P, system, beta, dT, mult, mod)
        sup, Ds_, us_, sup_u, k = cache[key]
        j = np.where((Ds_ == D) & np.isclose(us_, u))[0][0]
        res[name] = float(shortfall(N, us_[j:j + 1], sup_u, sup[j:j + 1], ss["d"], k, P["eta"])[0])
    return D, u, res


# ---------------------------------------------------------------- 一个季别
def run(args):
    crop, zi, si = args
    z, s, ss, P, system = setup(crop, zi, si)
    mod = RO if crop == "rice" else WM
    frost = bool(np.min(ss["Tpaddy"]) < 3.0)
    out = []
    designs = {b: design(crop, ss, P, system, b, mod) for b in BETAS}
    N_ref = designs[1.3][0]
    rng_seed = 20261002 + zi * 31 + si
    cache = {}
    for b in BETAS:
        tot, D, u, Fe = designs[b]
        rec = [("urea", round(100 * u)), (D, 100 - round(100 * u))]
        senv = dict(rec=rec, total=tot, FN=s["FN"], split_total=s["FN"], Y=s["Y"],
                    days=ss["L"] if crop == "rice" else s.get("days", ss["L"]))
        e = ENV.season_env(crop, z["id"], senv, np.random.default_rng(rng_seed))["scen"]["CRF"]
        row = dict(crop=crop, zone=z["id"], zone_name=z["name"], season=s["name"],
                   mode=s.get("est", s.get("mode")), frost=frost, beta=b, uptake=round(ss["uptake_total"], 1),
                   fert_demand=round(ss["fert_demand"], 1),
                   N_req=round(tot, 1), N_req_pct=round(100 * tot / N_ref, 1), D_req=D, urea_req=round(100 * u),
                   unreleased_N=round(tot * (1 - u) * (1 - Fe), 1),
                   NH3=e["NH3"][1], N2O=e["N2O"][1], LEACH=e["LEACH"][1], Nr=e["Nr"][1],
                   GHG=e["GHG_all"][1], DMG=e["DMG"][1])
        for lv in N_LEVELS:
            D2, u2, res = equal_n(crop, ss, P, system, b, N_ref * lv, mod, cache)
            tag = f"N{int(lv * 100)}"
            row[f"{tag}_D"], row[f"{tag}_urea"] = D2, round(100 * u2)
            for name, v in res.items():
                row[f"{tag}_{name}"] = round(100 * v / ss["uptake_total"], 2)   # 产量损失估计 %
            row[f"{tag}_worst"] = max(row[f"{tag}_{n}"] for n, *_ in STRESS)
            row[f"{tag}_mean"] = round(st.mean(row[f"{tag}_{n}"] for n, *_ in STRESS), 2)
        out.append(row)
    return out


# ---------------------------------------------------------------- 汇总
def summarise(rows, label):
    seasons = {(r["crop"], r["zone"], r["season"]) for r in rows}
    if not seasons:
        return None
    print(f"\n{label}（{len(seasons)} 季）")
    print("  β   | 等产：所需氮(相对β1.3)  NH3   N2O   淋溶  损害成本 | "
          "等氮100%：常年 最差  | 减氮15%：常年  5种平均  最差  | 常用D  尿素%")
    summ = {}
    for b in BETAS:
        rb = [r for r in rows if r["beta"] == b]
        m = lambda k: st.median(r[k] for r in rb)
        mean = lambda k: st.mean(r[k] for r in rb)
        summ[b] = dict(N=m("N_req_pct"), NH3=m("NH3"), N2O=m("N2O"), LEACH=m("LEACH"), DMG=m("DMG"),
                       y100=mean("N100_常年"), y100w=mean("N100_worst"),
                       y85=mean("N85_常年"), y85m=mean("N85_mean"), y85w=mean("N85_worst"),
                       D=m("N85_D"), u=m("N85_urea"))
        s = summ[b]
        print(f"  {b:<4}|      {s['N']:6.1f}%          {s['NH3']:5.1f} {s['N2O']:5.2f} {s['LEACH']:5.1f} "
              f"{s['DMG']:7.0f}  |      {s['y100']:4.2f}% {s['y100w']:5.2f}% |        {s['y85']:4.2f}%  "
              f"{s['y85m']:5.2f}%  {s['y85w']:5.2f}% | {s['D']:5.0f}  {s['u']:4.0f}")
    return summ


def plot(sums, path):
    import logging
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
    plt.rcParams["font.family"] = ["Microsoft YaHei", "PingFang SC", "Noto Sans CJK SC", "Source Han Sans SC",
                                   "WenQuanYi Zen Hei", "SimHei", "Arial Unicode MS", "sans-serif"]
    plt.rcParams["axes.unicode_minus"] = False
    INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
    C = {"水稻": "#2a78d6", "小麦·有冻土期": "#eb6834", "小麦·其他": "#1baf7a"}
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.3), facecolor=SURF, layout="constrained")
    panels = [("减氮 15% 时的产量损失估计（5 种条件平均，%）", "y85m"),
              ("减氮 15% 时最不利条件下的产量损失（%）", "y85w"),
              ("不断顿所需氮量（相对 β=1.3，%）", "N")]
    for ax, (title, key) in zip(axes, panels):
        ax.set_facecolor(SURF)
        for sp in ax.spines.values():
            sp.set_color(GRID)
        ax.tick_params(colors=INK2, labelsize=8.5)
        ax.grid(True, color=GRID, lw=0.7)
        for name, summ in sums.items():
            if summ:
                ax.plot(BETAS, [summ[b][key] for b in BETAS], "o-", color=C[name], lw=1.8, ms=4.5, label=name)
        ax.axvspan(2.0, 3.0, color="#e4e3df", alpha=0.6, lw=0)
        ax.set_title(title, loc="left", fontsize=10, color=INK)
        ax.set_xlabel("β（Weibull 形状参数）", color=INK2, fontsize=9)
    axes[0].legend(frameon=False, fontsize=9)
    fig.savefig(path, dpi=150, facecolor=SURF)
    plt.close(fig)


def main():
    jobs = [("rice", zi, si) for zi, z in enumerate(RZ) for si in range(len(z["seasons"]))]
    jobs += [("wheat", zi, si) for zi, z in enumerate(WZ) for si in range(len(z["seasons"]))]
    with Pool() as p:
        rows = [r for rs in p.map(run, jobs) for r in rs]
    with open(ROOT / "analysis" / "beta_scan_results.csv", "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print("产量损失估计 = 吸氮缺口 ÷ 作物总吸氮；氮、NH3、N2O、淋溶单位 kg N/ha，损害成本 USD/ha（中位数）")
    sums = {}
    sums["水稻"] = summarise([r for r in rows if r["crop"] == "rice"], "水稻")
    sums["小麦·有冻土期"] = summarise([r for r in rows if r["crop"] == "wheat" and r["frost"]], "小麦·有冻土期")
    sums["小麦·其他"] = summarise([r for r in rows if r["crop"] == "wheat" and not r["frost"]], "小麦·无冻土期")
    plot(sums, ROOT / "analysis" / "beta_scan.png")

    # 每季最好的 β：先看减氮 15% 时 5 种条件的平均产量损失，相差 < 0.2 个百分点的视为并列，再比所需氮量
    print("\n各季“最好的 β”分布（先比产量损失，并列时比所需氮量）")
    for crop, sel in (("水稻", lambda r: r["crop"] == "rice"), ("小麦·有冻土期", lambda r: r["crop"] == "wheat" and r["frost"]),
                      ("小麦·无冻土期", lambda r: r["crop"] == "wheat" and not r["frost"])):
        by = {}
        for r in rows:
            if sel(r):
                by.setdefault((r["zone"], r["season"]), []).append(r)
        best = []
        for rs in by.values():
            lo = min(r["N85_mean"] for r in rs)
            tie = [r for r in rs if r["N85_mean"] <= lo + 0.2]
            best.append(min(tie, key=lambda r: r["N_req"])["beta"])
        cnt = {b: best.count(b) for b in BETAS if best.count(b)}
        print(f"  {crop}：" + "，".join(f"β={b}: {c} 季" for b, c in cnt.items()) + f"；中位数 {st.median(best)}")


if __name__ == "__main__":
    main()
