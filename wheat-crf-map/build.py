"""小麦版：计算全部分区的推荐配方，输出 data.json，并生成地图页面 index.html。"""
import json
import math
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np

import model as M
from zones import ZONES, COUNTRIES

HERE = Path(__file__).parent


def rec_str(rec):
    return " + ".join(("尿素" if k == "urea" else f"CR{k}") + f" {p}%" for k, p in rec)


def best_recipe(ss, Ea=None, beta=None):
    if beta is not None:
        M.P["beta"] = beta
    best, best1, ideal, G, need, Fend, tau = M.optimize(ss, Ea=Ea)
    r2, t2 = M.round_recipe(best[1], best[2], G, need)
    r1, t1 = M.round_recipe(best1[1], best1[2], G, need)
    # 生产上组分越少越好：单一控释期的方案多用不超过 3% 的氮就选它
    if t1 <= t2 * 1.03:
        r2, t2 = r1, t1
    if beta is not None:
        M.P["beta"] = 1.3
    return dict(rec=r2, total=t2, simple=r1, simple_total=t1, ideal=ideal), G, need, Fend, tau


def split_urea(ss, mode):
    """对照：尿素分次施（同样要求冷年、常年、暖年都不断顿）。
    返青前施用的地区按 40/30/30（返青、拔节、旗叶）；其余按基肥 50% + 拔节 50%。"""
    n = ss["n"]
    if mode == "LW":
        fr, days = (0.4, 0.3, 0.3), [0, ss["tPI"], min(ss["tPI"] + 20, ss["tH"])]
    else:
        fr, days = (0.5, 0.5), [0, ss["tPI"]]
    base = M.urea_increment(n)
    worst = 0.0
    for dT in M.ROBUST_DT:
        _, need, _, _ = M.season_basis(ss, dT=dT)
        k = M.loss_rate(ss["Tpaddy"] + dT, M.P["k_loss"], ss["wet"])
        shape = np.zeros(n)
        for f, d0 in zip(fr, days):
            inc = np.zeros(n)
            inc[d0:] = base[: n - d0]
            shape += f * M.decay_conv(inc, k)
        worst = max(worst, float(np.max(need / np.maximum(shape, 1e-9))))
    return 5 * math.ceil(worst / 5), days


def run_season(args):
    zi, si = args
    z = ZONES[zi]
    s = z["seasons"][si]
    ss = M.season_setup(z, s)
    main, G, need, Fend, tau = best_recipe(ss)
    G0, need0, _, _ = M.season_basis(ss)
    sim = M.simulate(ss, main["rec"], main["total"], G0, need0, Fend)
    split_total, split_days = split_urea(ss, s["mode"])

    # 活化能敏感性：若实测 Ea 偏低/偏高，最优配方如何变化
    sens = {}
    for Ea in (38, 65):
        r, *_ = best_recipe(ss, Ea=Ea)
        sens[str(Ea)] = dict(rec=r["rec"], total=r["total"])
    # 形状敏感性：若做成 S 型（有滞后期）产品
    rS, *_ = best_recipe(ss, beta=2.5)
    # 年际温度波动：配方不变，释放随温度 ±1.5 °C 变化
    robust = {}
    for dT in (-1.5, 1.5):
        G2, need2, Fend2, _ = M.season_basis(ss, dT=dT)
        sup = sum(main["total"] * p / 100 * G2[k] for k, p in main["rec"])
        margin = sup - need2
        # 需求侧的“缓冲量”不算断顿：只看库是否低于 0
        buf = M.buffer_need(ss)
        pool = sup - (need2 - buf)
        robust[str(dT)] = dict(short_days=int(np.sum(pool < -0.5)),
                               min_buffer_pct=float(100 * np.min(margin / np.maximum(need2, 1e-6))),
                               unreleased_pct=float(100 * sum(main["total"] * p / 100 * (1 - Fend2[k])
                                                              for k, p in main["rec"] if k != "urea") / main["total"]))

    step = 2
    idx = list(range(0, ss["n"], step))
    if idx[-1] != ss["L"]:
        idx.append(ss["L"])
    tau_all = M.cum_tau(ss["Tpaddy"], M.P["Ea"], M.P["f_soil"])
    comps = {}
    for k, p in main["rec"]:
        kg = main["total"] * p / 100
        if k == "urea":
            cum = np.cumsum(M.urea_increment(ss["n"])) / (1 - M.P["urea_loss"])
        else:
            cum = M.release(tau_all, k)
        comps["尿素" if k == "urea" else f"CR{k}"] = [round(float(kg * cum[i]), 1) for i in idx]
    uptake_cum = np.cumsum(ss["d"])
    crop_cum = ss["Fcrop"] * ss["uptake_total"]
    c = M.CROP[s["type"]]
    n_cr = sum(p for k, p in main["rec"] if k != "urea") / 100
    kg_cr = main["total"] * n_cr / M.P["n_cru"]
    kg_u = main["total"] * (1 - n_cr) / M.P["n_urea"]
    return dict(
        zone=z["id"], idx=si, name=s["name"], est=s["mode"], est_label=M.EST[s["mode"]]["label"],
        type=s["type"], type_label=c["label"], note=s["note"], app=s["app"], days=ss["L"], sow_off=ss["sow_off"],
        cal_start=(ss["start"] + ss["sow_off"]) % 365, cal_days=ss["L"] - ss["sow_off"],
        start_doy=ss["start"], tPI=ss["tPI"], tH=ss["tH"], Y=s["Y"], INS=s["INS"], FN=s["FN"],
        nreq=c["nreq"], uptake_total=round(ss["uptake_total"]), fert_demand=round(ss["fert_demand"], 1),
        Tmean=round(float(ss["Tair"].mean()), 1), Tmin=round(float(ss["Tair"].min()), 1),
        Tmax=round(float(ss["Tair"].max()), 1), tau_end=round(float(tau_all[-1])),
        rec=main["rec"], rec_str=rec_str(main["rec"]), total=main["total"],
        simple=main["simple"], simple_str=rec_str(main["simple"]), simple_total=main["simple_total"],
        ideal=round(main["ideal"]) if main["ideal"] else None,
        kg_cr=round(kg_cr), kg_urea=round(kg_u), kg_product=round(kg_cr + kg_u),
        RE=round(sim["RE"], 3), unreleased=round(sim["unreleased_pct"], 1),
        split_total=split_total, split_days=split_days,
        sens={k: dict(rec_str=rec_str(v["rec"]), total=v["total"]) for k, v in sens.items()},
        stype=dict(rec_str=rec_str(rS["rec"]), total=rS["total"]),
        robust=robust,
        series=dict(
            t=idx,
            T=[round(float(ss["Tpaddy"][i]), 1) for i in idx],
            crop=[round(float(crop_cum[i]), 1) for i in idx],
            demand=[round(float(uptake_cum[i]), 1) for i in idx],
            comps=comps,
            pool=[round(float(sim["pool"][i]), 1) for i in idx],
            floor=[round(float(M.buffer_need(ss)[i]), 1) for i in idx],
        ),
    )


def main():
    jobs = [(zi, si) for zi, z in enumerate(ZONES) for si in range(len(z["seasons"]))]
    with Pool() as pool:
        out = pool.map(run_season, jobs)
    zones = []
    for zi, z in enumerate(ZONES):
        zs = [o for o in out if o["zone"] == z["id"]]
        zones.append(dict(id=z["id"], name=z["name"], country=z["country"], country_name=COUNTRIES[z["country"]],
                          station=z["station"], lat=z["lat"], lon=z["lon"], T=z["T"], crop=z["crop"], seasons=zs))
    params = dict(P=M.P, CROP=M.CROP, EST=M.EST, PRODUCTS=M.PRODUCTS)
    data = dict(zones=zones, params=params, countries=COUNTRIES)
    (HERE / "data.json").write_text(json.dumps(data, ensure_ascii=False, separators=(",", ":")))
    world = (HERE / "countries-50m.json").read_text()
    html = (HERE / "template.html").read_text()
    html = html.replace("/*__DATA__*/null", json.dumps(data, ensure_ascii=False, separators=(",", ":")))
    html = html.replace("/*__WORLD__*/null", world)
    (HERE / "index.html").write_text(html)
    # 汇总表
    print(f"{'分区':20s}{'季别':16s}{'配方':34s}{'总N':>5s}{'分次施尿素':>8s}{'当地常规':>6s}{'RE':>6s}")
    for z in zones:
        for s in z["seasons"]:
            print(f"{z['name'][:18]:20s}{s['name'][:14]:16s}{s['rec_str']:34s}{s['total']:5d}{s['split_total']:8d}{s['FN']:6d}{s['RE']:6.2f}")


if __name__ == "__main__":
    main()
