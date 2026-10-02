"""玉米版：计算全部分区的推荐配方，输出 data.json，并生成地图页面 index.html。"""
import json
import math
import sys
from multiprocessing import Pool
from pathlib import Path
from types import SimpleNamespace

import numpy as np

import maize_model as M
from maize_zones import ZONES, COUNTRIES, ENV_ZONE

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.parent / "common"))
import env as ENV  # noqa: E402
import field_release as FR  # noqa: E402
import scenarios as SC  # noqa: E402

ENV.ZONE_TABLE["maize"] = ENV_ZONE
STYPE_BETA = 2.0          # 玉米的 S 型情景取 β=2.0（β 扫描中玉米的最优值，见 ../analysis/beta_scan.py）
SPLIT = (0.4, 0.6)        # 对照：分次施尿素，基肥 40% + V10（大喇叭口前）追肥 60%


def adapter(ss):
    """S 型情景用到的模型接口（释放依赖本季的种植系统与水分，所以每季单独构造）。"""
    return SimpleNamespace(
        P=M.P,
        best_recipe=lambda ss_, Ea, beta: best_recipe(ss_, Ea=Ea, beta=beta),
        nominal_basis=lambda ss_: M.season_basis(ss_),
        season_basis=lambda ss_, dT: M.season_basis(ss_, dT=dT),
        simulate=M.simulate, buffer_need=M.buffer_need,
        cum_tau=lambda T, Ea, f: M.cum_tau(T, Ea, f, ss),
        release=M.release, urea_increment=M.urea_increment,
    )


def rec_str(rec):
    return " + ".join(("尿素" if k == "urea" else f"CR{k}") + f" {p}%" for k, p in rec)


def best_recipe(ss, Ea=None, beta=None):
    with SC.beta_set(M.P, M.P["beta"] if beta is None else beta):
        best, best1, ideal, G, need, Fend, tau = M.optimize(ss, Ea=Ea)
        r2, t2 = M.round_recipe(best[1], best[2], G, need)
        r1, t1 = M.round_recipe(best1[1], best1[2], G, need)
    pair, pair_total = r2, t2
    if t1 <= t2 * 1.03:            # 生产上组分越少越好：单一控释期多用不超过 3% 的氮就选它
        r2, t2 = r1, t1
    return dict(rec=r2, total=t2, simple=r1, simple_total=t1, ideal=ideal,
                pair=pair, pair_total=pair_total), G, need, Fend, tau


def split_urea(ss):
    """对照：尿素分两次施（同样要求冷年、常年、暖年都不断顿）。"""
    n = ss["n"]
    days = [0, ss["tPI"]]
    base = M.urea_increment(n)
    worst = 0.0
    for dT in M.ROBUST_DT:
        _, need, _, _ = M.season_basis(ss, dT=dT)
        k = M.loss_rate(M.coating_temp(ss, dT), M.P["k_loss"], ss["wet"])
        shape = np.zeros(n)
        for f, d0 in zip(SPLIT, days):
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
    split_total, split_days = split_urea(ss)

    sens = {}
    for Ea in (38, 65):
        r, *_ = best_recipe(ss, Ea=Ea)
        sens[str(Ea)] = dict(rec_str=rec_str(r["rec"]), total=r["total"])
    robust = {}
    buf = M.buffer_need(ss)
    for dT in (-1.5, 1.5):
        G2, need2, Fend2, _ = M.season_basis(ss, dT=dT)
        sup = sum(main["total"] * p / 100 * G2[k] for k, p in main["rec"])
        robust[str(dT)] = dict(
            short_days=int(np.sum(sup - (need2 - buf) < -0.5)),
            min_buffer_pct=float(100 * np.min((sup - need2) / np.maximum(need2, 1e-6))),
            unreleased_pct=float(100 * sum(main["total"] * p / 100 * (1 - Fend2[k])
                                           for k, p in main["rec"] if k != "urea") / main["total"]))

    idx = list(range(0, ss["n"], 2))
    if idx[-1] != ss["L"]:
        idx.append(ss["L"])
    tau_all = M.cum_tau(ss["Tpaddy"], M.P["Ea"], M.P["f_soil"], ss)
    alt = SC.s_type_alternatives(adapter(ss), ss, idx, STYPE_BETA, (38, 65), (-1.5, 1.5))
    comps = {}
    for k, p in main["rec"]:
        kg = main["total"] * p / 100
        cum = (np.cumsum(M.urea_increment(ss["n"])) / (1 - M.P["urea_loss"]) if k == "urea"
               else M.release(tau_all, k))
        comps["尿素" if k == "urea" else f"CR{k}"] = [round(float(kg * cum[i]), 1) for i in idx]
    c = M.CROP[s["type"]]
    n_cr = sum(p for k, p in main["rec"] if k != "urea") / 100
    kg_cr = main["total"] * n_cr / M.P["n_cru"]
    kg_u = main["total"] * (1 - n_cr) / M.P["n_urea"]
    est_label = f"{FR.SYSTEMS[s['system']]['label'].replace('玉米·', '')}，{FR.MOISTURE[s['moist']]['label'].split('（')[0]}"
    return dict(
        zone=z["id"], idx=si, name=s["name"], est=s["system"], est_label=est_label, moist=s["moist"],
        type=s["type"], type_label=c["label"], note=s["note"], app=s["app"], days=ss["L"], sow_off=0,
        cal_start=ss["start"], cal_days=ss["L"], start_doy=ss["start"], tPI=ss["tPI"], tH=ss["tH"],
        Y=s["Y"], INS=s["INS"], FN=s["FN"], nreq=c["nreq"], uptake_total=round(ss["uptake_total"]),
        fert_demand=round(ss["fert_demand"], 1),
        Tmean=round(float(ss["Tair"].mean()), 1), Tmin=round(float(ss["Tair"].min()), 1),
        Tmax=round(float(ss["Tair"].max()), 1), tau_end=round(float(tau_all[-1])),
        rec=main["rec"], rec_str=rec_str(main["rec"]), total=main["total"],
        simple=main["simple"], simple_str=rec_str(main["simple"]), simple_total=main["simple_total"],
        ideal=round(main["ideal"]) if main["ideal"] else None,
        kg_cr=round(kg_cr), kg_urea=round(kg_u), kg_product=round(kg_cr + kg_u),
        RE=round(sim["RE"], 3), unreleased=round(sim["unreleased_pct"], 1),
        split_total=split_total, split_days=split_days, sens=sens, alt=alt, robust=robust,
        series=dict(
            t=idx,
            T=[round(float(ss["Tpaddy"][i]), 1) for i in idx],
            crop=[round(float(ss["Fcrop"][i] * ss["uptake_total"]), 1) for i in idx],
            demand=[round(float(np.cumsum(ss["d"])[i]), 1) for i in idx],
            comps=comps,
            pool=[round(float(sim["pool"][i]), 1) for i in idx],
            floor=[round(float(buf[i]), 1) for i in idx],
        ),
    )


def main():
    jobs = [(zi, si) for zi, z in enumerate(ZONES) for si in range(len(z["seasons"]))]
    with Pool() as pool:
        out = pool.map(run_season, jobs)
    zones = []
    for z in ZONES:
        zs = [o for o in out if o["zone"] == z["id"]]
        zones.append(dict(id=z["id"], name=z["name"], country=z["country"], country_name=COUNTRIES[z["country"]],
                          station=z["station"], lat=z["lat"], lon=z["lon"], T=z["T"], crop=z["crop"], seasons=zs))
    params = dict(P=M.P, CROP=M.CROP, PRODUCTS=M.PRODUCTS)
    params["ENV"] = ENV.add_env("maize", zones)
    for key in ("s1", "s2"):
        ENV.add_env_alt("maize", zones, key)
    params["SCEN"] = dict(labels=SC.SCEN_LABELS, stype_beta=STYPE_BETA, lin_beta=M.P["beta"])
    data = dict(zones=zones, params=params, countries=COUNTRIES)
    blob = json.dumps(data, ensure_ascii=False, separators=(",", ":"))
    (HERE / "data.json").write_text(blob)
    html = (HERE / "template.html").read_text()
    html = html.replace("/*__DATA__*/null", blob)
    html = html.replace("/*__WORLD__*/null", (HERE / "countries-50m.json").read_text())
    common = HERE.parent / "common"
    html = html.replace("<!--__ENV_SECTION__-->", (common / "env_section.html").read_text())
    html = html.replace("/*__ENV_JS__*/", (common / "env.js").read_text())
    html = html.replace("/*__SCEN_JS__*/", (common / "scenario.js").read_text())
    (HERE / "index.html").write_text(html)
    print(f"{'产区':24s}{'配方':34s}{'总N':>5s}{'分次施尿素':>8s}{'当地常规':>6s}{'RE':>6s}  S型单档")
    for z in zones:
        for s in z["seasons"]:
            print(f"{z['name'][:22]:24s}{s['rec_str']:34s}{s['total']:5d}{s['split_total']:8d}{s['FN']:6d}"
                  f"{s['RE']:6.2f}  {s['alt']['s1']['rec_str'].replace('CR', 'S')} {s['alt']['s1']['total']}")


if __name__ == "__main__":
    main()
