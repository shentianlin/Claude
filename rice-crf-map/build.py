"""计算全部稻作区、稻季的推荐配方，输出 data.json、交互地图 index.html 和构建记录 BUILD_INFO.json。

用法
    python build.py            # 用全部 CPU 核并行计算
    python build.py --jobs 1   # 单进程（调试用）
"""
import argparse
import hashlib
import json
import os
import platform
import subprocess
import sys
from datetime import datetime, timezone
from multiprocessing import Pool
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import scipy

from crfmap.data import COUNTRIES, ZONES
from crfmap.demand import season_setup
from crfmap.optimize import best_recipe, buffer_need, season_basis, simulate, split_urea
from crfmap.params import CROP, EST, P, PARAMS_FILE, PRODUCTS, SCEN
from crfmap.release import cum_tau, release, urea_increment

HERE = Path(__file__).resolve().parent
COMMON = HERE.parent / "common"          # 环境效应、S 型情景模块（与小麦版共用）
sys.path.insert(0, str(COMMON))
import scenarios as SC                   # noqa: E402

# S 型情景用到的模型接口
ADAPTER = SimpleNamespace(
    P=P,
    # S 型滞后期对温度更敏感：S 型方案要求常年、偏冷、偏暖年份都不断顿（线性型方案本来就满足）
    best_recipe=lambda ss, Ea, beta: best_recipe(ss, Ea=Ea, beta=beta,
                                                 robust_dT=tuple(sorted(list(SCEN["robust_dT"]) + [0.0]))),
    nominal_basis=lambda ss: season_basis(ss),
    season_basis=lambda ss, dT: season_basis(ss, dT=dT),
    simulate=simulate, buffer_need=buffer_need, cum_tau=cum_tau, release=release, urea_increment=urea_increment,
)


def rec_str(rec):
    return " + ".join(("尿素" if k == "urea" else f"CR{k}") + f" {p}%" for k, p in rec)


def run_season(args):
    zi, si = args
    z = ZONES[zi]
    s = z["seasons"][si]
    ss = season_setup(z, s)
    main, G, need, Fend, tau = best_recipe(ss)
    sim = simulate(ss, main["rec"], main["total"], G, need, Fend)
    split_total, split_days = split_urea(ss, need)

    # 活化能敏感性：若实测 Ea 偏低/偏高，最优配方如何变化
    sens = {}
    for Ea in SCEN["ea_sensitivity"]:
        r, *_ = best_recipe(ss, Ea=Ea)
        sens[str(Ea)] = dict(rec=r["rec"], total=r["total"])
    # 年际温度波动：配方不变，释放随温度变化
    robust = {}
    for dT in SCEN["robust_dT"]:
        G2, need2, Fend2, _ = season_basis(ss, dT=dT)
        sup = sum(main["total"] * p / 100 * G2[k] for k, p in main["rec"])
        margin = sup - need2
        # 需求侧的“缓冲量”不算断顿：只看库是否低于 0
        buf = buffer_need(ss)
        pool = sup - (need2 - buf)
        robust[str(dT)] = dict(short_days=int(np.sum(pool < -0.5)),
                               min_buffer_pct=float(100 * np.min(margin / np.maximum(need2, 1e-6))),
                               unreleased_pct=float(100 * sum(main["total"] * p / 100 * (1 - Fend2[k])
                                                              for k, p in main["rec"] if k != "urea") / main["total"]))

    step = SCEN["series_step"]
    idx = list(range(0, ss["n"], step))
    if idx[-1] != ss["L"]:
        idx.append(ss["L"])
    tau_all = cum_tau(ss["Tpaddy"], P["Ea"], P["f_soil"])
    # 可选情景：S 型（有滞后期）产品的单档 / 两档方案
    alt = SC.s_type_alternatives(ADAPTER, ss, idx, SCEN["stype_beta"], SCEN["ea_sensitivity"], SCEN["robust_dT"])
    comps = {}
    for k, p in main["rec"]:
        kg = main["total"] * p / 100
        if k == "urea":
            cum = np.cumsum(urea_increment(ss["n"])) / (1 - P["urea_loss"])
        else:
            cum = release(tau_all, k)
        comps["尿素" if k == "urea" else f"CR{k}"] = [round(float(kg * cum[i]), 1) for i in idx]
    uptake_cum = np.cumsum(ss["d"])
    crop_cum = ss["Fcrop"] * ss["uptake_total"]
    c = CROP[s["type"]]
    n_cr = sum(p for k, p in main["rec"] if k != "urea") / 100
    kg_cr = main["total"] * n_cr / P["n_cru"]
    kg_u = main["total"] * (1 - n_cr) / P["n_urea"]
    return dict(
        zone=z["id"], idx=si, name=s["name"], est=s["est"], est_label=EST[s["est"]]["label"],
        type=s["type"], type_label=c["label"], note=s["note"], app=s["app"], days=s["days"],
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
        alt=alt,
        robust=robust,
        series=dict(
            t=idx,
            T=[round(float(ss["Tpaddy"][i]), 1) for i in idx],
            crop=[round(float(crop_cum[i]), 1) for i in idx],
            demand=[round(float(uptake_cum[i]), 1) for i in idx],
            comps=comps,
            pool=[round(float(sim["pool"][i]), 1) for i in idx],
            floor=[round(float(buffer_need(ss)[i]), 1) for i in idx],
        ),
    )


def _git(*args):
    try:
        return subprocess.run(["git", *args], cwd=HERE, capture_output=True, text=True, check=True).stdout.rstrip("\n")
    except (OSError, subprocess.CalledProcessError):
        return None


def build_info(outputs):
    """记录生成这次输出所用的源代码版本与运行环境，便于日后核对“某份 HTML 是哪个 commit 生成的”。"""
    generated = {"index.html", "data.json", "BUILD_INFO.json"}
    status = _git("status", "--porcelain", "--", ".", str(COMMON))
    dirty = None
    if status is not None:
        dirty = [line[3:] for line in status.splitlines()
                 if Path(line[3:].split(" -> ")[-1]).name not in generated and "__pycache__" not in line]
    return dict(
        source_commit=_git("rev-parse", "HEAD"),
        source_dirty_files=dirty,
        built_at_utc=datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
        params_file=PARAMS_FILE.name,
        outputs={name: hashlib.sha256(data).hexdigest() for name, data in outputs.items()},
    )


def main(argv=None):
    ap = argparse.ArgumentParser(description="生成稻季控释配方地图")
    ap.add_argument("--jobs", type=int, default=os.cpu_count(), help="并行进程数（默认 CPU 核数）")
    args = ap.parse_args(argv)

    jobs = [(zi, si) for zi, z in enumerate(ZONES) for si in range(len(z["seasons"]))]
    if args.jobs > 1:
        with Pool(args.jobs) as pool:
            out = pool.map(run_season, jobs)
    else:
        out = [run_season(j) for j in jobs]

    zones = []
    for z in ZONES:
        zs = [o for o in out if o["zone"] == z["id"]]
        zones.append(dict(id=z["id"], name=z["name"], country=z["country"], country_name=COUNTRIES[z["country"]],
                          station=z["station"], lat=z["lat"], lon=z["lon"], T=z["T"], crop=z["crop"], seasons=zs))
    params = dict(P=P, CROP=CROP, EST=EST, PRODUCTS=PRODUCTS)
    import env as ENV
    params["ENV"] = ENV.add_env("rice", zones)
    for key in ("s1", "s2"):
        ENV.add_env_alt("rice", zones, key)
    params["SCEN"] = dict(labels=SC.SCEN_LABELS, stype_beta=SCEN["stype_beta"], lin_beta=P["beta"])
    data = dict(zones=zones, params=params, countries=COUNTRIES)
    data_json = json.dumps(data, ensure_ascii=False, separators=(",", ":"))

    world = (HERE / "data" / "countries-50m.json").read_text()
    html = (HERE / "template.html").read_text()
    html = html.replace("/*__DATA__*/null", data_json)
    html = html.replace("/*__WORLD__*/null", world)
    html = html.replace("<!--__ENV_SECTION__-->", (COMMON / "env_section.html").read_text())
    html = html.replace("/*__ENV_JS__*/", (COMMON / "env.js").read_text())
    html = html.replace("/*__SCEN_JS__*/", (COMMON / "scenario.js").read_text())

    (HERE / "data.json").write_text(data_json)
    (HERE / "index.html").write_text(html)
    info = build_info({"index.html": html.encode(), "data.json": data_json.encode()})
    (HERE / "BUILD_INFO.json").write_text(json.dumps(info, ensure_ascii=False, indent=2) + "\n")

    print(f"{'分区':20s}{'季别':16s}{'配方':34s}{'总N':>5s}{'三次分施尿素':>8s}{'当地常规':>6s}{'RE':>6s}")
    for z in zones:
        for s in z["seasons"]:
            print(f"{z['name'][:18]:20s}{s['name'][:14]:16s}{s['rec_str']:34s}{s['total']:5d}{s['split_total']:8d}{s['FN']:6d}{s['RE']:6.2f}")
    print(f"\n{len(zones)} 个稻作区、{len(out)} 个稻季；源代码 commit {info['source_commit']}"
          + ("（有未提交改动：" + ", ".join(info["source_dirty_files"]) + "）" if info["source_dirty_files"] else ""))


if __name__ == "__main__":
    main()
