"""田间埋袋释放试验分析：用两套输入运行田间释放模型，计算误差，并做留一站点标定。

    运行 A（实测驱动）：埋袋深度的逐时土温 + 实测土壤水分（换算为基质势）+ 管理日志
    运行 B（气象驱动）：日气温 + 降水 + 管理日志，土温按种植系统预设估算，土壤水分用单层水量平衡估算
    B 与 A 的误差之差，就是“只用气象数据”的代价。

标定：用运行 A 的残差拟合 f_soil、水分敏感度 s、条施和表施系数、黏土和砂土系数（带先验，
没有对应数据的参数保持默认）。留一站点：每次去掉一个站点，用其余站点标定，再预测被去掉的站点。

用法
    python validation/analyze.py validation/example            # 数据目录内需有 7 个 CSV（见 README）
    python validation/analyze.py 数据目录 -o 结果目录
"""
import argparse
import csv
import datetime as dt
import logging
import math
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

sys.path.insert(0, str(Path(__file__).resolve().parent))
import vmodel as VM  # noqa: E402

FR = VM.FR

# 合格标准
RMSE_OK = 8.0           # 百分点
D80_ABS, D80_REL = 10.0, 0.15
B_OVER_A = 1.5
# 缺测处理
MIN_HOURS = 12          # 一天至少 12 个土温读数才算有效
MAX_GAP_MOIST = 21      # 土壤水分缺测 ≤ 21 天时线性插值（取袋时测一次的称重法也能用）
OBS_SD = 4.0            # 标定时假定的观测误差（百分点）

# 标定参数：名称、默认值、先验标准差、是否取对数、上下界
CAL = [
    ("f_soil", FR.DEFAULTS["f_soil"], 0.3, True, 0.3, 2.0),
    ("moisture_sens", FR.DEFAULTS["moisture_sens"], 0.15, False, 0.0, 0.9),
    ("band", FR.PLACEMENT["band"]["factor"], 0.2, True, 0.4, 1.5),
    ("surface", FR.PLACEMENT["surface"]["factor"], 0.5, True, 0.8, 6.0),
    ("clay", FR.TEXTURE["clay"]["factor"], 0.15, True, 0.6, 1.8),
    ("sand", FR.TEXTURE["sand"]["factor"], 0.15, True, 0.6, 1.8),
]
CAL_LABEL = {"f_soil": "土壤基础系数 f_soil", "moisture_sens": "水分敏感度 s", "band": "条施/侧深施系数",
             "surface": "表面撒施系数", "clay": "黏土系数", "sand": "砂土系数"}
DEFAULT_P = {k: v for k, v, *_ in CAL}

PLACEMENT_ALIAS = {"incorporated": "incorporated", "混施": "incorporated", "混施入土": "incorporated",
                   "撒施混土": "incorporated", "band": "band", "条施": "band", "侧深施": "band",
                   "surface": "surface", "表施": "surface", "表面撒施": "surface"}
TEXTURE_ALIAS = {"clay": "clay", "黏土": "clay", "粘土": "clay", "loam": "loam", "壤土": "loam",
                 "sand": "sand", "砂土": "sand", "沙土": "sand"}
EVENT_ALIAS = {"flood": "flood", "淹水": "flood", "上水": "flood", "复水": "flood",
               "drain": "drain", "排水": "drain", "落干": "drain", "晒田": "drain", "烤田": "drain",
               "irrigate": "irrigate", "灌溉": "irrigate", "film_on": "film_on", "覆膜": "film_on",
               "film_off": "film_off", "揭膜": "film_off", "harvest": "harvest", "收获": "harvest"}
PL_LABEL = {k: v["label"] for k, v in FR.PLACEMENT.items()}
PL_SHORT = {"incorporated": "混施", "band": "条施", "surface": "表施"}

# 各表的列：规范名 → 可用的列名（中英文都认）
COLS = {
    "sites": dict(site=["站点", "site_id", "site"], name=["名称", "name"], lat=["纬度", "lat"], lon=["经度", "lon"],
                  system=["种植系统", "system"], texture=["质地", "texture"],
                  apply=["施肥日期", "apply_date"], depth=["埋袋深度_cm", "depth_cm"],
                  theta_r=["theta_r"], theta_s=["theta_s"], alpha=["alpha_per_kPa", "alpha"], n_vg=["n_vg"]),
    "products": dict(product=["产品", "product"], lot=["批号", "lot"], D=["D_天", "D"], beta=["beta", "β"],
                     f0=["f0", "f₀"], Ea=["Ea_kJ", "Ea"], damage=["破损率", "damage"]),
    "bags": dict(site=["站点", "site_id", "site"], product=["产品", "product"], placement=["施肥位置", "placement"],
                 bag=["袋号", "bag_id"], date=["取样日期", "date"], n_init=["初始氮_g", "initial_N_g"],
                 n_rem=["剩余氮_g", "remaining_N_g"], rel=["累计释放%", "release_pct"]),
    "soil_temp": dict(site=["站点", "site_id", "site"], placement=["施肥位置", "placement"],
                      time=["时间", "datetime", "time"], T=["土温_C", "soil_temp_C", "T"]),
    "soil_moisture": dict(site=["站点", "site_id", "site"], time=["时间", "datetime", "time", "日期", "date"],
                          theta=["含水量_m3m3", "theta"], psi=["基质势_kPa", "psi_kPa"]),
    "weather": dict(site=["站点", "site_id", "site"], date=["日期", "date"], tmax=["最高气温_C", "tmax"],
                    tmin=["最低气温_C", "tmin"], tmean=["平均气温_C", "tmean"], precip=["降水_mm", "precip_mm"]),
    "management": dict(site=["站点", "site_id", "site"], date=["日期", "date"], event=["事件", "event"],
                       amount=["水量_mm", "amount_mm"]),
}
REQUIRED = {"sites": ["site", "system", "texture", "apply", "lat"], "products": ["product", "D"],
            "bags": ["site", "product", "placement", "date"], "soil_temp": ["site", "time", "T"],
            "soil_moisture": ["site", "time"], "weather": ["site", "date"], "management": ["site", "date", "event"]}
OPTIONAL_FILES = {"soil_temp", "soil_moisture", "management"}


class DataError(Exception):
    pass


# ---------------------------------------------------------------- 读取
def parse_time(s):
    s = s.strip().replace("/", "-").replace("T", " ")
    date, _, clock = s.partition(" ")
    y, m, d = (int(x) for x in date.split("-"))
    if not clock:
        return dt.datetime(y, m, d)
    parts = [int(float(x)) for x in clock.split(":")] + [0, 0]
    return dt.datetime(y, m, d, parts[0], parts[1])


def num(s):
    s = (s or "").strip().replace("%", "")
    return float(s) if s else None


def read_table(folder, kind):
    path = folder / f"{kind}.csv"
    if not path.exists():
        if kind in OPTIONAL_FILES:
            return []
        raise DataError(f"缺少 {path.name}")
    with open(path, encoding="utf-8-sig", newline="") as fh:
        reader = csv.reader(fh)
        header = [h.strip() for h in next(reader)]
        index = {}
        for key, names in COLS[kind].items():
            for nm in names:
                if nm in header:
                    index[key] = header.index(nm)
                    break
        miss = [COLS[kind][k][0] for k in REQUIRED[kind] if k not in index]
        if miss:
            raise DataError(f"{path.name} 缺少列：{'、'.join(miss)}")
        rows = []
        for ln, row in enumerate(reader, start=2):
            if not any(c.strip() for c in row) or row[0].strip().startswith("#"):
                continue
            rec = {k: (row[i].strip() if i < len(row) else "") for k, i in index.items()}
            rec["_line"] = f"{path.name} 第 {ln} 行"
            rows.append(rec)
    return rows


def load(folder):
    T = {k: read_table(folder, k) for k in COLS}
    sites = {}
    for r in T["sites"]:
        sys_, tex = r["system"], TEXTURE_ALIAS.get(r["texture"])
        if sys_ not in FR.SYSTEMS:
            raise DataError(f"{r['_line']}：种植系统“{sys_}”不认识，可选 {', '.join(FR.SYSTEMS)}")
        if tex is None:
            raise DataError(f"{r['_line']}：质地“{r['texture']}”不认识，可选 黏土/壤土/砂土")
        vg = VM.VG_DEFAULT[tex]
        if all(num(r.get(k)) is not None for k in ("theta_r", "theta_s", "alpha", "n_vg")):
            vg = tuple(num(r[k]) for k in ("theta_r", "theta_s", "alpha", "n_vg"))
        sites[r["site"]] = dict(id=r["site"], name=r.get("name") or r["site"], lat=num(r["lat"]),
                                system=sys_, texture=tex, apply=parse_time(r["apply"]).date(), vg=vg,
                                vg_measured=vg is not VM.VG_DEFAULT[tex])
    products = {}
    for r in T["products"]:
        def frac(v, dflt):
            v = num(v)
            return dflt if v is None else (v / 100 if v > 1 else v)
        products[r["product"]] = dict(
            D=num(r["D"]), beta=num(r.get("beta")) or FR.DEFAULTS["beta"],
            f0=frac(r.get("f0"), FR.DEFAULTS["f0"]), Ea=num(r.get("Ea")) or FR.DEFAULTS["Ea"],
            damage=frac(r.get("damage"), FR.DEFAULTS["damage"]), lot=r.get("lot", ""))

    # 取袋数据
    rel_vals = [num(r.get("rel")) for r in T["bags"] if num(r.get("rel")) is not None]
    rel_pct = bool(rel_vals) and max(rel_vals) > 1.5
    obs = defaultdict(lambda: defaultdict(list))
    for r in T["bags"]:
        s = sites.get(r["site"])
        if s is None:
            raise DataError(f"{r['_line']}：站点 {r['site']} 不在 sites.csv 中")
        if r["product"] not in products:
            raise DataError(f"{r['_line']}：产品 {r['product']} 不在 products.csv 中")
        pl = PLACEMENT_ALIAS.get(r["placement"])
        if pl is None:
            raise DataError(f"{r['_line']}：施肥位置“{r['placement']}”不认识，可选 混施/条施/侧深施/表施")
        ni, nr, rel = num(r.get("n_init")), num(r.get("n_rem")), num(r.get("rel"))
        if ni and nr is not None:
            f = 1 - nr / ni
        elif rel is not None:
            f = rel / 100 if rel_pct else rel
        else:
            raise DataError(f"{r['_line']}：需要“初始氮_g + 剩余氮_g”或“累计释放%”")
        day = (parse_time(r["date"]).date() - s["apply"]).days
        if day < 0:
            raise DataError(f"{r['_line']}：取样日期早于施肥日期")
        obs[(r["site"], r["product"], pl)][day].append(f)

    weather = defaultdict(dict)
    for r in T["weather"]:
        tmax, tmin, tmean = num(r.get("tmax")), num(r.get("tmin")), num(r.get("tmean"))
        if tmean is None and tmax is not None and tmin is not None:
            tmean = (tmax + tmin) / 2
        weather[r["site"]][parse_time(r["date"]).date()] = dict(
            tmax=tmax, tmin=tmin, tmean=tmean, precip=num(r.get("precip")))
    soilT = defaultdict(lambda: defaultdict(list))
    for r in T["soil_temp"]:
        pl = PLACEMENT_ALIAS.get(r.get("placement", ""), "")
        v = num(r["T"])
        if v is not None:
            soilT[(r["site"], pl)][parse_time(r["time"]).date()].append(v)
    moist = defaultdict(lambda: defaultdict(list))
    for r in T["soil_moisture"]:
        s = sites.get(r["site"])
        if s is None:
            continue
        psi, th = num(r.get("psi")), num(r.get("theta"))
        if psi is None and th is not None:
            psi = float(VM.vg_psi(th / 100 if th > 1 else th, *s["vg"]))
        if psi is not None:
            moist[r["site"]][parse_time(r["time"]).date()].append(-abs(psi))
    mgmt = defaultdict(list)
    for r in T["management"]:
        ev = EVENT_ALIAS.get(r["event"])
        if ev is None:
            raise DataError(f"{r['_line']}：事件“{r['event']}”不认识，可选 淹水/排水/灌溉/覆膜/揭膜/收获")
        mgmt[r["site"]].append((parse_time(r["date"]).date(), ev, num(r.get("amount")) or 0.0))
    return sites, products, obs, weather, soilT, moist, mgmt


# ---------------------------------------------------------------- 逐日环境
def fill(values):
    """含 None 的序列线性插值，两端外推取最近值；全空返回 None。"""
    v = np.array([np.nan if x is None else x for x in values], float)
    ok = ~np.isnan(v)
    if not ok.any():
        return None
    idx = np.arange(len(v))
    return np.interp(idx, idx[ok], v[ok])


def site_environment(s, n, weather, moist, mgmt, notes):
    dates = [s["apply"] + dt.timedelta(days=k) for k in range(n)]
    w = weather.get(s["id"], {})
    tmean = fill([w.get(d, {}).get("tmean") for d in dates])
    if tmean is None:
        raise DataError(f"站点 {s['id']} 没有气温数据（weather.csv）")
    tmax = fill([w.get(d, {}).get("tmax") for d in dates])
    tmin = fill([w.get(d, {}).get("tmin") for d in dates])
    pr = [w.get(d, {}).get("precip") for d in dates]
    n_missing_w = sum(d not in w for d in dates)
    if n_missing_w:
        notes.append(f"{s['id']}：气象数据缺 {n_missing_w} 天，气温已插值、降水按 0 计")
    precip = np.array([x or 0.0 for x in pr])

    # 管理事件
    sysd = FR.SYSTEMS[s["system"]]
    ev = sorted(mgmt.get(s["id"], []))
    has_water_events = any(e in ("flood", "drain") for _, e, _ in ev)
    if has_water_events:
        # 第一个水分事件之前的状态：第一个事件是排水，说明之前在淹水
        first = next(e for _, e, _ in ev if e in ("flood", "drain"))
        state = first == "drain"
    else:
        state = sysd["moisture"] == "flooded"
    film_state = s["system"] == "maize_film"
    flooded, film, irrig = np.zeros(n, bool), np.zeros(n, bool), np.zeros(n)
    film_events = any(e in ("film_on", "film_off") for _, e, _ in ev)
    for k, d in enumerate(dates):
        for dd, e, amt in ev:
            if dd == d:
                if e == "flood":
                    state = True
                elif e == "drain":
                    state = False
                elif e == "film_on":
                    film_state = True
                elif e == "film_off":
                    film_state = False
                elif e == "irrigate":
                    irrig[k] += amt
        flooded[k] = state
        film[k] = film_state if film_events else (s["system"] == "maize_film" and k < 60)
    if "switch" in sysd and not has_water_events:
        flooded[sysd["switch"][0]:] = True

    # 运行 B 的基质势
    kc = VM.SYSTEM_KC[s["system"]]
    if tmax is not None and tmin is not None and not (sysd["moisture"] == "awd" and not has_water_events):
        doy = np.array([d.timetuple().tm_yday for d in dates])
        et0 = VM.et0_hargreaves(tmax, tmin, s["lat"], doy)
        psiB = VM.bucket(precip, irrig, et0, tmean, kc, s["vg"], flooded, film)
        modeB = "水量平衡"
    else:
        psiB = FR.matric_potential(sysd["moisture"], n)
        if "switch" in sysd:
            psiB[sysd["switch"][0]:] = FR.matric_potential(sysd["switch"][1], n - sysd["switch"][0])
        psiB[flooded] = 0.0
        modeB = "种植系统预设"

    # 运行 A 的基质势：实测（缺测 ≤ MAX_GAP_MOIST 天插值），淹水日取 0，其余用 B
    m = moist.get(s["id"], {})
    psiA, n_meas = psiB.copy(), 0
    meas_days = [k for k, d in enumerate(dates) if d in m]
    if meas_days:
        vals = {k: float(np.mean(m[dates[k]])) for k in meas_days}
        for a, b in zip(meas_days, meas_days[1:] + [None]):
            psiA[a] = vals[a]
            if b is not None and b - a <= MAX_GAP_MOIST:
                lg = np.interp(np.arange(a, b + 1), [a, b], np.log10(-np.array([vals[a], vals[b]]) + 1e-3))
                psiA[a:b + 1] = -(10 ** lg)
        n_meas = len(meas_days)
    psiA[flooded] = 0.0
    return dict(dates=dates, tmean=tmean, precip=precip, flooded=flooded, film=film, psiA=psiA, psiB=psiB,
                modeB=modeB, n_moist=n_meas, has_moist=bool(meas_days))


def build_groups(sites, products, obs, weather, soilT, moist, mgmt, notes):
    by_site = defaultdict(list)
    for key in obs:
        by_site[key[0]].append(key)
    groups, envs = [], {}
    for sid, keys in by_site.items():
        s = sites[sid]
        n = max(max(obs[k]) for k in keys) + 1
        env = site_environment(s, n, weather, moist, mgmt, notes)
        envs[sid] = env
        for key in sorted(keys):
            _, prod, pl = key
            p = products[prod]
            aT_B, Ts_B = VM.air_to_factor(env["tmean"], s["system"], pl, p["Ea"])
            logger = soilT.get((sid, pl)) or {}
            shift = 0.0
            if not logger and soilT.get((sid, "")):
                logger = soilT[(sid, "")]
                if pl == "surface":
                    shift = FR.PLACEMENT["surface"]["dT"]
                    notes.append(f"{sid}：表施处理没有单独的地表温度记录，借用埋深土温 +{shift:g} °C，误差会偏大")
            aT_A, TsA, n_logger = aT_B.copy(), np.full(n, np.nan), 0
            for k, d in enumerate(env["dates"]):
                hrs = logger.get(d, [])
                if len(hrs) >= MIN_HOURS:
                    aT_A[k] = VM.day_factor_hourly(np.array(hrs) + shift, p["Ea"])
                    TsA[k] = float(np.mean(hrs)) + shift
                    n_logger += 1
            days = np.array(sorted(obs[key]))
            vals = [np.array(obs[key][d]) for d in days]
            groups.append(dict(
                key=key, site=sid, product=prod, placement=pl, texture=s["texture"], prod=p, n=n,
                aT={"A": aT_A, "B": aT_B}, x={"A": VM.dryness(env["psiA"]), "B": VM.dryness(env["psiB"])},
                TsA=TsA, TsB=Ts_B, n_logger=n_logger,
                days=days, mean=np.array([v.mean() for v in vals]),
                sd=np.array([v.std(ddof=1) if len(v) > 1 else np.nan for v in vals]),
                nbag=np.array([len(v) for v in vals])))
            if days[0] == 0:
                f0e = p["f0"] + (1 - p["f0"]) * p["damage"]
                d0 = float(np.mean(vals[0]))
                if abs(d0 - f0e) > 0.05:
                    notes.append(f"{sid} {prod} {PL_SHORT[pl]}：第 0 天对照袋释放 {d0 * 100:.1f}%，与初期溶出 + 破损率 "
                                 f"{f0e * 100:.1f}% 相差超过 5 个百分点，请检查制袋、冲洗操作或破损率")
            if n_logger < n:
                miss = n - n_logger
                notes.append(f"{sid} {prod} {PL_SHORT[pl]}：{miss} 天没有有效土温记录，运行 A 这些天用气温估算"
                             + ("（整季都没有，运行 A 与 B 的温度相同）" if n_logger == 0 else ""))
    return groups, envs


# ---------------------------------------------------------------- 模型与误差
def predict(g, run, P):
    pl_f = 1.0 if g["placement"] == "incorporated" else P[g["placement"]]
    tex_f = 1.0 if g["texture"] == "loam" else P[g["texture"]]
    rate = g["aT"][run] * (1 - P["moisture_sens"] * g["x"][run]) * pl_f * tex_f * P["f_soil"]
    p = g["prod"]
    return VM.release_curve(rate, p["D"], p["beta"], p["f0"], p["damage"])


def metrics(g, F):
    m = g["days"] > 0
    pred = F[g["days"]]
    err = (pred - g["mean"])[m] * 100
    rmse = float(np.sqrt(np.mean(err ** 2))) if m.any() else float("nan")
    bias = float(np.mean(err)) if m.any() else float("nan")
    out = dict(rmse=rmse, bias=bias)
    for p in (0.5, 0.8):
        o = VM.day_at(g["mean"], p, g["days"])
        e = VM.day_at(F, p)
        out[f"d{int(p * 100)}_obs"], out[f"d{int(p * 100)}_pred"] = o, e
    o, e = out["d80_obs"], out["d80_pred"]
    if o is None:
        out["d80_ok"] = None
    else:
        out["d80_ok"] = e is not None and abs(e - o) <= max(D80_ABS, D80_REL * o)
    out["rmse_ok"] = rmse <= RMSE_OK
    out["pass"] = out["rmse_ok"] and out["d80_ok"] is not False
    return out


# ---------------------------------------------------------------- 标定
def support(groups, with_sites=False):
    """每个参数有多少次取样数据可用来标定（with_sites=True 时同时返回来自哪些站点）。"""
    sup = dict(f_soil=0, moisture_sens=0, band=0, surface=0, clay=0, sand=0)
    src = {k: set() for k in sup}
    for g in groups:
        k = int((g["days"] > 0).sum())
        names = ["f_soil"] + [x for x in (g["placement"], g["texture"]) if x in sup]
        if g["x"]["A"][: g["days"].max()].mean() > 0.05:
            names.append("moisture_sens")
        for nm in names:
            sup[nm] += k
            src[nm].add(g["site"])
    return (sup, src) if with_sites else sup


def calibrate(groups):
    sup = support(groups)
    free = [c for c in CAL if sup[c[0]] > 0]

    def unpack(z):
        P = dict(DEFAULT_P)
        for (name, *_r, is_log, _lo, _hi), v in zip(free, z):
            P[name] = math.exp(v) if is_log else v
        return P

    def resid(z):
        P = unpack(z)
        r = []
        for g in groups:
            F = predict(g, "A", P)
            m = g["days"] > 0
            r.append((F[g["days"]] - g["mean"])[m] * 100 / OBS_SD)
        for (name, d0, sd, is_log, *_), v in zip(free, z):
            r.append([(v - math.log(d0)) / sd if is_log else (v - d0) / sd])
        return np.concatenate(r)

    z0 = [math.log(d) if lg else d for _, d, _, lg, *_ in free]
    lo = [math.log(a) if lg else a for _, _, _, lg, a, _ in free]
    hi = [math.log(b) if lg else b for _, _, _, lg, _, b in free]
    if not free:
        return dict(DEFAULT_P), sup
    sol = least_squares(resid, z0, bounds=(lo, hi))
    return unpack(sol.x), sup


# ---------------------------------------------------------------- 主流程
def analyse(folder, out):
    notes = []
    sites, products, obs, weather, soilT, moist, mgmt = load(folder)
    if not obs:
        raise DataError("bags.csv 中没有数据")
    groups, envs = build_groups(sites, products, obs, weather, soilT, moist, mgmt, notes)
    site_ids = sorted({g["site"] for g in groups})

    P_all, sup_all = calibrate(groups)
    folds = {}
    if len(site_ids) > 1:
        for sid in site_ids:
            folds[sid] = calibrate([g for g in groups if g["site"] != sid])
    else:
        notes.append("只有一个站点，无法做留一站点验证")

    sets = [("默认", lambda g: DEFAULT_P)]
    if folds:
        sets.append(("留一站点标定", lambda g: folds[g["site"]][0]))
    sets.append(("全部数据标定", lambda g: P_all))

    rows, preds = [], {}
    for g in groups:
        for sname, getP in sets:
            for run in ("A", "B"):
                F = predict(g, run, getP(g))
                preds[(g["key"], sname, run)] = F
                rows.append(dict(g=g, set=sname, run=run, **metrics(g, F)))

    out.mkdir(parents=True, exist_ok=True)
    write_csvs(out, groups, envs, sites, rows, preds, sets, P_all, sup_all, folds, site_ids)
    summary = summarise(rows, site_ids, [s for s, _ in sets])
    diag = diagnostics(groups, envs, sites)
    plot_sites(out, groups, sites, preds, sets)
    plot_scatter(out, groups, preds, [s for s, _ in sets])
    report(out, folder, sites, groups, summary, diag, P_all, sup_all, folds, notes, rows)
    print_summary(summary, P_all, sup_all, notes, out)
    return 0


def fmt(v, nd=1):
    return "" if v is None or (isinstance(v, float) and math.isnan(v)) else f"{v:.{nd}f}"


def write_csvs(out, groups, envs, sites, rows, preds, sets, P_all, sup_all, folds, site_ids):
    with open(out / "metrics.csv", "w", encoding="utf-8-sig", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["站点", "产品", "施肥位置", "参数", "运行", "取样次数", "RMSE_百分点", "偏差_百分点",
                    "d50_实测", "d50_预测", "d80_实测", "d80_预测", "d80_误差_天", "RMSE合格", "d80合格", "合格"])
        for r in rows:
            g = r["g"]
            e = (r["d80_pred"] - r["d80_obs"]) if r["d80_obs"] is not None and r["d80_pred"] is not None else None
            w.writerow([g["site"], g["product"], PL_SHORT[g["placement"]], r["set"], r["run"],
                        int((g["days"] > 0).sum()), fmt(r["rmse"]), fmt(r["bias"]), fmt(r["d50_obs"]),
                        fmt(r["d50_pred"]), fmt(r["d80_obs"]), fmt(r["d80_pred"]), fmt(e),
                        "是" if r["rmse_ok"] else "否",
                        "—" if r["d80_ok"] is None else ("是" if r["d80_ok"] else "否"),
                        "是" if r["pass"] else "否"])
    snames = [s for s, _ in sets]
    with open(out / "observed_vs_predicted.csv", "w", encoding="utf-8-sig", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["站点", "产品", "施肥位置", "天数", "日期", "实测均值%", "标准差%", "袋数"]
                   + [f"{s}_{r}%" for s in snames for r in ("A", "B")])
        for g in groups:
            ap = sites[g["site"]]["apply"]
            for i, d in enumerate(g["days"]):
                w.writerow([g["site"], g["product"], PL_SHORT[g["placement"]], int(d),
                            (ap + dt.timedelta(days=int(d))).isoformat(), fmt(g["mean"][i] * 100),
                            fmt(g["sd"][i] * 100), int(g["nbag"][i])]
                           + [fmt(preds[(g["key"], s, r)][d] * 100) for s in snames for r in ("A", "B")])
    with open(out / "calibration.csv", "w", encoding="utf-8-sig", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["参数", "说明", "默认值", "全部数据标定", "数据量（取样次数）"]
                   + [f"去掉{s}" for s in site_ids if s in folds])
        for name, d0, *_ in CAL:
            w.writerow([name, CAL_LABEL[name], fmt(d0, 3),
                        fmt(P_all[name], 3) if sup_all[name] else "无数据，保持默认", sup_all[name]]
                       + [fmt(folds[s][0][name], 3) if folds[s][1][name] else "默认" for s in site_ids if s in folds])
    with open(out / "daily_inputs.csv", "w", encoding="utf-8-sig", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["站点", "施肥位置", "天数", "日期", "气温_C", "实测土温_C", "估算土温_C",
                    "基质势A_kPa", "基质势B_kPa", "淹水", "覆膜"])
        seen = set()
        for g in groups:
            k = (g["site"], g["placement"])
            if k in seen:
                continue
            seen.add(k)
            env = envs[g["site"]]
            for d in range(g["n"]):
                w.writerow([g["site"], PL_SHORT[g["placement"]], d, env["dates"][d].isoformat(),
                            fmt(env["tmean"][d]), fmt(g["TsA"][d]), fmt(g["TsB"][d]),
                            fmt(env["psiA"][d]), fmt(env["psiB"][d]),
                            int(env["flooded"][d]), int(env["film"][d])])


def summarise(rows, site_ids, snames):
    S = {}
    for sname in snames:
        for run in ("A", "B"):
            rs = [r for r in rows if r["set"] == sname and r["run"] == run]
            errs = []
            for r in rs:
                g = r["g"]
                errs.append(r["rmse"] ** 2 * int((g["days"] > 0).sum()))
            nobs = sum(int((r["g"]["days"] > 0).sum()) for r in rs)
            d80 = [abs(r["d80_pred"] - r["d80_obs"]) for r in rs
                   if r["d80_obs"] is not None and r["d80_pred"] is not None]
            S[(sname, run)] = dict(n=len(rs), rmse=math.sqrt(sum(errs) / max(nobs, 1)),
                                   d80=float(np.median(d80)) if d80 else None,
                                   d80_missed=sum(1 for r in rs if r["d80_obs"] is not None and r["d80_pred"] is None),
                                   passed=sum(r["pass"] for r in rs))
            for sid in site_ids:
                rr = [r for r in rs if r["g"]["site"] == sid]
                k = sum(int((r["g"]["days"] > 0).sum()) for r in rr)
                S[(sname, run, sid)] = math.sqrt(sum(r["rmse"] ** 2 * int((r["g"]["days"] > 0).sum())
                                                     for r in rr) / max(k, 1))
    return S


def diagnostics(groups, envs, sites):
    """实测土温 − 估算土温，用来改进运行 B 的土温预设。"""
    D = []
    seen = set()
    for g in groups:
        k = (g["site"], g["placement"])
        if k in seen or not np.isfinite(g["TsA"]).any():
            continue
        seen.add(k)
        ok = np.isfinite(g["TsA"])
        env = envs[g["site"]]
        D.append(dict(site=g["site"], placement=g["placement"], days=int(ok.sum()),
                      air=float(env["tmean"][ok].mean()), meas=float(g["TsA"][ok].mean()),
                      est=float(g["TsB"][ok].mean()),
                      diff=float((g["TsA"][ok] - g["TsB"][ok]).mean()),
                      sd=float((g["TsA"][ok] - g["TsB"][ok]).std())))
    return D


# ---------------------------------------------------------------- 输出
def report(out, folder, sites, groups, S, diag, P_all, sup_all, folds, notes, rows):
    snames = [s for s in ("默认", "留一站点标定", "全部数据标定") if (s, "A") in S]
    L = ["# 田间埋袋释放试验分析报告", "",
         f"数据目录：`{folder}`；生成时间：{dt.datetime.now():%Y-%m-%d %H:%M}", "",
         f"共 {len(sites)} 个站点、{len(groups)} 个“站点 × 产品 × 施肥位置”组合。"
         "运行 A 用实测土温和土壤水分，运行 B 只用气温和降水。", "",
         "合格标准：RMSE ≤ 8 个百分点；d80（田间释放 80% 的天数）误差 ≤ ±10 天或 ±15%（取大者）；"
         "运行 B 的 RMSE 不超过运行 A 的 1.5 倍。", "",
         "## 1. 总体误差", "",
         "| 参数 | 运行 | 组合数 | 合并 RMSE（百分点） | d80 误差中位数（天） | 预测未达 80% | 合格组合 |",
         "|---|---|---|---|---|---|---|"]
    for sname in snames:
        for run in ("A", "B"):
            s = S[(sname, run)]
            L.append(f"| {sname} | {run} | {s['n']} | {s['rmse']:.1f} | {fmt(s['d80'])} | "
                     f"{s['d80_missed']} | {s['passed']}/{s['n']} |")
    L += ["", "“留一站点标定”是关键结果：每个站点的预测都只用了其他站点标定的参数，代表模型推广到新地区时的表现。"
          "“全部数据标定”会偏乐观，只用来给出最终参数。", "",
          "## 2. 各站点：只用气象数据的代价", "",
          "| 站点 | 名称 | 参数 | RMSE A | RMSE B | B/A | B/A ≤ 1.5 |", "|---|---|---|---|---|---|---|"]
    site_ids = sorted({g["site"] for g in groups})
    for sid in site_ids:
        for sname in snames[:2]:
            a, b = S[(sname, "A", sid)], S[(sname, "B", sid)]
            ratio = b / a if a > 0 else float("inf")
            L.append(f"| {sid} | {sites[sid]['name']} | {sname} | {a:.1f} | {b:.1f} | {ratio:.2f} | "
                     f"{'是' if ratio <= B_OVER_A or b <= RMSE_OK else '否'} |")
    L += ["", "B 的 RMSE 本身已 ≤ 8 个百分点时，即使 B/A 大于 1.5 也视为可接受。", "",
          "## 3. 标定参数", "",
          "| 参数 | 默认 | 全部数据标定 | 留一站点范围 | 数据量（取样次数） | 数据来源站点 |",
          "|---|---|---|---|---|---|"]
    _, src = support(groups, with_sites=True)
    single = []
    for name, d0, *_ in CAL:
        vals = [f[0][name] for f in folds.values() if f[1][name]]
        rng = f"{min(vals):.2f}–{max(vals):.2f}" if vals else "—"
        L.append(f"| {CAL_LABEL[name]} | {d0:.2f} | "
                 f"{f'{P_all[name]:.3f}' if sup_all[name] else '无数据，保持默认'} | {rng} | {sup_all[name]} | "
                 f"{'、'.join(sorted(src[name])) or '—'} |")
        if len(src[name]) == 1:
            single.append(CAL_LABEL[name])
    L += ["", "去掉某个站点后，如果某个参数在其余站点没有数据，这一折就保持默认值，所以该站点的留一预测会偏差较大；"
          "这正说明推广到没有同类数据的地区时的风险。"]
    if single:
        L += ["", f"**注意**：{'、'.join(single)} 只有 1 个站点的数据，标定值会和该站点的其他条件"
              "（气候、水分、土温估算误差）混在一起，不能单独解释。要分开这些影响，需要至少 2 个站点、"
              "或在同一站点内设对照处理（如同一地块混施与表施并排）。"]
    L += ["", "标定值稳定后，可写入 `common/field_release.py` 的 DEFAULTS、PLACEMENT、TEXTURE。", "",
          "## 4. 土温：实测与估算", "",
          "| 站点 | 施肥位置 | 有效天数 | 平均气温 | 实测土温 | 估算土温 | 实测 − 估算 | 逐日差标准差 |",
          "|---|---|---|---|---|---|---|---|"]
    for d in diag:
        L.append(f"| {d['site']} | {PL_SHORT[d['placement']]} | {d['days']} | {d['air']:.1f} | {d['meas']:.1f} | "
                 f"{d['est']:.1f} | {d['diff']:+.1f} | {d['sd']:.1f} |")
    L += ["", "“实测 − 估算”系统偏离超过约 1 °C 时，应修正 `field_release.py` 中对应种植系统的土温偏移（offset）。", "",
          "## 5. 各组合明细（留一站点标定；只有一个站点时为默认参数）", "",
          "| 站点 | 产品 | 位置 | 运行 | RMSE | 偏差 | d80 实测 | d80 预测 | 合格 |", "|---|---|---|---|---|---|---|---|---|"]
    show = snames[1] if len(snames) == 3 else snames[0]
    for r in rows:
        if r["set"] != show:
            continue
        g = r["g"]
        L.append(f"| {g['site']} | {g['product']} | {PL_SHORT[g['placement']]} | {r['run']} | {r['rmse']:.1f} | "
                 f"{r['bias']:+.1f} | {fmt(r['d80_obs'], 0) or '未达'} | {fmt(r['d80_pred'], 0) or '未达'} | "
                 f"{'是' if r['pass'] else '否'} |")
    if notes:
        L += ["", "## 6. 数据提示", ""] + [f"- {n}" for n in dict.fromkeys(notes)]
    L += ["", "## 输出文件", "",
          "- `metrics.csv`：每个组合 × 参数 × 运行的误差和是否合格",
          "- `observed_vs_predicted.csv`：每次取样的实测值与各情景预测值",
          "- `calibration.csv`：标定参数（含每一折留一站点的值）",
          "- `daily_inputs.csv`：逐日气温、实测与估算土温、基质势、淹水和覆膜状态，可用于检查输入",
          "- `fig_站点.png`：各站点释放曲线；`fig_scatter.png`：预测值与实测值散点图"]
    (out / "report.md").write_text("\n".join(L) + "\n", encoding="utf-8")


def print_summary(S, P_all, sup_all, notes, out):
    print("\n参数            运行  组合  合并RMSE  d80误差中位数  合格")
    for (k, v) in S.items():
        if len(k) == 2:
            print(f"{k[0]:<12}  {k[1]:>3}  {v['n']:>4}  {v['rmse']:>8.1f}  {fmt(v['d80']):>12}  "
                  f"{v['passed']}/{v['n']}")
    print("\n全部数据标定参数：" + "，".join(
        f"{k}={P_all[k]:.3f}" if sup_all[k] else f"{k}=默认" for k in DEFAULT_P))
    if notes:
        print(f"\n数据提示 {len(set(notes))} 条，见 report.md")
    print(f"\n结果已写入 {out}")


def _style():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
    plt.rcParams["font.family"] = ["Microsoft YaHei", "PingFang SC", "Noto Sans CJK SC", "Source Han Sans SC",
                                   "WenQuanYi Zen Hei", "SimHei", "Arial Unicode MS", "sans-serif"]
    plt.rcParams["axes.unicode_minus"] = False
    return plt


INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
C_A, C_B, C_CAL = "#2a78d6", "#eb6834", "#1baf7a"


def _ax(ax):
    ax.set_facecolor(SURF)
    for sp in ax.spines.values():
        sp.set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=8)
    ax.grid(True, color=GRID, lw=0.7)


def plot_sites(out, groups, sites, preds, sets):
    plt = _style()
    cal = "留一站点标定" if any(s == "留一站点标定" for s, _ in sets) else "全部数据标定"
    for sid in sorted({g["site"] for g in groups}):
        gs = [g for g in groups if g["site"] == sid]
        ncol = min(3, len(gs))
        nrow = math.ceil(len(gs) / ncol)
        fig, axes = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 3.3 * nrow + 0.8), facecolor=SURF, squeeze=False,
                                 layout="constrained")
        for ax in axes.flat[len(gs):]:
            ax.set_visible(False)
        for ax, g in zip(axes.flat, gs):
            _ax(ax)
            t = np.arange(g["n"] + 1)
            ax.plot(t, preds[(g["key"], "默认", "A")] * 100, color=C_A, lw=1.6, label="A 实测驱动·默认参数")
            ax.plot(t, preds[(g["key"], "默认", "B")] * 100, color=C_B, lw=1.6, label="B 气象驱动·默认参数")
            ax.plot(t, preds[(g["key"], cal, "A")] * 100, color=C_CAL, lw=1.6, ls=(0, (4, 2)),
                    label=f"A · {cal}")
            ax.errorbar(g["days"], g["mean"] * 100, yerr=np.nan_to_num(g["sd"] * 100), fmt="o", ms=4.5,
                        color=INK, mec=SURF, mew=0.8, elinewidth=0.8, capsize=2, zorder=4, label="实测（均值 ± 标准差）")
            ax.axhline(80, color=INK2, lw=0.7, ls=(0, (2, 3)))
            ax.set_ylim(0, 105)
            ax.set_title(f"{g['product']} · {PL_SHORT[g['placement']]}", loc="left", fontsize=10.5, color=INK)
            ax.set_xlabel("施肥后天数", color=INK2, fontsize=8.5)
            ax.set_ylabel("累计释放（%）", color=INK2, fontsize=8.5)
        h, lab = axes.flat[0].get_legend_handles_labels()
        fig.legend(h, lab, loc="outside upper right", ncol=2, frameon=False, fontsize=8.5)
        fig.suptitle(f"{sid} {sites[sid]['name']}", x=0.01, ha="left", fontsize=12.5, color=INK, fontweight="bold")
        fig.savefig(out / f"fig_{sid}.png", dpi=150)
        plt.close(fig)


def plot_scatter(out, groups, preds, snames):
    plt = _style()
    show = [s for s in ("默认", "留一站点标定") if s in snames]
    fig, axes = plt.subplots(1, len(show), figsize=(4.6 * len(show), 4.6), facecolor=SURF, squeeze=False)
    for ax, sname in zip(axes.flat, show):
        _ax(ax)
        for run, c in (("A", C_A), ("B", C_B)):
            xs, ys = [], []
            for g in groups:
                m = g["days"] > 0
                xs += list(g["mean"][m] * 100)
                ys += list(preds[(g["key"], sname, run)][g["days"][m]] * 100)
            ax.plot(xs, ys, "o", ms=3.5, color=c, alpha=0.75, mec="none", label=f"运行 {run}")
        ax.plot([0, 100], [0, 100], color=INK2, lw=0.8)
        for d in (-8, 8):
            ax.plot([0, 100], [d, 100 + d], color=INK2, lw=0.6, ls=(0, (2, 3)))
        ax.set_xlim(0, 100)
        ax.set_ylim(0, 100)
        ax.set_aspect("equal")
        ax.set_title(sname, loc="left", fontsize=10.5, color=INK)
        ax.set_xlabel("实测累计释放（%）", color=INK2, fontsize=8.5)
        ax.set_ylabel("预测累计释放（%）", color=INK2, fontsize=8.5)
        ax.legend(frameon=False, fontsize=8.5, loc="upper left")
    fig.tight_layout()
    fig.savefig(out / "fig_scatter.png", dpi=150)
    plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser(description="田间埋袋释放试验分析")
    ap.add_argument("folder", help="数据目录，内含 sites/products/bags/weather 等 CSV")
    ap.add_argument("-o", "--out", help="结果目录（默认：数据目录下的 output/）")
    a = ap.parse_args(argv)
    folder = Path(a.folder)
    out = Path(a.out) if a.out else folder / "output"
    try:
        return analyse(folder, out)
    except DataError as e:
        print(f"数据有误：{e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
