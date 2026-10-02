"""生成一套【模拟】示例数据，用来演示 analyze.py 和检查它能否找回已知参数。不是实测数据。

做法：设定一组“真实”参数（与模型默认值有差别），用合成气象、合成的逐时土温和土壤水分，
按田间释放模型算出“真实”释放过程，再加上取袋误差（每袋标准差 2.5 个百分点）。
运行 A 用的土温、水分就是这些合成“实测”值；运行 B 只用合成气象。
真实土温与模型预设的土温偏移、昼夜温差不同，所以运行 B 会有额外误差。

    python validation/make_example.py        # 写入 validation/example/
"""
import csv
import datetime as dt
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import vmodel as VM  # noqa: E402

FR = VM.FR
OUT = Path(__file__).resolve().parent / "example"
RNG = np.random.default_rng(20261002)

TRUE = dict(f_soil=0.80, moisture_sens=0.35, band=0.80, surface=2.0, clay=1.10, sand=0.95)

PRODUCTS = {   # 假设的实验室定标结果
    "CR40": dict(lot="A2501", D=40, beta=1.25, f0=0.035, Ea=48.0, damage=0.02),
    "CR90": dict(lot="A2502", D=90, beta=1.30, f0=0.030, Ea=52.0, damage=0.02),
    "CR150": dict(lot="A2503", D=150, beta=1.35, f0=0.025, Ea=55.0, damage=0.015),
}

# 每个站点：月均温、月降水与雨日、真实土温偏移与昼夜温差半幅、管理事件（天数, 事件, 水量）
SITES = [
    dict(id="HB", name="（模拟）湖北·中稻", lat=30.6, lon=114.3, system="rice_flooded", texture="clay",
         apply=dt.date(2025, 6, 5), days=120, products=["CR40", "CR90"], placements=["incorporated", "band"],
         T12=[4, 6.5, 11, 17, 22, 26, 29, 28.5, 24, 18, 11.5, 6],
         P12=[45, 60, 95, 135, 160, 220, 190, 120, 75, 75, 55, 30], R12=[8, 9, 12, 13, 13, 14, 12, 10, 8, 9, 8, 6],
         offset=1.6, amp=2.0, moisture_sensor=None,
         events=[(0, "flood", 0), (35, "drain", 0), (45, "flood", 0), (100, "drain", 0)]),
    dict(id="HN", name="（模拟）河南·灌溉冬小麦", lat=34.7, lon=113.6, system="wheat_irrigated", texture="loam",
         apply=dt.date(2024, 10, 12), days=240, products=["CR90", "CR150"], placements=["incorporated", "surface"],
         T12=[0.5, 3.5, 9, 16, 21.5, 26.5, 27.5, 26, 21.5, 16, 8.5, 2.5],
         P12=[10, 13, 30, 40, 50, 70, 155, 125, 70, 40, 25, 10], R12=[3, 4, 6, 6, 7, 8, 12, 11, 8, 6, 4, 3],
         offset=1.5, amp=5.0, moisture_sensor="gravimetric",
         events=[(0, "irrigate", 60), (50, "irrigate", 75), (150, "irrigate", 75), (190, "irrigate", 70)]),
    dict(id="SX", name="（模拟）陕西·旱作冬小麦", lat=35.2, lon=107.7, system="wheat_rainfed", texture="sand",
         apply=dt.date(2024, 9, 28), days=255, products=["CR90", "CR150"], placements=["incorporated"],
         T12=[-4, -1, 5, 11.5, 16.5, 20.5, 22.5, 21, 16, 10, 3.5, -2.5],
         P12=[6, 9, 25, 38, 55, 70, 115, 115, 95, 50, 18, 5], R12=[3, 4, 7, 8, 9, 10, 12, 12, 12, 9, 5, 3],
         offset=1.2, amp=5.5, moisture_sensor="fdr", events=[]),
    dict(id="GS", name="（模拟）甘肃·地膜玉米", lat=37.9, lon=102.6, system="maize_film", texture="loam",
         apply=dt.date(2025, 4, 20), days=160, products=["CR90", "CR150"], placements=["band"],
         T12=[-7.5, -3, 4, 11.5, 17, 21, 23, 21.5, 16.5, 9, 1.5, -5.5],
         P12=[1, 2, 5, 8, 15, 20, 35, 35, 25, 8, 2, 1], R12=[1, 1, 2, 3, 4, 5, 7, 7, 5, 3, 1, 1],
         offset=None, amp=6.5, moisture_sensor="fdr",
         events=[(0, "film_on", 0), (0, "irrigate", 40), (45, "irrigate", 90), (75, "irrigate", 90),
                 (100, "irrigate", 90), (125, "irrigate", 80)]),
]

LINEAR_DAYS = [0, 3, 7, 14, 21, 28, 42, 56, 70, 90, 110]
WHEAT_DAYS = [0, 3, 7, 14, 21, 28, 42, 56, 75, 105, 135, 150, 165, 180, 195, 210, 225, 240, 255]


def monthly(arr, dates):
    m = np.array([d.month - 1 for d in dates])
    return np.asarray(arr, float)[m]


def synth_weather(s, n):
    dates = [s["apply"] + dt.timedelta(days=k) for k in range(n)]
    doy = np.array([d.timetuple().tm_yday for d in dates])
    base = FR.daily_temps(np.array(s["T12"], float), int(doy[0]) - 1, n)
    noise = np.zeros(n)
    for k in range(1, n):
        noise[k] = 0.7 * noise[k - 1] + RNG.normal(0, 1.8)
    tmean = base + noise
    dtr = np.clip(RNG.normal(10 if s["system"] != "rice_flooded" else 7.5, 2, n), 3, 18)
    p_wet = monthly(s["R12"], dates) / 30.4
    mean_amt = monthly(s["P12"], dates) / np.maximum(monthly(s["R12"], dates), 0.5)
    wet = RNG.random(n) < p_wet
    precip = np.where(wet, RNG.exponential(1.0, n) * mean_amt, 0.0)
    dtr = np.where(wet, dtr * 0.6, dtr)
    return dates, doy, tmean, tmean + dtr / 2, tmean - dtr / 2, np.round(precip, 1)


def true_offset(s, k):
    if s["system"] == "maize_film":                       # 地膜增温：前 50 天约 +4 °C，之后冠层遮阴逐渐降到 +1
        return 4.0 if k < 50 else max(1.0, 4.0 - (k - 50) * 0.06)
    return s["offset"]


def true_soil_hourly(s, tmean, placement):
    """“真实”逐时土温：气温平滑（土壤热惯性）+ 偏移 + 昼夜波动；冻结期按积雪、土层缓冲。"""
    n = len(tmean)
    sm = np.copy(tmean)
    for k in range(1, n):
        sm[k] = 0.55 * tmean[k] + 0.45 * sm[k - 1]
    surface = placement == "surface"
    hours = np.arange(24)
    out = np.zeros((n, 24))
    for k in range(n):
        if sm[k] >= 0:
            Tm = sm[k] + true_offset(s, k) + (2.5 if surface else 0)
            amp = 9.0 if surface else s["amp"]
        else:
            Tm = sm[k] * (0.6 if surface else 0.4) + 0.5
            amp = 3.0 if surface else 1.0
        out[k] = Tm + amp * np.sin(2 * np.pi * (hours - 9) / 24) + RNG.normal(0, 0.15, 24)
    return np.round(out, 2)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rows = {k: [] for k in ("sites", "products", "bags", "soil_temp", "soil_moisture", "weather", "management")}
    for name, p in PRODUCTS.items():
        rows["products"].append([name, p["lot"], p["D"], p["beta"], round(p["f0"] * 100, 2), p["Ea"],
                                 round(p["damage"] * 100, 2)])
    for s in SITES:
        n = s["days"] + 1
        vg = VM.VG_DEFAULT[s["texture"]]
        rows["sites"].append([s["id"], s["name"], s["lat"], s["lon"], s["system"],
                              {"clay": "黏土", "loam": "壤土", "sand": "砂土"}[s["texture"]],
                              s["apply"].isoformat(), 5 if s["system"].startswith("rice") else 8, "", "", "", "",
                              "模拟数据，仅作演示"])
        dates, doy, tmean, tmax, tmin, precip = synth_weather(s, n)
        for k in range(n):
            rows["weather"].append([s["id"], dates[k].isoformat(), round(tmax[k], 1), round(tmin[k], 1),
                                    round(tmean[k], 1), precip[k]])
        ev = {(d, e): a for d, e, a in s["events"]}
        for d, e, a in s["events"]:
            label = {"flood": "淹水", "drain": "排水", "irrigate": "灌溉", "film_on": "覆膜"}[e]
            rows["management"].append([s["id"], dates[d].isoformat(), label, a or "", ""])
        # 真实水分：比分析用的水量平衡更浅的土层、更大的蒸散
        flooded = np.zeros(n, bool)
        state = False
        irrig = np.zeros(n)
        for k in range(n):
            if (k, "flood") in ev:
                state = True
            if (k, "drain") in ev:
                state = False
            flooded[k] = state
            irrig[k] = ev.get((k, "irrigate"), 0.0)
        film = np.full(n, s["system"] == "maize_film")
        et0 = VM.et0_hargreaves(tmax, tmin, s["lat"], doy)
        psi_true = VM.bucket(precip, irrig, et0, tmean, VM.SYSTEM_KC[s["system"]], vg, flooded, film,
                             depth=140.0, et_scale=1.15)
        theta_true = VM.vg_theta(psi_true, *vg)
        if s["moisture_sensor"] == "fdr":
            for k in range(n):
                for h in (8, 20):
                    rows["soil_moisture"].append([s["id"], f"{dates[k].isoformat()} {h:02d}:00",
                                                  round(float(theta_true[k]) + RNG.normal(0, 0.008), 3), ""])
        days_list = WHEAT_DAYS if s["system"].startswith("wheat") else LINEAR_DAYS
        days_list = [d for d in days_list if d <= s["days"]]
        if s["moisture_sensor"] == "gravimetric":       # 只在取袋日用称重法测一次
            for d in days_list:
                rows["soil_moisture"].append([s["id"], f"{dates[d].isoformat()} 10:00",
                                              round(float(theta_true[d]) + RNG.normal(0, 0.015), 3), ""])
        x_true = VM.dryness(psi_true)
        for pl in s["placements"]:
            hourly = true_soil_hourly(s, tmean, pl)
            for k in range(n):
                for h in range(24):
                    rows["soil_temp"].append([s["id"], {"incorporated": "混施", "band": "条施", "surface": "表施"}[pl],
                                              f"{dates[k].isoformat()} {h:02d}:00", hourly[k, h]])
            for prod in s["products"]:
                p = PRODUCTS[prod]
                aT = np.array([VM.day_factor_hourly(hourly[k], p["Ea"]) for k in range(n)])
                pl_f = 1.0 if pl == "incorporated" else TRUE[pl]
                tex_f = 1.0 if s["texture"] == "loam" else TRUE[s["texture"]]
                rate = aT * (1 - TRUE["moisture_sens"] * x_true) * pl_f * tex_f * TRUE["f_soil"]
                F = VM.release_curve(rate, p["D"], p["beta"], p["f0"], p["damage"])
                for d in days_list:
                    for rep in range(1, 4):
                        f = float(np.clip(F[d] + RNG.normal(0, 0.025), 0, 1))
                        n_init = round(RNG.uniform(3.4, 3.6), 4)          # 约 7.5 g 包膜尿素
                        rows["bags"].append([s["id"], prod, {"incorporated": "混施", "band": "条施",
                                                             "surface": "表施"}[pl],
                                             f"{s['id']}-{prod}-{pl[0].upper()}-{d:03d}-{rep}", dates[d].isoformat(),
                                             n_init, round(n_init * (1 - f), 4), "", "", "", "", "模拟"])
                    if F[d] > 0.92:
                        break

    header = {
        "sites": ["站点", "名称", "纬度", "经度", "种植系统", "质地", "施肥日期", "埋袋深度_cm",
                  "theta_r", "theta_s", "alpha_per_kPa", "n_vg", "备注"],
        "products": ["产品", "批号", "D_天", "beta", "f0", "Ea_kJ", "破损率"],
        "bags": ["站点", "产品", "施肥位置", "袋号", "取样日期", "初始氮_g", "剩余氮_g", "累计释放%",
                 "完整粒", "空壳", "破损粒", "备注"],
        "soil_temp": ["站点", "施肥位置", "时间", "土温_C"],
        "soil_moisture": ["站点", "时间", "含水量_m3m3", "基质势_kPa"],
        "weather": ["站点", "日期", "最高气温_C", "最低气温_C", "平均气温_C", "降水_mm"],
        "management": ["站点", "日期", "事件", "水量_mm", "备注"],
    }
    for k, rs in rows.items():
        with open(OUT / f"{k}.csv", "w", encoding="utf-8-sig", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(header[k])
            w.writerows(rs)
    (OUT / "说明.txt").write_text(
        "本目录全部是 make_example.py 生成的【模拟数据】，不是实测数据，只用于演示分析流程。\n"
        "生成时设定的“真实”参数：" + "，".join(f"{k}={v}" for k, v in TRUE.items()) + "\n"
        "分析结果中的“全部数据标定”参数越接近这些值，说明标定方法越可靠。\n", encoding="utf-8")
    print(f"已写入 {OUT}（模拟数据）")


if __name__ == "__main__":
    main()
