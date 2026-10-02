"""田间埋袋验证用的逐日释放模型：用逐日（或逐时）的实测或估算环境数据驱动 common/field_release.py。

与 field_release.simulate 的区别：
  - simulate 用 12 个月均温和种植系统预设，适合地图和模拟器；
  - 这里直接接受逐日数据：实测土温（逐时）、实测或水量平衡估算的基质势、管理事件（淹水、排水、覆膜）。
释放公式、温度因子、水分因子、施肥位置与质地系数都直接调用 field_release，保证与地图模型一致。

第 d 天的累计释放 F[d] 指施肥后经过 d 整天的值（F[0] = 初期溶出），与取袋日“施肥后第 d 天”对应。
"""
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "common"))
import field_release as FR  # noqa: E402

# van Genuchten 默认参数（Carsel & Parrish 1988；黏土取黏壤土、砂土取砂壤土的值），α 单位 1/kPa
VG_DEFAULT = {
    "clay": (0.095, 0.41, 0.019 * 10.2, 1.31),
    "loam": (0.078, 0.43, 0.036 * 10.2, 1.56),
    "sand": (0.065, 0.41, 0.075 * 10.2, 1.89),
}

# 旱地水量平衡用的作物系数（季节平均，粗略值）
SYSTEM_KC = {
    "rice_flooded": 1.05, "rice_awd": 1.05, "rice_dsr": 1.0,
    "wheat_irrigated": 0.75, "wheat_rainfed": 0.75,
    "maize_open": 0.8, "maize_film": 0.8, "notill_residue": 0.7,
}
BUCKET_MM = 200.0     # 水量平衡土层厚度（mm）
FILM_ET = 0.6         # 覆膜期间蒸散按 60% 计
SNOW_MELT = 3.0       # 积雪融化 mm/(°C·天)

LOG33, LOG1500 = math.log10(33), math.log10(1500)


# ---------------------------------------------------------------- 土壤水分
def vg_theta(psi, tr, ts, a, n):
    h = np.maximum(-np.asarray(psi, float), 0.0)
    m = 1 - 1 / n
    return tr + (ts - tr) / (1 + (a * h) ** n) ** m


def vg_psi(theta, tr, ts, a, n):
    m = 1 - 1 / n
    se = np.clip((np.asarray(theta, float) - tr) / (ts - tr), 1e-4, 1.0)
    return -((se ** (-1 / m) - 1) ** (1 / n)) / a


def dryness(psi):
    """基质势 → 干旱指数 x（0：−33 kPa 及更湿；1：−1500 kPa 及更干）。水分因子 = 1 − s·x，与 FR.moisture_factor 相同。"""
    lg = np.log10(np.maximum(-np.asarray(psi, float), 1e-9))
    return np.clip((lg - LOG33) / (LOG1500 - LOG33), 0.0, 1.0)


def ra_mm(lat, doy):
    """FAO-56 地外辐射，折算为蒸发当量 mm/天。"""
    phi = math.radians(lat)
    doy = np.asarray(doy, float)
    dr = 1 + 0.033 * np.cos(2 * np.pi * doy / 365)
    dec = 0.409 * np.sin(2 * np.pi * doy / 365 - 1.39)
    ws = np.arccos(np.clip(-math.tan(phi) * np.tan(dec), -1, 1))
    ra = 24 * 60 / np.pi * 0.0820 * dr * (ws * math.sin(phi) * np.sin(dec) + math.cos(phi) * np.cos(dec) * np.sin(ws))
    return 0.408 * np.maximum(ra, 0)


def et0_hargreaves(tmax, tmin, lat, doy):
    tmean = (tmax + tmin) / 2
    return np.maximum(0.0023 * ra_mm(lat, doy) * (tmean + 17.8) * np.sqrt(np.maximum(tmax - tmin, 0)), 0)


def bucket(precip, irrig, et0, tmean, kc, vg, flooded, film, depth=BUCKET_MM, et_scale=1.0):
    """单层水量平衡 → 逐日基质势（kPa）。flooded 为 True 的日子 ψ=0、土层饱和。"""
    tr, ts, a, n = vg
    fc, wp = float(vg_theta(-33, *vg)), float(vg_theta(-1500, *vg))
    W, snow = fc * depth, 0.0
    psi = np.zeros(len(precip))
    for d in range(len(precip)):
        if flooded[d]:
            W = ts * depth
            psi[d] = 0.0
            continue
        p = precip[d]
        if tmean[d] < 0:
            snow, p = snow + p, 0.0
        elif snow > 0:
            melt = min(snow, SNOW_MELT * tmean[d])
            snow, p = snow - melt, p + melt
        W = min(W + p + irrig[d], ts * depth)
        if W > fc * depth:
            W -= 0.5 * (W - fc * depth)
        ks = np.clip((W - 0.7 * wp * depth) / (0.5 * (fc - wp) * depth), 0, 1)
        W -= kc * et_scale * et0[d] * (FILM_ET if film[d] else 1.0) * ks
        W = max(W, (tr + 0.005) * depth)
        psi[d] = float(vg_psi(W / depth, *vg))
    return psi


# ---------------------------------------------------------------- 温度
def day_factor_hourly(temps, Ea):
    """一天内逐时土温 → 当天平均温度因子（Arrhenius × 冻结因子）。"""
    T = np.asarray(temps, float)
    return float((FR.arrhenius(T, Ea) * FR.freeze_ramp(T)).mean())


def air_to_factor(Tair, system, placement, Ea):
    """气温 → 包膜处土温（种植系统预设 + 施肥位置）→ 温度因子；与 FR.simulate 的做法相同。"""
    Ts = FR.soil_temperature(np.asarray(Tair, float), system, placement)
    amp = FR.PLACEMENT[placement]["amp"]
    amp = FR.SYSTEMS[system]["amp"] if amp is None else amp
    return FR.temp_factor(Ts, amp, Ea), Ts


# ---------------------------------------------------------------- 释放
def release_curve(rate, D, beta, f0, damage):
    """逐日速率（25 °C 等效天/天）→ F[0..n]。"""
    tau = np.concatenate([[0.0], np.cumsum(rate)])
    f0e = f0 + (1 - f0) * damage
    return FR.release(tau, D, f0e, beta)


def day_at(F, p, days=None):
    """F 首次达到 p 的天数（线性插值）；没达到返回 None。"""
    F = np.asarray(F, float)
    days = np.arange(len(F)) if days is None else np.asarray(days, float)
    idx = np.where(F >= p)[0]
    if not len(idx):
        return None
    i = idx[0]
    if i == 0:
        return float(days[0])
    return float(days[i - 1] + (p - F[i - 1]) / (F[i] - F[i - 1]) * (days[i] - days[i - 1]))
