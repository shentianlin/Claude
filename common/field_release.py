"""包膜尿素田间释放模型：把 25 °C 静水释放曲线换算成不同种植系统下的田间释放。

    F(τ) = f0' + (1 − f0')·[1 − exp(−(τ/λ)^β)]，λ 由 F(D) = 80% 求出（与地图模型相同）
    dτ/dt = a_T · a_W · a_P · a_S · f_soil                （每天累计的 25 °C 等效天数）

各因子与文献依据
  a_T 温度：Arrhenius，按包膜所在位置的温度计算，并对昼夜波动积分（Arrhenius 是凸函数，
      同样的日均温，昼夜温差越大释放越快）；土温 0–3 °C 线性降到 0（冻结后水分无法扩散进颗粒）。
      温度是田间释放的首要因素（Golden 2011；Pang 2026；El-Hout 2021 用积温预测 PCU 释放）。
      恒温与波动温度下实测差别很大（Ransom 2020）。
  a_W 土壤水分：按基质势 ψ。−33 kPa 以上（含淹水）取 1；更干时按 log10(−ψ) 线性下降，
      到 −1500 kPa 时为 1 − s。s 为产品的“水分敏感度”（与膜材有关），默认 0.25：20 °C 下 CR60
      在 −1000 kPa 时比 −10 kPa 推迟约 27 天，与 Verburg 2020 的“最多推迟约 30 天”相当。作物生长适宜湿度范围内几乎无影响
      （Sentek 2023；Ransom 2020）。
  a_P 施肥位置：混施入土 = 1；条施、侧深施 = 0.85（颗粒集中，扩散受限，Janke 2022）；
      表面撒施不覆土 = 2.5（Ransom 2020：表施的 45–180 天产品都在 35 天内释放完），且表面温度、
      昼夜温差按地表取值。表面撒施的系数不确定性很大。
  a_S 土壤质地：黏土 1.15、壤土 1.0、砂土 0.95（Golden 2011：黏质土释放更快）。证据有限，默认壤土。
  f_soil 土壤中相对静水的基础系数 0.9（与水稻地图一致；小麦地图的 0.8 ≈ 0.9 × 旱地水分因子）。
  破损率 d：搬运、掺混造成的破损颗粒在施后立即全部溶出，计入初期溶出：f0' = f0 + (1 − f0)·d。

种植系统预设（SYSTEMS）规定：包膜处土温相对气温的偏移（可随时间变化）、昼夜温差半幅、
水分状态（基质势随时间的变化）。同一作物的不同栽培方式会落在不同预设上。
"""
import math

import numpy as np

R_GAS = 8.314

# 水分状态：基质势在 wet 与 dry（kPa，负值）之间按周期 period 天余弦摆动；flood=True 表示淹水（ψ=0）
MOISTURE = {
    "flooded": dict(label="淹水", flood=True),
    "awd": dict(label="干湿交替灌溉（落干到约 −20 kPa）", wet=-1.0, dry=-20.0, period=10),
    "irrigated": dict(label="灌溉旱地（−10 ~ −80 kPa）", wet=-10.0, dry=-80.0, period=10),
    "humid": dict(label="雨养·湿润（−10 ~ −300 kPa）", wet=-10.0, dry=-300.0, period=14),
    "semiarid": dict(label="雨养·半干旱（−50 ~ −1200 kPa）", wet=-50.0, dry=-1200.0, period=20),
}

# 种植系统预设
#   offset 土温 − 气温（°C）：{"base": x} 为常数；{"early": a, "until": d1, "late": b, "ramp": r}
#          表示前 d1 天为 a，之后在 r 天内线性过渡到 b（地膜覆盖、秸秆覆盖等）
#   amp 昼夜温差半幅（°C，包膜所在深度），moisture 水分状态，
#   switch 可选：[第几天, 新水分状态]，例如旱直播稻播后约 30 天淹水
SYSTEMS = {
    "rice_flooded": dict(label="水稻·淹水（移栽或水直播）", offset={"base": 1.0}, amp=2.5, moisture="flooded"),
    "rice_awd": dict(label="水稻·干湿交替灌溉", offset={"base": 1.0}, amp=3.5, moisture="awd"),
    "rice_dsr": dict(label="水稻·旱直播，约 30 天后淹水", offset={"base": 1.0}, amp=4.0, moisture="humid",
                     switch=[30, "flooded"]),
    "wheat_irrigated": dict(label="小麦·灌溉", offset={"base": 1.0}, amp=4.0, moisture="irrigated"),
    "wheat_rainfed": dict(label="小麦·旱作", offset={"base": 1.0}, amp=5.0, moisture="semiarid"),
    "maize_open": dict(label="玉米·露地", offset={"base": 1.0}, amp=4.5, moisture="humid"),
    "maize_film": dict(label="玉米·地膜覆盖", offset={"early": 3.5, "until": 60, "late": 1.0, "ramp": 40},
                       amp=6.0, moisture="irrigated"),
    "notill_residue": dict(label="免耕·秸秆覆盖", offset={"early": -1.0, "until": 45, "late": 0.0, "ramp": 0},
                           amp=3.0, moisture="humid"),
}


def offset_at(spec, d):
    if "base" in spec:
        return spec["base"]
    if d < spec["until"]:
        return spec["early"]
    if spec["ramp"] <= 0:
        return spec["late"]
    f = min(1.0, (d - spec["until"]) / spec["ramp"])
    return spec["early"] + (spec["late"] - spec["early"]) * f


PLACEMENT = {
    "incorporated": dict(label="撒施后混入土壤", factor=1.0, dT=0.0, amp=None),
    "band": dict(label="条施 / 侧深施", factor=0.85, dT=0.0, amp=None),
    "surface": dict(label="表面撒施不覆土", factor=2.5, dT=2.0, amp=10.0),
}
TEXTURE = {"clay": dict(label="黏土", factor=1.15), "loam": dict(label="壤土", factor=1.0),
           "sand": dict(label="砂土", factor=0.95)}

DEFAULTS = dict(Ea=50.0, f0=0.03, beta=1.3, f_soil=0.90, moisture_sens=0.25, damage=0.02)


# ---------------------------------------------------------------- 基本函数
def daily_temps(T12, start_doy, n):
    mid = np.array([15.2 + 30.44 * m for m in range(12)])
    xs = np.concatenate([mid - 365, mid, mid + 365])
    ys = np.concatenate([T12, T12, T12])
    return np.interp((start_doy + np.arange(n)) % 365, xs, ys)


def arrhenius(T, Ea):
    return np.exp(Ea * 1000 / R_GAS * (1 / 298.15 - 1 / (np.asarray(T, float) + 273.15)))


def freeze_ramp(T):
    return np.clip(np.asarray(T, float) / 3.0, 0.0, 1.0)


def temp_factor(Tmean, amp, Ea, n_hours=24):
    """对一天内的正弦温度波动积分 Arrhenius × 冻结因子；amp 为半幅。"""
    h = np.sin(2 * np.pi * (np.arange(n_hours) + 0.5) / n_hours)
    T = np.asarray(Tmean, float)[..., None] + amp * h
    return (arrhenius(T, Ea) * freeze_ramp(T)).mean(axis=-1)


def matric_potential(regime, n):
    m = MOISTURE[regime]
    if m.get("flood"):
        return np.zeros(n)
    phase = (1 - np.cos(2 * np.pi * np.arange(n) / m["period"])) / 2
    return -np.exp(np.log(-m["wet"]) + (np.log(-m["dry"]) - np.log(-m["wet"])) * phase)


def moisture_factor(psi, s):
    """基质势 ψ（kPa，≤0）→ 释放速率倍数。"""
    x = np.clip((np.log10(np.maximum(-np.asarray(psi, float), 1e-9)) - math.log10(33)) /
                (math.log10(1500) - math.log10(33)), 0.0, 1.0)
    return 1 - s * x


def lam(D, f0, beta):
    return D / (math.log((1 - f0) / 0.2)) ** (1 / beta)


def release(tau, D, f0, beta):
    return f0 + (1 - f0) * (1 - np.exp(-(np.maximum(tau, 0) / lam(D, f0, beta)) ** beta))


def soil_temperature(Tair, system, placement):
    """包膜所在位置的日均温。气温低于 0 °C 时按积雪、土层缓冲取 0.4·气温 + 0.5（与小麦模型一致）。"""
    sysd, pl = SYSTEMS[system], PLACEMENT[placement]
    off = np.array([offset_at(sysd["offset"], d) for d in range(len(Tair))]) + pl["dT"]
    T = np.where(Tair >= 0, Tair + off, Tair * 0.4 + 0.5)
    return T


# ---------------------------------------------------------------- 模拟
def simulate(T12, start_doy, days, D, beta=None, f0=None, Ea=None, system="rice_flooded",
             placement="incorporated", texture="loam", moisture_sens=None, damage=None, f_soil=None):
    """返回逐日结果与汇总。"""
    beta = DEFAULTS["beta"] if beta is None else beta
    f0 = DEFAULTS["f0"] if f0 is None else f0
    Ea = DEFAULTS["Ea"] if Ea is None else Ea
    s = DEFAULTS["moisture_sens"] if moisture_sens is None else moisture_sens
    dmg = DEFAULTS["damage"] if damage is None else damage
    fs = DEFAULTS["f_soil"] if f_soil is None else f_soil
    sysd, pl = SYSTEMS[system], PLACEMENT[placement]
    n = days + 1
    Tair = daily_temps(T12, start_doy, n)
    Ts = soil_temperature(Tair, system, placement)
    amp = pl["amp"] if pl["amp"] is not None else sysd["amp"]
    aT = temp_factor(Ts, amp, Ea)
    psi = matric_potential(sysd["moisture"], n)
    if "switch" in sysd:
        d0, reg = sysd["switch"]
        psi[d0:] = matric_potential(reg, n - d0)
    aW = moisture_factor(psi, s)
    aP = pl["factor"]
    aS = TEXTURE[texture]["factor"]
    rate = aT * aW * aP * aS * fs
    tau = np.cumsum(rate)
    f0e = f0 + (1 - f0) * dmg
    F = release(tau, D, f0e, beta)

    def day_at(p):
        idx = np.where(F >= p)[0]
        return int(idx[0]) if len(idx) else None

    return dict(Tair=Tair, Tsoil=Ts, aT=aT, aW=aW, rate=rate, tau=tau, F=F, psi=psi,
                f0_eff=f0e, d50=day_at(0.5), d80=day_at(0.8), F_end=float(F[-1]),
                mean_rate=float(rate.mean()))
