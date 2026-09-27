"""包膜尿素温度—释放模型。

“控释期”D = 25 °C 静水中累计释放 80% 的天数（与 GB/T 23348、HG/T 4215 定义一致）。
温度只改变释放速率、不改变曲线形状（Hara 2000）。用 Arrhenius 方程把逐日温度折算成
“25 °C 等效天数” τ：
    a(T) = exp[ Ea/R · (1/298.15 − 1/(T+273.15)) ]，   τ(t) = f_soil · Σ a(T_d)
曲线形状：带初期溶出的 Weibull 函数
    F(τ) = f0 + (1−f0)·[1 − exp(−(τ/λ)^β)]，由 F(D)=0.80 解出 λ。
"""
import math

import numpy as np

from .params import P

R_GAS = 8.314


def arrhenius(T, Ea):
    """温度 T（°C）下相对 25 °C 的释放速率倍数。"""
    return np.exp(Ea * 1000 / R_GAS * (1 / 298.15 - 1 / (np.asarray(T) + 273.15)))


def lam(D, f0, beta):
    """Weibull 尺度参数，使 F(D) = 0.80。"""
    return D / (math.log((1 - f0) / 0.2)) ** (1 / beta)


def release(tau, D, f0=None, beta=None):
    """控释期 D 的产品在 25 °C 等效天数 tau 时的累计释放比例。"""
    f0 = P["f0"] if f0 is None else f0
    beta = P["beta"] if beta is None else beta
    tau = np.maximum(tau, 0)
    return f0 + (1 - f0) * (1 - np.exp(-(tau / lam(D, f0, beta)) ** beta))


def cum_tau(Tpaddy, Ea, f_soil):
    """第 t 天结束时累计的 25 °C 等效天数。"""
    return np.cumsum(arrhenius(Tpaddy, Ea) * f_soil)


def urea_increment(n, prm=P):
    """1 kg 尿素 N 每天进入有效氮库的量（扣除初始损失）。"""
    t = np.arange(n + 1)
    cum = 1 - np.exp(-t / prm["urea_tau"])
    return np.diff(cum) * (1 - prm["urea_loss"])


def cr_increment(tau, D):
    """1 kg 包膜尿素 N 每天释放的量。"""
    cum = np.concatenate([[0.0], release(tau, D)])
    return np.diff(cum)
