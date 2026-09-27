"""气温：12 个月均温周期插值成逐日气温；日期换算。"""
import numpy as np


def daily_temps(T12, start_doy, n):
    """由 12 个月均温（各月中点）周期线性插值，得到从 start_doy 起 n 天的逐日气温。"""
    mid = np.array([15.2 + 30.44 * m for m in range(12)])
    xs = np.concatenate([mid - 365, mid, mid + 365])
    ys = np.concatenate([T12, T12, T12])
    d = (start_doy + np.arange(n)) % 365
    return np.interp(d, xs, ys)


def doy(mmdd):
    """'MM-DD' → 一年中的第几天（从 0 起，不计闰年）。"""
    m, d = map(int, mmdd.split("-"))
    return [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334][m - 1] + d - 1
