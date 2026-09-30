"""线性型与 S 型包膜尿素在 25 °C 静水中的释放示意图（控释期 D=90 天）。

F(t) = f0 + (1−f0)·[1 − exp(−(t/λ)^β)]，λ 由 F(D)=80% 求出，f0=3%（初期溶出）。
运行：python analysis/release_shapes.py  →  analysis/release_shapes.png
"""
import math
from pathlib import Path

import logging
import matplotlib
matplotlib.use("Agg")
logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams["font.family"] = ["WenQuanYi Zen Hei", "sans-serif"]
plt.rcParams["axes.unicode_minus"] = False

D, F0 = 90, 0.03
CURVES = [  # (β, 名称, 颜色, 线型)
    (1.3, "线性型  β=1.3（现有产品）", "#2a78d6", "-"),
    (2.5, "S 型  β=2.5", "#eb6834", "--"),
    (3.5, "S 型  β=3.5", "#1baf7a", ":"),
]
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"


def lam(b):
    return D / (math.log((1 - F0) / 0.2)) ** (1 / b)


def F(t, b):
    return F0 + (1 - F0) * (1 - np.exp(-(t / lam(b)) ** b))


def t_at(p, b):
    return lam(b) * (-math.log(1 - (p - F0) / (1 - F0))) ** (1 / b)


t = np.linspace(0, 180, 1801)
fig, (a1, a2) = plt.subplots(1, 2, figsize=(12.5, 5.2), dpi=200, facecolor=SURF)
for ax in (a1, a2):
    ax.set_facecolor(SURF)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=9)
    ax.grid(True, color=GRID, linewidth=0.8)
    ax.set_xlim(0, 180)
    ax.set_xticks(range(0, 181, 30))
    ax.set_xlabel("浸泡天数（25 °C 静水）", color=INK2, fontsize=10)

# 左：累计释放
for b, name, c, ls in CURVES:
    a1.plot(t, 100 * F(t, b), color=c, lw=2.2, ls=ls, solid_capstyle="round")
a1.set_ylim(0, 102)
a1.set_ylabel("累计释放（%）", color=INK2, fontsize=10)
a1.axhline(80, color=INK2, lw=0.9, ls=(0, (2, 3)))
a1.axvline(D, color=INK2, lw=0.9, ls=(0, (2, 3)))
a1.plot([D], [80], "o", ms=7, color=INK, mec=SURF, mew=1.5, zorder=5)
a1.annotate("三条曲线都在第 90 天释放 80%\n（同为“控释期 90 天”）", (D, 80), (100, 58),
            fontsize=9, color=INK, arrowprops=dict(arrowstyle="-", color=INK2, lw=0.8))
# 直接标注
for (b, name, c, ls), (x, dy) in zip(CURVES, [(30, 5), (52, -4), (40, -3)]):
    y = 100 * F(np.array([x]), b)[0]
    a1.annotate(name, (x, y), (x + (-28 if b == 1.3 else 6), y + (8 if b == 1.3 else -12 if b == 2.5 else -9)),
                fontsize=9, color=INK, arrowprops=dict(arrowstyle="-", color=c, lw=1))
# 滞后期（β=3.5 释放不到 10% 的阶段）
lag = t_at(0.10, 3.5)
a1.axvspan(0, lag, color="#1baf7a", alpha=0.08, lw=0)
a1.text(lag / 2, 96, f"滞后期\n（β=3.5 前 {lag:.0f} 天\n释放不到 10%）", ha="center", va="top", fontsize=8.5, color=INK2)
a1.set_title("累计释放曲线", loc="left", fontsize=12, color=INK, pad=10)

# 右：日释放速率
for b, name, c, ls in CURVES:
    r = np.gradient(100 * F(t, b), t)
    a2.plot(t, r, color=c, lw=2.2, ls=ls)
    pk = lam(b) * ((b - 1) / b) ** (1 / b)
    rp = np.interp(pk, t, r)
    a2.plot([pk], [rp], "o", ms=6, color=c, mec=SURF, mew=1.5, zorder=5)
    lx, ly = {1.3: (3, 1.36), 2.5: (98, 1.48), 3.5: (84, 1.70)}[b]
    a2.annotate(f"{name.split('（')[0]}\n高峰第 {pk:.0f} 天", (pk, rp), (lx, ly),
                fontsize=9, color=INK, arrowprops=dict(arrowstyle="-", color=c, lw=1))
a2.set_ylim(0, None)
a2.set_ylabel("日释放速率（% / 天）", color=INK2, fontsize=10)
a2.set_title("日释放速率", loc="left", fontsize=12, color=INK, pad=10)

fig.suptitle("线性型与 S 型包膜尿素在 25 °C 静水中的释放（示意，控释期 D = 90 天）",
             x=0.01, ha="left", fontsize=13.5, color=INK, fontweight="bold")
fig.text(0.01, 0.005, "F(t) = f₀ + (1 − f₀)·[1 − exp(−(t/λ)^β)]，f₀ = 3%，λ 由 F(90 天) = 80% 求出。"
         "β 越大滞后越明显、释放越集中；其他控释期的曲线按 D 等比例横向伸缩。",
         fontsize=8.5, color=INK2)
fig.tight_layout(rect=(0, 0.03, 1, 0.94))
out = Path(__file__).with_suffix(".png")
fig.savefig(out, facecolor=SURF)
print(out)
