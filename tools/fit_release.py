"""包膜尿素静水释放曲线拟合工具。

把各档产品的“天数—累计释放%”数据填进 CSV，本工具逐档拟合
    F(t) = f0 + (1 − f0)·[1 − exp(−(t/λ)^β)]
输出控释期 D（累计释放 80% 的天数）、形状参数 β、初期溶出 f0、拟合优度，
并画出实测点与拟合曲线的对比图。同一档位有多个温度的数据时，另按 Arrhenius 方程估算表观活化能 Ea。

用法
    python tools/fit_release.py 数据.csv                 # 结果写到 数据.csv 所在目录下的 fit_output/
    python tools/fit_release.py 数据.csv -o 结果目录
    python tools/fit_release.py 数据.csv --fix-f0 3      # 固定初期溶出为 3%，只拟合 λ、β

CSV 两种格式（UTF-8 或 Excel 另存的“CSV UTF-8”均可）
  长表：每行一个数据点，列名 grade/档位、day/天数、release/累计释放（%），可选 temp/温度（°C，缺省 25）
        grade,day,release,temp
        CR90,1,2.8,25
        CR90,7,9.5,25
  宽表：第一列为天数，其余每列一个档位（列名即档位名），可选第二列为温度
        day,CR60,CR90,CR120
        1,3.1,2.8,2.5
        7,14.2,9.5,7.1
  空格子会被跳过；累计释放可以写百分数（2.8）或小数（0.028）。
"""
import argparse
import csv
import logging
import math
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

R_GAS = 8.314
ALIASES = {
    "grade": {"grade", "档位", "产品", "型号", "product"},
    "day": {"day", "days", "天数", "天", "时间", "time"},
    "release": {"release", "release_pct", "累计释放", "累计释放%", "累计释放率", "释放率", "释放%", "cumulative"},
    "temp": {"temp", "temp_c", "温度", "温度°c", "温度℃", "temperature"},
}


# ---------------------------------------------------------------- 读取数据
def _norm(h):
    return h.strip().lower().replace("（", "(").replace("）", ")").replace(" ", "")


def _key(h):
    n = _norm(h)
    for k, names in ALIASES.items():
        if n in {_norm(x) for x in names} or n.split("(")[0] in {_norm(x) for x in names}:
            return k
    return None


def _num(x):
    x = (x or "").strip().replace("%", "")
    if x == "":
        return None
    return float(x)


def read_data(path):
    """返回 {(档位, 温度): [(天数, 释放%), ...]}。"""
    with open(path, newline="", encoding="utf-8-sig") as f:
        rows = [r for r in csv.reader(f) if any(c.strip() for c in r)]
    if not rows:
        sys.exit("CSV 是空的。")
    header, body = rows[0], rows[1:]
    keys = [_key(h) for h in header]
    data = defaultdict(list)
    if "grade" in keys and "release" in keys and "day" in keys:            # 长表
        ig, idd, ir = keys.index("grade"), keys.index("day"), keys.index("release")
        it = keys.index("temp") if "temp" in keys else None
        for r in body:
            d, v = _num(r[idd]), _num(r[ir])
            if d is None or v is None:
                continue
            T = _num(r[it]) if it is not None and it < len(r) else None
            data[(r[ig].strip(), 25.0 if T is None else T)].append((d, v))
    elif keys and keys[0] == "day":                                        # 宽表
        it = 1 if len(keys) > 1 and keys[1] == "temp" else None
        cols = [i for i in range(1, len(header)) if i != it]
        for r in body:
            d = _num(r[0])
            if d is None:
                continue
            T = _num(r[it]) if it is not None else None
            for i in cols:
                v = _num(r[i]) if i < len(r) else None
                if v is not None:
                    data[(header[i].strip(), 25.0 if T is None else T)].append((d, v))
    else:
        sys.exit("看不懂表头：长表需要 grade/档位、day/天数、release/累计释放 三列；宽表第一列须为 day/天数。")
    out = {}
    for k, pts in data.items():
        pts.sort()
        t = np.array([p[0] for p in pts], float)
        y = np.array([p[1] for p in pts], float)
        if y.max() <= 1.5:                     # 写成了小数
            y = y * 100
        out[k] = (t, y)
    return out


# ---------------------------------------------------------------- 模型与拟合
def model(t, f0, lam, beta):
    return 100 * (f0 + (1 - f0) * (1 - np.exp(-(np.maximum(t, 0) / lam) ** beta)))


def t_at(p, f0, lam, beta):
    """累计释放达到 p（0–1）的天数。"""
    y = (p - f0) / (1 - f0)
    if y <= 0:
        return 0.0
    return lam * (-math.log(1 - y)) ** (1 / beta)


def fit_one(t, y, fix_f0=None):
    ok = t > 0
    t, y = t[ok], y[ok]
    if len(t) < 3 + (fix_f0 is None):
        raise ValueError(f"数据点太少（{len(t)} 个），至少需要 {3 + (fix_f0 is None)} 个")
    # 初值：λ 取累计释放到 ~63% 的天数附近
    target = min(60.0, 0.8 * y.max())
    lam0 = float(np.interp(target, y, t)) if y.max() > target else float(t.max())
    f00 = min(max(y[0] / 100 * 0.8, 0.0), 0.2) if fix_f0 is None else fix_f0 / 100
    best = None
    for b0 in (0.8, 1.3, 2.0, 3.0):
        if fix_f0 is None:
            fun = lambda p: model(t, p[0], p[1], p[2]) - y
            r = least_squares(fun, [f00, lam0, b0], bounds=([0, 1e-3, 0.3], [0.5, 1e5, 10]))
            f0, lam, beta = r.x
        else:
            f0c = fix_f0 / 100
            fun = lambda p: model(t, f0c, p[0], p[1]) - y
            r = least_squares(fun, [lam0, b0], bounds=([1e-3, 0.3], [1e5, 10]))
            (lam, beta), f0 = r.x, f0c
        if best is None or r.cost < best[0]:
            best = (r.cost, f0, lam, beta)
    _, f0, lam, beta = best
    pred = model(t, f0, lam, beta)
    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    D = t_at(0.80, f0, lam, beta)
    D_obs = float(np.interp(80, y, t)) if y.max() >= 80 and np.all(np.diff(y) >= -1) else None
    return dict(
        f0=f0, lam=lam, beta=beta, D=D, D_obs=D_obs,
        t10=t_at(0.10, f0, lam, beta), t50=t_at(0.50, f0, lam, beta),
        peak=lam * ((beta - 1) / beta) ** (1 / beta) if beta > 1 else 0.0,
        r2=1 - ss_res / ss_tot if ss_tot > 0 else float("nan"),
        rmse=math.sqrt(ss_res / len(t)), n=len(t), t_max=float(t.max()), y_max=float(y.max()),
    )


def shape_label(beta):
    if beta < 1.0:
        return "前快后慢"
    if beta < 1.8:
        return "线性型"
    return "S 型"


def notes(r):
    msg = []
    if r["y_max"] < 80:
        msg.append(f"最高只测到 {r['y_max']:.0f}%，D 为外推值")
    if r["r2"] < 0.98:
        msg.append("拟合一般，曲线可能不是单一 Weibull 形状")
    if r["rmse"] > 3:
        msg.append(f"残差偏大（RMSE {r['rmse']:.1f} 个百分点）")
    if r["D_obs"] and abs(r["D_obs"] - r["D"]) > 0.05 * r["D"]:
        msg.append(f"实测插值 D={r['D_obs']:.0f} 天与拟合值差 >5%")
    return "；".join(msg)


def grade_sort_key(g):
    digits = "".join(ch for ch in g if ch.isdigit())
    return (int(digits) if digits else 10 ** 9, g)


# ---------------------------------------------------------------- 作图
def plot(data, fits, eas, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
    plt.rcParams["font.family"] = ["Microsoft YaHei", "PingFang SC", "Noto Sans CJK SC", "Source Han Sans SC",
                                   "WenQuanYi Zen Hei", "SimHei", "Arial Unicode MS", "sans-serif"]
    plt.rcParams["axes.unicode_minus"] = False
    INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
    SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]

    grades = sorted({g for g, _ in data}, key=grade_sort_key)
    temps = sorted({T for _, T in data})
    tcolor = {T: SERIES[i % len(SERIES)] for i, T in enumerate(temps)}
    n = len(grades)
    ncol = min(4, n)
    nrow = math.ceil(n / ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.0 * ncol, 3.4 * nrow + 0.9), dpi=160,
                             facecolor=SURF, squeeze=False)
    for ax in axes.flat[n:]:
        ax.set_visible(False)
    for ax, g in zip(axes.flat, grades):
        ax.set_facecolor(SURF)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        for sp in ("left", "bottom"):
            ax.spines[sp].set_color(GRID)
        ax.tick_params(colors=INK2, labelsize=8)
        ax.grid(True, color=GRID, lw=0.7)
        tmax = max(max(data[(gg, T)][0].max(), fits[(gg, T)]["D"] * 1.25 if (gg, T) in fits else 0)
                   for gg, T in data if gg == g)
        tt = np.linspace(0, tmax * 1.05, 400)
        lines = []
        for T in temps:
            if (g, T) not in data:
                continue
            t, y = data[(g, T)]
            c = tcolor[T]
            ax.plot(t, y, "o", ms=5, color=c, mec=SURF, mew=1, zorder=3)
            if (g, T) in fits:
                r = fits[(g, T)]
                ax.plot(tt, model(tt, r["f0"], r["lam"], r["beta"]), color=c, lw=1.8)
                lines.append(f"{T:g} °C：D={r['D']:.0f} 天  β={r['beta']:.2f}  f₀={100 * r['f0']:.1f}%  R²={r['r2']:.3f}")
        ax.axhline(80, color=INK2, lw=0.7, ls=(0, (2, 3)))
        ax.set_ylim(0, 102)
        ax.set_xlim(0, tmax * 1.05)
        title = g + (f"   Ea≈{eas[g]['Ea']:.0f} kJ/mol" if g in eas else "")
        ax.set_title(title, loc="left", fontsize=10.5, color=INK, pad=6)
        ax.text(0.98, 0.04, "\n".join(lines), transform=ax.transAxes, ha="right", va="bottom",
                fontsize=7.2, color=INK, bbox=dict(boxstyle="round,pad=0.3", fc=SURF, ec=GRID, lw=0.6))
        ax.set_xlabel("天数", color=INK2, fontsize=8.5)
        ax.set_ylabel("累计释放（%）", color=INK2, fontsize=8.5)
    if len(temps) > 1:
        from matplotlib.lines import Line2D
        handles = [Line2D([], [], color=tcolor[T], marker="o", lw=1.8, label=f"{T:g} °C") for T in temps]
        fig.legend(handles=handles, loc="upper right", ncol=len(temps), frameon=False, fontsize=9,
                   labelcolor=INK)
    fig.suptitle("包膜尿素静水释放：实测点（圆点）与 Weibull 拟合曲线（实线）", x=0.01, ha="left",
                 fontsize=12.5, color=INK, fontweight="bold")
    fig.text(0.01, 0.005, "F(t) = f₀ + (1 − f₀)·[1 − exp(−(t/λ)^β)]；D 为拟合曲线达到 80% 的天数；虚线为 80%。",
             fontsize=8, color=INK2)
    fig.tight_layout(rect=(0, 0.02, 1, 0.95))
    fig.savefig(path, facecolor=SURF)
    plt.close(fig)


# ---------------------------------------------------------------- 主程序
def main(argv=None):
    ap = argparse.ArgumentParser(description="拟合包膜尿素静水释放曲线（D、β、f₀、Ea）")
    ap.add_argument("csv", help="数据文件（CSV）")
    ap.add_argument("-o", "--out", help="输出目录（默认：数据文件旁的 fit_output/）")
    ap.add_argument("--fix-f0", type=float, default=None, metavar="百分数",
                    help="固定初期溶出 f₀（%%），只拟合 λ 与 β；早期取样少时建议使用")
    args = ap.parse_args(argv)

    src = Path(args.csv)
    out = Path(args.out) if args.out else src.parent / "fit_output"
    out.mkdir(parents=True, exist_ok=True)
    data = read_data(src)

    fits, errors = {}, {}
    for k, (t, y) in data.items():
        try:
            fits[k] = fit_one(t, y, args.fix_f0)
        except ValueError as e:
            errors[k] = str(e)

    # 多温度 → Arrhenius 表观活化能：1/λ ∝ exp(−Ea/RT)
    eas = {}
    by_grade = defaultdict(list)
    for (g, T), r in fits.items():
        by_grade[g].append((T, r))
    for g, lst in by_grade.items():
        if len({T for T, _ in lst}) >= 2:
            x = np.array([1 / (T + 273.15) for T, _ in lst])
            yv = np.array([math.log(1 / r["lam"]) for _, r in lst])
            slope, icpt = np.polyfit(x, yv, 1)
            Ea = -slope * R_GAS / 1000
            ts = sorted(T for T, _ in lst)
            betas = [r["beta"] for _, r in sorted(lst)]
            eas[g] = dict(Ea=Ea, temps=ts, betas=betas)

    rows = []
    for (g, T) in sorted(fits, key=lambda k: (grade_sort_key(k[0]), k[1])):
        r = fits[(g, T)]
        rows.append({
            "档位": g, "温度°C": f"{T:g}", "数据点数": r["n"],
            "D_拟合(天)": f"{r['D']:.1f}", "D_实测插值(天)": f"{r['D_obs']:.1f}" if r["D_obs"] else "",
            "β": f"{r['beta']:.3f}", "曲线类型": shape_label(r["beta"]),
            "f0(%)": f"{100 * r['f0']:.2f}", "λ(天)": f"{r['lam']:.2f}",
            "释放10%(天)": f"{r['t10']:.1f}", "释放50%(天)": f"{r['t50']:.1f}",
            "t50/t80": f"{r['t50'] / r['D']:.3f}", "日释放高峰(天)": f"{r['peak']:.1f}",
            "R²": f"{r['r2']:.4f}", "RMSE(百分点)": f"{r['rmse']:.2f}",
            "Ea(kJ/mol)": f"{eas[g]['Ea']:.1f}" if g in eas else "",
            "提示": notes(r),
        })
    if fits:
        with open(out / "fit_results.csv", "w", newline="", encoding="utf-8-sig") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
        plot(data, fits, eas, out / "fit_curves.png")

    # 终端汇总
    print(f"读取 {len(data)} 条曲线；拟合成功 {len(fits)} 条。\n")
    print(f"{'档位':<10}{'温度':>6}{'D(天)':>8}{'β':>7}{'f0%':>7}{'R²':>8}{'RMSE':>7}  类型   提示")
    for r in rows:
        print(f"{r['档位']:<10}{r['温度°C']:>6}{r['D_拟合(天)']:>8}{r['β']:>7}{r['f0(%)']:>7}{r['R²']:>8}"
              f"{r['RMSE(百分点)']:>7}  {r['曲线类型']:<5}  {r['提示']}")
    for k, e in errors.items():
        print(f"{k[0]}（{k[1]:g} °C）未拟合：{e}")
    if eas:
        print("\n表观活化能（按各温度拟合的 λ 做 Arrhenius 回归）：")
        for g, e in sorted(eas.items(), key=lambda kv: grade_sort_key(kv[0])):
            bs = "、".join(f"{T:g} °C β={b:.2f}" for T, b in zip(e["temps"], e["betas"]))
            spread = max(e["betas"]) - min(e["betas"])
            verdict = ("β 基本不随温度变化，符合模型“温度只改变快慢、不改变形状”的假设" if spread <= 0.3 else
                       "β 随温度变化明显，模型“温度只改变快慢”的假设不成立，需按温度分别处理")
            print(f"  {g}: Ea ≈ {e['Ea']:.0f} kJ/mol（{bs}；{verdict}）")
    b25 = [fits[k]["beta"] for k in fits if abs(k[1] - 25) < 0.5]
    lin = [b for b in b25 if b < 1.8]
    sty = [b for b in b25 if b >= 1.8]
    if lin:
        print(f"\n25 °C 线性型档位（{len(lin)} 个）β 中位数：{float(np.median(lin)):.2f}"
              f" —— 可填入 rice-crf-map/params.toml [release] 的 beta")
    if sty:
        print(f"25 °C S 型档位（{len(sty)} 个）β 中位数：{float(np.median(sty)):.2f}"
              f" —— 可填入 rice-crf-map/params.toml [scenarios] 的 stype_beta")
    if fits:
        print(f"\n结果表：{out / 'fit_results.csv'}\n对比图：{out / 'fit_curves.png'}")
    else:
        print("\n没有曲线拟合成功，未生成结果表和对比图。每条曲线至少需要 4 个天数 > 0 的数据点"
              "（用 --fix-f0 时 3 个）。")
        return 1


if __name__ == "__main__":
    sys.exit(main())
