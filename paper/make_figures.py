"""Journal-format figures for manuscript 1 (Field Crops Research / Elsevier artwork guidelines).

Widths: single column 90 mm, 1.5 column 140 mm, double column 190 mm. Sans-serif (Arial-metric
Liberation Sans), 6.5–8 pt text, vector PDF plus 600 dpi PNG for every figure.

Inputs are the repository's own outputs, so rerun the models first if parameters change:
    python rice-crf-map/build.py && python wheat-crf-map/build.py && python maize-crf-map/build.py
    python analysis/mixed_shapes.py && python analysis/beta_scan.py
    python paper/make_figures.py            # writes paper/figures/ and paper/tables/
"""
import csv
import json
import logging
import statistics as st
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "paper" / "figures"
TAB = ROOT / "paper" / "tables"
sys.path.insert(0, str(ROOT / "analysis"))
sys.path.insert(0, str(ROOT / "common"))

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import PolyCollection  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
MM = 1 / 25.4
W1, W15, W2 = 90 * MM, 140 * MM, 190 * MM
plt.rcParams.update({
    "font.family": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "font.size": 7, "axes.titlesize": 7.5, "axes.labelsize": 7, "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
    "legend.fontsize": 6.5, "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5, "lines.linewidth": 1.2, "axes.spines.top": False,
    "axes.spines.right": False, "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.dpi": 600,
    "legend.frameon": False, "axes.titlelocation": "left", "axes.titleweight": "bold",
})
INK, GREY, LIGHT = "#1a1a1a", "#6b6b6b", "#d9d9d9"
# colour-blind-safe categorical set (Okabe–Ito)
C = dict(rice="#0072B2", wheat_frost="#D55E00", wheat_other="#E69F00", maize_temp="#009E73",
         maize_trop="#CC79A7", urea="#56B4E9", lin="#0072B2", s="#D55E00", grey=GREY)
CROP_NAME = {"rice": "Rice", "wheat": "Wheat", "maize": "Maize"}


def save(fig, name):
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f"{name}.pdf")
    fig.savefig(OUT / f"{name}.png", dpi=600)
    plt.close(fig)
    print("wrote", name)


def panel(ax, letter, title=""):
    ax.set_title(f"{letter}  {title}".rstrip(), loc="left", fontsize=7.5, fontweight="bold", pad=4)


def load_map(crop):
    return json.load(open(ROOT / f"{crop}-crf-map" / "data.json", encoding="utf-8"))


def read_csv(path):
    with open(path, encoding="utf-8-sig") as fh:
        return list(csv.DictReader(fh))


# ---------------------------------------------------------------- Fig. 1 worked example
def fig1():
    import field_release as FR
    d = load_map("rice")
    z = next(z for z in d["zones"] if z["id"] == "cn-myz")
    s = z["seasons"][0]
    t = np.array(s["series"]["t"])
    fig, axes = plt.subplots(1, 3, figsize=(W2, 52 * MM), layout="constrained")
    # a: static water
    ax = axes[0]
    days = np.linspace(0, 160, 321)
    for beta, col, lab in ((1.3, C["lin"], "Linear, β = 1.3"), (2.5, C["s"], "S-shaped, β = 2.5")):
        ax.plot(days, 100 * FR.release(days, 90, 0.03, beta), color=col, label=lab)
    ax.axhline(80, color=LIGHT, lw=0.6, zorder=0)
    ax.axvline(90, color=LIGHT, lw=0.6, zorder=0)
    ax.text(93, 74, "80% at D = 90 d", color=GREY, fontsize=6, va="top")
    ax.set(xlabel="Days in water at 25 °C", ylabel="Cumulative N release (%)", ylim=(0, 102), xlim=(0, 160))
    ax.legend(loc="upper left")
    panel(ax, "a", "Release in static water")
    # b: field supply vs demand, linear recipe and S alternative
    ax = axes[1]
    sup_lin = np.sum([np.array(v) for v in s["series"]["comps"].values()], axis=0)
    sup_s = np.sum([np.array(v) for v in s["alt"]["s1"]["series_comps"].values()], axis=0)
    ax.plot(t, s["series"]["demand"], color=INK, lw=1.4, label="Crop demand on fertiliser N")
    ax.plot(t, sup_lin, color=C["lin"], label=f"Released: {rec_en(s['rec'])} ({s['total']} kg N)")
    ax.plot(t, sup_s, color=C["s"], label=f"Released: {rec_en(s['alt']['s1']['rec'], 'S')} ({s['alt']['s1']['total']} kg N)")
    for x, lab in ((s["tPI"], "PI"), (s["tH"], "Heading")):
        ax.axvline(x, color=LIGHT, lw=0.6, zorder=0)
        ax.text(x + 1.5, 4, lab, color=GREY, fontsize=6)
    ax.set(xlabel="Days after transplanting", ylabel="Cumulative N (kg ha$^{-1}$)", xlim=(0, s["days"]))
    ax.legend(loc="upper left", fontsize=6)
    panel(ax, "b", "Field release vs crop demand")
    # c: pool vs floor
    ax = axes[2]
    ax.fill_between(t, s["series"]["pool"], color=C["lin"], alpha=0.18, lw=0)
    ax.plot(t, s["series"]["pool"], color=C["lin"], label="Soil fertiliser-N pool, linear")
    ax.plot(t, s["alt"]["s1"]["series_pool"], color=C["s"], label="Soil fertiliser-N pool, S-shaped")
    ax.plot(t, s["series"]["floor"], color=INK, lw=0.9, ls=(0, (3, 2)), label="No-shortfall floor")
    ax.set(xlabel="Days after transplanting", ylabel="N in pool (kg ha$^{-1}$)", xlim=(0, s["days"]),
           ylim=(0, 1.55 * max(s["series"]["pool"])))
    ax.legend(loc="upper right", fontsize=6)
    panel(ax, "c", "Pool against the floor")
    fig.suptitle("Wuhan, single-season hybrid indica rice (normal year)", x=0.005, ha="left", fontsize=7, color=GREY)
    save(fig, "Fig1_worked_example")


def rec_en(rec, pfx="CR"):
    return " + ".join(("urea" if k == "urea" else f"{pfx}{k}") + f" {p}%" for k, p in rec)


# ---------------------------------------------------------------- Fig. 2 global maps
def world_polys():
    topo = json.load(open(ROOT / "rice-crf-map" / "data" / "countries-50m.json"))
    sx, sy = topo["transform"]["scale"]
    tx, ty = topo["transform"]["translate"]
    arcs = []
    for arc in topo["arcs"]:
        a = np.cumsum(np.array(arc, float), axis=0)
        arcs.append(np.column_stack([a[:, 0] * sx + tx, a[:, 1] * sy + ty]))

    def ring(idx):
        pts = []
        for i in idx:
            a = arcs[i] if i >= 0 else arcs[~i][::-1]
            pts.append(a if not pts else a[1:])
        return np.vstack(pts)

    polys = []
    for g in topo["objects"]["countries"]["geometries"]:
        if g["type"] == "Polygon":
            polys += [ring(r) for r in g["arcs"][:1]]
        elif g["type"] == "MultiPolygon":
            polys += [ring(p[0]) for p in g["arcs"]]
    # drop rings crossing the antimeridian, which draw as horizontal smears
    return [p for p in polys if np.ptp(p[:, 0]) < 300]


def fig2():
    polys = world_polys()
    fig, axes = plt.subplots(3, 1, figsize=(W2, 205 * MM), layout="constrained")
    cmap = plt.get_cmap("viridis_r")
    norm = matplotlib.colors.Normalize(30, 130)
    for ax, crop, letter in zip(axes, ("rice", "wheat", "maize"), "abc"):
        ax.add_collection(PolyCollection(polys, facecolor="#eeeeee", edgecolor="#bdbdbd", lw=0.2))
        d = load_map(crop)
        xs, ys, cs, ss_, n = [], [], [], [], 0
        for z in d["zones"]:
            for k, s in enumerate(z["seasons"]):
                cr = max(((k2, p) for k2, p in s["rec"] if k2 != "urea"), key=lambda r: r[1])[0]
                jitter = 1.2 * k
                xs.append(z["lon"] + jitter); ys.append(z["lat"] - jitter)
                cs.append(cr); ss_.append(s["total"]); n += 1
        sc = ax.scatter(xs, ys, c=cs, s=np.array(ss_) * 0.18, cmap=cmap, norm=norm, edgecolor=INK, linewidth=0.35,
                        zorder=3)
        ax.set(xlim=(-180, 180), ylim=(-50, 72), aspect="equal", xticks=[], yticks=[])
        for sp in ax.spines.values():
            sp.set_visible(False)
        panel(ax, letter, f"{CROP_NAME[crop]} ({len(d['zones'])} zones, {n} seasons)")
    cb = fig.colorbar(sc, ax=axes, orientation="horizontal", fraction=0.025, pad=0.01, aspect=50, shrink=0.6)
    cb.set_label("Main control period in the recommended blend, D (days at 25 °C)")
    handles = [Line2D([], [], marker="o", ls="", markerfacecolor="white", markeredgecolor=INK, markeredgewidth=0.35,
                      markersize=np.sqrt(v * 0.18), label=f"{v} kg N ha$^{{-1}}$") for v in (60, 150, 250)]
    axes[0].legend(handles=handles, loc="lower left", title="Total N", title_fontsize=6.5, fontsize=6)
    save(fig, "Fig2_global_recipes")


# ---------------------------------------------------------------- Fig. 3 grades vs shape
def fig3():
    fig, axes = plt.subplots(1, 2, figsize=(W15, 60 * MM), layout="constrained",
                             gridspec_kw=dict(width_ratios=[1, 1.5]))
    rng = np.random.default_rng(1)
    ax = axes[0]
    vals = []
    for crop in ("rice", "wheat", "maize"):
        d = load_map(crop)
        vals.append([100 * (1 - s["total"] / s["split_total"]) for z in d["zones"] for s in z["seasons"]])
    box(ax, vals, [C["rice"], C["wheat_other"], C["maize_temp"]], rng)
    ax.set_xticks(range(3), ["Rice", "Wheat", "Maize"])
    ax.set_ylabel("N saved vs. equally safe split urea (%)")
    ax.axhline(0, color=LIGHT, lw=0.6, zorder=0)
    panel(ax, "a", "Urea + one linear grade")
    ax = axes[1]
    rows = read_csv(ROOT / "analysis" / "mixed_shapes_results.csv")
    labels = ["2 linear\ngrades", "1 S-shaped\ngrade", "Linear +\nS-shaped"]
    keys = ["B", "C", "D"]
    pos, data, cols = [], [], []
    for j, crop in enumerate(("rice", "wheat")):
        r = [x for x in rows if x["crop"] == crop]
        for i, k in enumerate(keys):
            data.append([100 * (float(x["A"]) - float(x[k])) / float(x["A"]) for x in r if x[k]])
            pos.append(i * 2.6 + j * 1.0)
            cols.append(C["rice"] if crop == "rice" else C["wheat_other"])
    box(ax, data, cols, rng, positions=pos)
    ax.set_xticks([i * 2.6 + 0.5 for i in range(3)], labels)
    ax.axhline(0, color=LIGHT, lw=0.6, zorder=0)
    ax.set_ylabel("Further N saved vs. urea + one linear grade (%)")
    ax.legend(handles=[Line2D([], [], marker="s", ls="", color=C["rice"], label="Rice"),
                       Line2D([], [], marker="s", ls="", color=C["wheat_other"], label="Wheat")], loc="upper left")
    panel(ax, "b", "More grades or a different shape (β = 2.5)?")
    save(fig, "Fig3_grades_vs_shape")


def box(ax, data, cols, rng, positions=None):
    positions = list(range(len(data))) if positions is None else positions
    bp = ax.boxplot(data, positions=positions, widths=0.6, showfliers=False, patch_artist=True,
                    medianprops=dict(color=INK, lw=1), whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6),
                    boxprops=dict(lw=0.6))
    for patch, c in zip(bp["boxes"], cols):
        patch.set_facecolor(matplotlib.colors.to_rgba(c, 0.25))
        patch.set_edgecolor(c)
    for x, v, c in zip(positions, data, cols):
        ax.scatter(x + rng.uniform(-0.18, 0.18, len(v)), v, s=4, color=c, alpha=0.7, lw=0, zorder=3)


# ---------------------------------------------------------------- Fig. 4 β scan
GROUPS = [("Rice", lambda r: r["crop"] == "rice", C["rice"]),
          ("Wheat, frozen soil", lambda r: r["crop"] == "wheat" and r["frost"] == "True", C["wheat_frost"]),
          ("Wheat, other", lambda r: r["crop"] == "wheat" and r["frost"] == "False", C["wheat_other"]),
          ("Maize, temperate", lambda r: r["crop"] == "maize" and r["trop"] == "False", C["maize_temp"]),
          ("Maize, tropical", lambda r: r["crop"] == "maize" and r["trop"] == "True", C["maize_trop"])]


def beta_rows():
    rows = read_csv(ROOT / "analysis" / "beta_scan_results.csv")
    betas = sorted({float(r["beta"]) for r in rows})
    return rows, betas


def fig4():
    rows, betas = beta_rows()
    fig, axes = plt.subplots(1, 3, figsize=(W2, 58 * MM), layout="constrained")
    specs = [("N85_mean", "Yield-loss proxy at 85% N,\nmean of 5 conditions (%)", "a", "Equal N: average risk", np.mean),
             ("N85_worst", "Yield-loss proxy at 85% N,\nworst condition (%)", "b", "Equal N: worst case", np.mean),
             ("N_req_pct", "N needed for no shortfall\n(% of β = 1.3)", "c", "Equal yield: N needed", np.median)]
    for ax, (key, ylab, letter, title, agg) in zip(axes, specs):
        for name, sel, col in GROUPS:
            ys = [agg([float(r[key]) for r in rows if sel(r) and float(r["beta"]) == b]) for b in betas]
            ax.plot(betas, ys, "o-", color=col, ms=2.8, lw=1.1, label=name)
        ax.axvspan(2.0, 2.5, color="#f0f0f0", zorder=0, lw=0)
        ax.axvline(1.3, color=GREY, lw=0.6, ls=(0, (2, 2)))
        ax.set(xlabel="Release shape β", ylabel=ylab, xticks=[1, 2, 3, 4, 5])
        panel(ax, letter, title)
    axes[0].text(1.33, axes[0].get_ylim()[1] * 0.97, "current\nlinear", fontsize=6, color=GREY, va="top")
    axes[2].legend(loc="upper center", fontsize=6)
    save(fig, "Fig4_beta_scan")


# ---------------------------------------------------------------- Fig. 5 mechanism: robustness of high β
def fig5():
    import beta_scan as B
    crop, zi, si = "maize", 0, 0
    z, s, ss, P, system = B.setup(crop, zi, si)
    mod = B.WM
    betas = (1.3, 2.0, 4.0)        # 2.0 = optimum for maize in the β scan
    dTs = ((-3.0, "#0072B2", "3 °C cooler"), (0.0, INK, "Normal year"), (3.0, "#D55E00", "3 °C warmer"))
    fig, axes = plt.subplots(1, 3, figsize=(W2, 54 * MM), layout="constrained", sharey=True)
    for ax, beta, letter in zip(axes, betas, "abc"):
        tot, D, u, _ = B.design(crop, ss, P, system, beta, mod)
        for dT, col, lab in dTs:
            k = B.k_loss(crop, ss, P, dT)
            tau = B.field_tau(ss, system, dT)
            inc, _ = B.cr_inc(tau, D, beta)
            supply = tot * (u * mod.urea_increment(ss["n"], P) + (1 - u) * inc)
            pool, p, short = [], 0.0, 0.0
            for t in range(ss["n"]):
                p = p * (1 - k[t]) + supply[t]
                draw = ss["d"][t] / P["eta"]
                got = min(p, draw)
                p -= got
                short += ss["d"][t] - got * P["eta"]
                pool.append(p)
            ax.plot(np.arange(ss["n"]), pool, color=col, lw=1.1,
                    label=f"{lab}: shortfall {100 * short / ss['uptake_total']:.1f}%")
        for x, lab in ((ss["tPI"], "V10"), (ss["tH"], "Silking")):
            ax.axvline(x, color=LIGHT, lw=0.6, zorder=0)
            ax.text(x + 1.5, 1.5, lab, color=GREY, fontsize=6)
        ax.set(xlabel="Days after sowing", xlim=(0, ss["L"]))
        ax.legend(loc="upper right", fontsize=5.8, handlelength=1.2)
        panel(ax, letter, f"β = {beta}: urea {round(100 * u)}% + D{D}, {tot:.0f} kg N")
    axes[0].set_ylabel("Fertiliser-N pool (kg ha$^{-1}$)")
    fig.suptitle("NE China spring maize; each recipe designed to be safe within ±1.5 °C", x=0.005, ha="left",
                 fontsize=7, color=GREY)
    save(fig, "Fig5_robustness")


# ---------------------------------------------------------------- Fig. 6 environment
def fig6():
    rows, _ = beta_rows()
    best = {"Rice": 2.5, "Wheat": 2.5, "Maize": 2.0}
    sel = {"Rice": lambda r: r["crop"] == "rice", "Wheat": lambda r: r["crop"] == "wheat",
           "Maize": lambda r: r["crop"] == "maize"}
    metrics = [("N_req", "Total N"), ("NH3", "NH$_3$"), ("N2O", "N$_2$O"), ("LEACH", "Leaching +\nrunoff"),
               ("GHG", "GHG"), ("DMG", "Damage\ncost")]
    fig, ax = plt.subplots(figsize=(W15, 58 * MM), layout="constrained")
    width = 0.26
    for j, (crop, col) in enumerate(zip(best, (C["rice"], C["wheat_other"], C["maize_temp"]))):
        by = {}
        for r in rows:
            if sel[crop](r):
                by.setdefault((r["zone"], r["season"]), {})[float(r["beta"])] = r
        meds, lo, hi = [], [], []
        for key, _ in metrics:
            ch = [100 * (float(v[best[crop]][key]) / float(v[1.3][key]) - 1) for v in by.values()
                  if float(v[1.3][key]) > 0]
            q = np.percentile(ch, [25, 50, 75])
            meds.append(q[1]); lo.append(q[1] - q[0]); hi.append(q[2] - q[1])
        x = np.arange(len(metrics)) + (j - 1) * width
        ax.bar(x, meds, width * 0.92, color=col, label=f"{crop} (β = {best[crop]})")
        ax.errorbar(x, meds, yerr=[lo, hi], fmt="none", ecolor=INK, elinewidth=0.6, capsize=1.5)
    ax.axhline(0, color=INK, lw=0.6)
    ax.set_xticks(range(len(metrics)), [m for _, m in metrics])
    ax.set_ylabel("Change vs. β = 1.3 at equal yield (%)\nmedian and interquartile range")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=3)
    ax.set_ylim(top=1)
    panel(ax, "", "Environmental effects of the optimal release shape")
    save(fig, "Fig6_environment")


# ---------------------------------------------------------------- Fig. S1 validation demo (simulated data)
def figS1():
    rows = read_csv(ROOT / "validation" / "example" / "output" / "observed_vs_predicted.csv")
    rows = [r for r in rows if int(r["天数"]) > 0]
    fig, axes = plt.subplots(1, 2, figsize=(W15, 70 * MM), layout="constrained")
    for ax, (key, title), letter in zip(axes, (("默认", "Literature parameters"), ("留一站点标定", "Leave-one-site-out calibration")), "ab"):
        for run, col, lab in (("A", C["lin"], "Run A: measured soil T and moisture"),
                              ("B", C["s"], "Run B: air temperature and rainfall")):
            x = [float(r["实测均值%"]) for r in rows]
            y = [float(r[f"{key}_{run}%"]) for r in rows]
            rmse = float(np.sqrt(np.mean((np.array(y) - np.array(x)) ** 2)))
            ax.scatter(x, y, s=5, color=col, alpha=0.75, lw=0, label=f"{lab} (RMSE {rmse:.1f})")
        ax.plot([0, 100], [0, 100], color=INK, lw=0.6)
        for dd in (-8, 8):
            ax.plot([0, 100], [dd, 100 + dd], color=GREY, lw=0.5, ls=(0, (2, 2)))
        ax.set(xlim=(0, 100), ylim=(0, 100), aspect="equal", xlabel="Observed cumulative release (%)",
               ylabel="Predicted cumulative release (%)")
        ax.legend(loc="lower right", fontsize=5.8, handletextpad=0.2)
        panel(ax, letter, title)
    fig.suptitle("Simulated burial-bag data (4 synthetic sites) demonstrating the validation workflow",
                 x=0.005, ha="left", fontsize=7, color=GREY)
    save(fig, "FigS1_validation_demo")


# ---------------------------------------------------------------- tables
def tables():
    TAB.mkdir(parents=True, exist_ok=True)
    rows, betas = beta_rows()
    out = []
    for crop in ("rice", "wheat", "maize"):
        d = load_map(crop)
        S = [s for z in d["zones"] for s in z["seasons"]]
        Ds = [k for s in S for k, _ in s["rec"] if k != "urea"]
        ureas = [dict(s["rec"]).get("urea", 0) for s in S]
        rb = [r for r in rows if r["crop"] == crop]
        opt = min(betas, key=lambda b: st.mean(float(r["N85_mean"]) for r in rb if float(r["beta"]) == b))
        nmin = min(betas, key=lambda b: st.median(float(r["N_req_pct"]) for r in rb if float(r["beta"]) == b))
        out.append({
            "Crop": CROP_NAME[crop], "Zones": len(d["zones"]), "Seasons": len(S),
            "Control periods used (d)": f"{min(Ds)}–{max(Ds)}", "Urea share (%)": f"{min(ureas)}–{max(ureas)}",
            "Median total N (kg/ha)": round(st.median(s["total"] for s in S)),
            "N saved vs split urea (median %)": round(100 * st.median(1 - s["total"] / s["split_total"] for s in S), 1),
            "N saved by S-type single grade (median %)": round(100 * st.median(1 - s["alt"]["s1"]["total"] / s["total"] for s in S), 1),
            "β with lowest mean yield-loss proxy": opt, "β with lowest N need": nmin,
        })
    with open(TAB / "Table2_summary_by_crop.csv", "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]))
        w.writeheader()
        w.writerows(out)
    print("wrote Table2_summary_by_crop.csv")
    for r in out:
        print(r)


if __name__ == "__main__":
    for f in (fig1, fig2, fig3, fig4, fig5, fig6, figS1, tables):
        f()
