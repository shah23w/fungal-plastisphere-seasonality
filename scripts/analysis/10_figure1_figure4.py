"""Figure 1, the seasonal and thermal variation in fungal enrichment, and Figure 4,
the cross-study forest plot.

    python3 10_figure1_figure4.py [input_dir] [out_dir]
"""
import csv, math, os, sys
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
D = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "..", "tables")
OUT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, "..", "..", "figures")
os.makedirs(OUT, exist_ok=True)


def rd(f):
    return list(csv.DictReader(open(os.path.join(D, f)), delimiter="\t"))


src = rd("TableS1_blocks.tsv")
gen = rd("TableS15_fungal_genera.tsv")
spc = rd("TableS16_species_coverage.tsv")
cs = rd("TableS3_crossstudy.tsv")
SEASF = os.path.join(D, "crossstudy_season.tsv")
SEAS = {}
for r in csv.DictReader(open(SEASF), delimiter="\t"):
    SEAS[r["study"]] = (r["season"], r["multi_season"])
COVNAME = {"water_temp_C": "Water temperature", "precip_mm_day": "Precipitation",
           "velocity_m_s": "Water velocity", "turbidity_NTU": "Turbidity",
           "salinity_psu": "Salinity", "pH": "pH", "DO_mg_L": "Dissolved oxygen",
           "COD_mg_L": "Chemical oxygen demand", "TN_ug_L": "Total nitrogen",
           "TP_ug_L": "Total phosphorus", "SO4_ug_L": "Sulfate",
           "total_MPs_per_L": "Microplastic concentration"}
cov = []
for r in csv.reader(open(os.path.join(D, "TableS9_covariates.tsv")), delimiter="\t"):
    if len(r) >= 7 and r[0] == "univariate_bound":
        try:
            cov.append((r[1], float(r[2]), float(r[3]), float(r[5])))
        except ValueError:
            pass
    elif len(r) == 6 and r[0] not in ("covariate", "") and r[1]:
        try:
            cov.append((r[0], float(r[1]), float(r[2]), float(r[4])))
        except ValueError:
            pass

SE = ["Winter", "Spring", "Summer", "Autumn"]
sites = sorted({r["site"] for r in src})
MP_C, NP_C, RW_C = "#c1443a", "#2f7d6f", "#4a6fa5"
INK, GREY = "#22252a", "#b9bcc2"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7.5,
                     "axes.linewidth": 0.7, "axes.spines.top": False,
                     "axes.spines.right": False, "legend.frameon": False})


def hl(v):
    w = [(v[i] + v[j]) / 2 for i in range(len(v)) for j in range(i, len(v))]
    return float(np.median(w))


def srank_p(v):
    v = [x for x in v if x != 0]
    return stats.wilcoxon(v, method="exact").pvalue if len(v) > 1 else 1.0


fig = plt.figure(figsize=(7.4, 7.3))
gs = gridspec.GridSpec(3, 2, figure=fig, hspace=0.46, wspace=0.52,
                       left=0.085, right=0.975, top=0.975, bottom=0.045,
                       height_ratios=[1, 1, 1])

A = fig.add_subplot(gs[1, 0])
sl = []
for s in sites:
    g = sorted((r for r in src if r["site"] == s), key=lambda r: float(r["temp"]))
    x = [float(r["temp"]) for r in g]
    y = [float(r["d"]) for r in g]
    _lr = stats.linregress(x, y)
    A.plot(x, [_lr.intercept + _lr.slope * _v for _v in x], "-", color=INK,
           lw=0.9, alpha=0.55, zorder=3, solid_capstyle="round")
    A.plot(x, y, "-", color=GREY, lw=0.5, alpha=0.65, zorder=1)
    sl.append(_lr.slope)
A.scatter([float(r["temp"]) for r in src], [float(r["d"]) for r in src],
          s=9, c=MP_C, lw=0, alpha=0.75, zorder=2)
xs = np.array([min(float(r["temp"]) for r in src), max(float(r["temp"]) for r in src)])
m = float(np.mean(sl))
gmv = float(np.mean([float(r["d"]) for r in src]))
gt = float(np.mean([float(r["temp"]) for r in src]))
A.axhline(0, c=GREY, lw=0.7, ls="--", zorder=0)
tt = stats.ttest_1samp(sl, 0.0)
A.legend(handles=[
    Line2D([], [], ls="-", lw=0.9, c=INK, alpha=0.55, label="Per-site fit"),
    Line2D([], [], marker="o", ls="", ms=3, mfc=MP_C, mec="none",
           label="Site $\\times$ season block")],
    fontsize=6, loc="upper left", frameon=False,
    handletextpad=0.5, borderpad=0.2, labelspacing=0.25)
A.set_xlabel("Water temperature (°C)")
A.set_ylabel("Fungal enrichment, MP − NP (log$_2$)")
A.set_title("c", loc="left",
            fontweight="bold", fontsize=11)

B = fig.add_subplot(gs[0, 0])
raw_p = {z: srank_p([float(r["d"]) for r in src if r["season"] == z]) for z in SE}
adj_p, _run = {}, 0.0
for _i, _z in enumerate(sorted(SE, key=lambda q: raw_p[q])):
    _run = max(_run, (len(SE) - _i) * raw_p[_z])
    adj_p[_z] = min(1.0, _run)
for i, s in enumerate(SE):
    v = [float(r["d"]) for r in src if r["season"] == s]
    xj = i + (np.arange(len(v)) - len(v) / 2) * 0.035
    B.scatter(xj, v, s=11, c=MP_C, lw=0, alpha=0.7)
    h = hl(v)
    p = adj_p[s]
    B.plot([i - 0.26, i + 0.26], [h, h], "-", c=INK, lw=2)
    star = ("***" if p < 0.001 else "**" if p < 0.01
            else "*" if p < 0.05 else "ns")
    B.text(i, max(v) + 0.9, star, ha="center",
           va="bottom", fontsize=(9 if star != "ns" else 6.6),
           c=(INK if p < 0.05 else "#8a8f98"))
B.axhline(0, c=GREY, lw=0.7, ls="--")
B.set_xticks(range(4))
B.set_xticklabels(SE, fontsize=7)
B.set_ylabel("Fungal enrichment, MP − NP (log$_2$)")
B.set_title("a", loc="left",
            fontweight="bold", fontsize=11)

C = fig.add_subplot(gs[0, 1])
for key, col, lab in (("log2_MP", MP_C, "microplastic"),
                      ("log2_NP", NP_C, "natural particle"),
                      ("log2_RW", RW_C, "river water")):
    xs2, ys = [], []
    for s in SE:
        g = [r for r in src if r["season"] == s and r.get(key)]
        if not g:
            continue
        xs2.append(SE.index(s))
        ys.append(float(np.median([float(r[key]) for r in g])))
    C.plot(xs2, ys, "-o", ms=4.5, c=col, lw=1.5, label=lab)
C.set_xticks(range(4)); C.set_xticklabels(SE, fontsize=7); C.set_xlim(-0.35, 3.35)
C.set_ylabel("Median fungal read fraction (log$_2$)")
C.set_title("b", loc="left",
            fontweight="bold", fontsize=11)
C.legend(["Microplastic", "Natural particle", "River water"],
         fontsize=6, loc="upper left", frameon=False,
         handlelength=1.3, borderpad=0.2, labelspacing=0.25)

Dx = fig.add_subplot(gs[2, 0])
agg = {}
for r in gen:
    k = r["genus"]
    a = 0 if r["arm"] == "plastic" else 1
    agg.setdefault(k, [0, 0])[a] += int(r["reads"])
top = sorted(agg, key=lambda g: -agg[g][0])[:8]
yv = np.arange(len(top))
Dx.barh(yv - 0.19, [max(agg[g][0], 0.6) for g in top], 0.38, color=MP_C,
        label="Microplastic")
Dx.barh(yv + 0.19, [max(agg[g][1], 0.6) for g in top], 0.38, color=NP_C,
        label="Natural particle")
Dx.set_yticks(yv)
Dx.set_yticklabels(top, fontstyle="italic", fontsize=7)
Dx.invert_yaxis()
Dx.set_xscale("log")
Dx.set_xlim(0.6, None)
Dx.tick_params(axis="x", which="minor", length=0)
Dx.set_xlabel("Fungal reads")
Dx.set_title("e", loc="left",
             fontweight="bold", fontsize=11)
Dx.legend(fontsize=6, loc="lower right", frameon=False,
          handlelength=1.1, borderpad=0.2, labelspacing=0.25)

E = fig.add_subplot(gs[2, 1])
for r in spc:
    o = float(r["index_coverage"]) * 100
    e = float(r["poisson_expected"]) * 100
    hit = r["species"].startswith("Pichia kudriav")
    E.scatter([e], [o], s=(26 if hit else 13), c=(MP_C if hit else GREY),
              lw=(0.5 if hit else 0), edgecolors=INK, zorder=(3 if hit else 2))
E.plot([0, 100], [0, 100], ls=":", c=INK, lw=0.8)
E.set_xlim(0, 105)
E.set_ylim(0, 105)
E.set_xlabel("Expected coverage (%)")
E.set_ylabel("Observed coverage (%)")
E.set_title("f", loc="left", fontweight="bold", fontsize=11)
best = max(spc, key=lambda r: float(r["index_coverage"]))
E.legend(handles=[
    Line2D([], [], marker="o", ls="", ms=5, mfc=MP_C, mec=INK, mew=0.5,
           label="$P.\\ kudriavzevii$"),
    Line2D([], [], marker="o", ls="", ms=4, mfc=GREY, mec="none",
           label="Other species")],
    fontsize=6, loc="lower right", frameon=False,
    handletextpad=0.4, borderpad=0.2, labelspacing=0.25)

F = fig.add_subplot(gs[1, 1])
RAIN = [("antecedent_rain_1d_mm", 2.8750, 15.7520, 0.5911),
        ("antecedent_rain_3d_mm", 1.1747, 1.1279, 0.007296),
        ("antecedent_rain_7d_mm", 1.5837, 1.5584, 0.008524)]
_have = {c[0] for c in cov}
for _r in RAIN:
    if _r[0] not in _have:
        cov.append(_r); print("bounds file lacked %s - added from Table S9" % _r[0])
for _k, _v in (("antecedent_rain_1d_mm", "Antecedent rainfall (1 d)"),
               ("antecedent_rain_3d_mm", "Antecedent rainfall (3 d)"),
               ("antecedent_rain_7d_mm", "Antecedent rainfall (7 d)")):
    COVNAME[_k] = _v
rat = sorted([(n, obs / md, p) for (n, obs, md, p) in cov if md > 0],
             key=lambda c: abs(c[1]))
yv = np.arange(len(rat))
for i, (n, r, p) in enumerate(rat):
    hot = n == "water_temp_C"
    F.plot([0, r], [i, i], "-", c=GREY, lw=1.0, alpha=0.8, zorder=1)
    F.plot([r], [i], "o", ms=(6 if hot else 4.2),
           c=(MP_C if hot else RW_C if n.startswith("antecedent_rain") else INK),
           zorder=3)
F.axvline(0, c=GREY, lw=0.7, ls="--")
for _b in (-1.0, 1.0):
    F.axvline(_b, c=INK, lw=0.6, ls=":", alpha=0.85, zorder=0)
F.set_yticks(yv)
F.set_yticklabels([COVNAME.get(c[0], c[0]) for c in rat], fontsize=6.4)
F.legend(handles=[
    Line2D([], [], marker="o", ls="", ms=5, c=MP_C, label="Water temperature"),
    Line2D([], [], marker="o", ls="", ms=4, c=RW_C, label="Antecedent rainfall"),
    Line2D([], [], marker="o", ls="", ms=4, c=INK, label="Other variable"),
    ],
    fontsize=5.6, loc="upper left", frameon=False,
    handletextpad=0.5, borderpad=0.2, labelspacing=0.28)
F.set_xlim(-1.35, 1.15)
F.set_xticks([-1.0, -0.5, 0.0, 0.5, 1.0])
F.set_xlabel("Observed slope / detectable bound")
F.set_title("d", loc="left",
            fontweight="bold", fontsize=11)

fig4 = plt.figure(figsize=(7.0, 4.2))
G = fig4.add_subplot(111)
def dl_pool(rws):
    y = np.array([float(r["diff_log2"]) for r in rws])
    se = np.array([float(r["se"]) for r in rws])
    w = 1.0 / se ** 2
    fe = float((w * y).sum() / w.sum())
    Q = float((w * (y - fe) ** 2).sum()); dfq = len(y) - 1
    C = float(w.sum() - (w ** 2).sum() / w.sum())
    t2 = max(0.0, (Q - dfq) / C)
    i2 = (max(0.0, (Q - dfq) / Q) * 100.0) if Q > 0 else 0.0
    wr = 1.0 / (se ** 2 + t2)
    est = float((wr * y).sum() / wr.sum())
    return est, math.sqrt(1.0 / wr.sum()), i2
POOLED_OUT = []
POOL_ROWS = []
rows = sorted([r for r in cs if r.get("system", "").strip()],
              key=lambda r: float(r["diff_log2"]))
yv = np.arange(len(rows))
for i, r in enumerate(rows):
    lo, hi = float(r["ci_lo"]), float(r["ci_hi"])
    pooled_ok = int(r["n_surface_control"]) >= 3
    c = MP_C if "Hong Kong" in r["system"] else (INK if pooled_ok else GREY)
    G.plot([lo, hi], [i, i], "-", c=c, lw=1.2)
    G.plot([float(r["diff_log2"])], [i], "o", ms=4.5, c=c)
_pa = [r for r in rows if int(r["n_surface_control"]) >= 3]
_pe = [r for r in _pa if "Hong Kong" not in r["system"]]
# Both pooled diamonds are drawn in the figure's own ink. They are distinguished
# from the study rows by shape and by their own y-axis labels, so no extra colour
# is needed and none is introduced.
POOL_EXT = POOL_ALL = INK
for _j, (_rw, _lb, _cl) in enumerate(((_pe, "Random effects, external only", POOL_EXT),
                                      (_pa, "Random effects, all nine", POOL_ALL))):
    _e, _se, _i2 = dl_pool(_rw)
    _yy = -1.15 - 0.95 * _j
    _lo, _hi = _e - 1.96 * _se, _e + 1.96 * _se
    G.add_patch(plt.Polygon([[_lo, _yy], [_e, _yy + 0.3], [_hi, _yy], [_e, _yy - 0.3]],
                            closed=True, facecolor=_cl, edgecolor=_cl, lw=0.8, zorder=4))
    POOL_ROWS.append((_yy, "%s ($I^2$ = %.0f%%)" % (_lb, _i2)))
    POOLED_OUT.append((_lb, _e, 2 ** _e, 2 ** _lo, 2 ** _hi, _i2, len(_rw)))
G.axvline(0, c=GREY, lw=0.8, ls="--")
G.set_yticks(yv)
G.set_yticks(list(yv) + [p[0] for p in POOL_ROWS])
G.set_yticklabels([r["study"] for r in rows] + [p[1] for p in POOL_ROWS], fontsize=6.6)
G.axhline(-0.62, c=INK, lw=0.6)
G.set_ylim(-1.15 - 0.95 * 2 - 0.65, len(rows) - 0.45)
G.legend(handles=[
    Line2D([], [], marker="o", ls="-", ms=4.5, c=MP_C, label="This dataset"),
    Line2D([], [], marker="o", ls="-", ms=4.5, c=INK,
           label="Other study, $n$ $\\geq$ 3 controls"),
    Line2D([], [], marker="o", ls="-", ms=4.5, c=GREY,
           label="Other study, $n$ < 3 controls"),
    ],
    fontsize=6.6, loc="upper left", frameon=False,
    title_fontsize=6.6, alignment="left",
    handletextpad=0.5, borderpad=0.2, labelspacing=0.3)
G.set_xlabel("Fungal fraction, plastic − control (log$_2$)")

ns = sum(1 for r in rows if SEAS.get(r["study"], ("?", "no"))[1] == "no")
nsum = sum(1 for r in rows if SEAS.get(r["study"], ("?", "no"))[1] == "no"
           and "Summer" in SEAS.get(r["study"], ("?", "no"))[0])

fig4.subplots_adjust(left=0.17, right=0.97, top=0.93, bottom=0.17)
for ext, kw in (("png", {"dpi": 400}), ("pdf", {}),
                ("tiff", {"dpi": 400, "pil_kwargs": {"compression": "tiff_lzw"}})):
    fig.savefig(os.path.join(OUT, "Figure1.%s" % ext), bbox_inches="tight", **kw)
    fig4.savefig(os.path.join(OUT, "Figure4.%s" % ext), bbox_inches="tight", **kw)
print("wrote Figure1.png/.pdf/.tiff  (6 panels)")
print("wrote Figure4.png/.pdf/.tiff  (cross-study forest, was panel 1g)")
print("")
print("=== VERIFICATION (must match the manuscript) ===")
print("c  mean within-site slope %+.4f log2/C, %d/%d sites positive, P=%.4f"
      % (m, sum(1 for v in sl if v > 0), len(sl), tt.pvalue))
print("c  global regression line drawn: NO (per-site segments only)")
print("a  Holm-adjusted P and plotted symbol by season:")
for _z in SE:
    _p = adj_p[_z]
    _st = ("***" if _p < 0.001 else "**" if _p < 0.01 else "*" if _p < 0.05 else "ns")
    print("      %-7s nominal %.4f  Holm %.4f  -> %s" % (_z, raw_p[_z], _p, _st))
print("      expected: Winter ns (0.0905), Summer * (0.0105), Spring/Autumn ns")
print("a  bars are Hodges-Lehmann estimates:")
for _z in SE:
    _h = hl([float(r["d"]) for r in src if r["season"] == _z])
    print("      %-7s HL %+.4f log2 (%.2f-fold)" % (_z, _h, 2 ** _h))

print("d  observed/bound, ordered by that ratio:")
for _n, _r, _p in rat:
    print("      %-28s %6.3f  P=%.4f" % (COVNAME.get(_n, _n), _r, _p))
print("Fig4  pooled random-effects estimates:")
for _lb, _e, _f, _fl, _fh, _i2, _k in POOLED_OUT:
    print("      %-32s k=%d  %.3f-fold (95%% CI %.3f-%.3f)  I2=%.0f%%"
          % (_lb, _k, _f, _fl, _fh, _i2))
print("      expected: external only 0.834, all nine 0.893, I2 78/79%")
print("Fig4  row labels are accessions: %s ..." % ", ".join(r["study"] for r in rows[:3]))
print("")
print("Figure 1e: top genus %s" % top[0])
print("Figure 1f: max coverage %s %.1f%%" % (best["species"], float(best["index_coverage"]) * 100))
print("Figure 4: %d studies, %d of %d single-timepoint sampled summer" % (len(rows), nsum, ns))

