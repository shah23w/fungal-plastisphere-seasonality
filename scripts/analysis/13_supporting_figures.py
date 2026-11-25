"""Supporting figures S1 to S12, from the deposited tables.

    python3 13_supporting_figures.py [table_dir] [out_dir]

Exit status is 0 only if every plotted value matches the value in the paper.
"""
import csv, math, os, sys
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
D = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "..", "tables")
OUT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, "..", "..", "figures")
os.makedirs(OUT, exist_ok=True)


def T(n):
    p = os.path.join(D, n)
    if not os.path.exists(p):
        raise SystemExit("missing %s" % p)
    return list(csv.DictReader(open(p), delimiter="\t"))


MP, NP, RW = "#c1443a", "#2f7d6f", "#4a6fa5"
INK, GREY, AMBER = "#22252a", "#b9bcc2", "#d9a43b"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7.5,
                     "axes.linewidth": 0.7, "axes.spines.top": False,
                     "axes.spines.right": False, "legend.frameon": False})
SE = ["Winter", "Spring", "Summer", "Autumn"]
DEG = "\u00b0"
FUNGI, BACT = "#9b3fb5", "#2f7d6f"
fails = []
nchk = [0]


def num(fmt, v):
    """Format a number with a typographic minus, so in-plot text matches the
    minus sign matplotlib draws on the axes."""
    return (fmt % v).replace("-", "\u2212")


def chk(name, got, want, tol=0.0):
    nchk[0] += 1
    ok = (got == want) if isinstance(want, str) else abs(got - want) <= tol
    print("   %-52s %-14s want %-14s %s"
          % (name, "%.4f" % got if not isinstance(got, str) else got,
             "%.4f" % want if not isinstance(want, str) else want,
             "OK" if ok else "MISMATCH"))
    if not ok:
        fails.append(name)


def save(fig, stem):
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(OUT, "%s.%s" % (stem, ext)),
                    dpi=400 if ext == "png" else None, bbox_inches="tight")
    plt.close(fig)
    print("wrote %s" % stem)



# ================================================================ S1, S2, S4
# ================================================ S8, S10, S11, S12
# ------------------------------------------------------------------- FigS1
cr = {r["season"]: r for r in T("TableS5_classification_rate.tsv")}
db = {r["quantity"]: float(r["value"]) for r in T("TableS6_denominator_bias.tsv")}
fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.1))
A = ax[0]
w = 0.28
pts1 = [float(cr[s]["diff_of_medians_points"]) for s in SE]
pts2 = [float(cr[s]["mean_of_block_diffs_points"]) for s in SE]
lg = [float(cr[s]["mean_of_block_log2_rate_ratio"]) for s in SE]
x = np.arange(4)
A.bar(x - w, pts1, w, color=GREY, label="Difference of medians (left axis)")
A.bar(x, pts2, w, color=NP, label="Mean of block differences (left axis)")
A.axhline(0, c=INK, lw=0.7)
A.set_xticks(x); A.set_xticklabels(SE)
A.set_ylabel("MP minus NP, classified reads\n(percentage points)")
A.set_title("a", loc="left", fontweight="bold", fontsize=11)
A2 = A.twinx()
A2.bar(x + w, lg, w, color=MP, label="Mean block log$_2$ ratio (right axis)")
# All text stays black. The two scales are told apart by which side they sit on
# and by the legend, not by the colour of the axis label.
A2.set_ylabel("MP minus NP, classified reads\n(log$_2$ ratio)", color=INK)
A2.tick_params(axis="y", colors=INK)
A2.spines["right"].set_visible(True)
A2.spines["right"].set_color(INK)
A2.spines["top"].set_visible(False)
lim = max(abs(np.array(lg))) * 1.45
A2.set_ylim(-lim, lim)
lim1 = max(max(abs(np.array(pts1))), max(abs(np.array(pts2)))) * 1.45
A.set_ylim(-lim1, lim1)
h1, l1 = A.get_legend_handles_labels()
h2, l2 = A2.get_legend_handles_labels()
A.legend(h1 + h2, l1 + l2, fontsize=6.0, loc="lower left", handletextpad=0.5,
         labelspacing=0.25, borderpad=0.2)
B = ax[1]
B.errorbar([0], [db["withinsite_slope_dlog2rate_per_C"]],
           yerr=[[db["withinsite_slope_dlog2rate_per_C"] - db["withinsite_slope_CI_low"]],
                 [db["withinsite_slope_CI_high"] - db["withinsite_slope_dlog2rate_per_C"]]],
           fmt="o", ms=7, c=MP, lw=1.3, capsize=5)
B.axhline(0, c=INK, lw=0.9, ls="--")
B.set_xlim(-0.6, 0.6); B.set_xticks([])
B.set_ylabel("Temperature slope of the classification-rate\ncontrast "
             "(log$_2$ units °C$^{-1}$)")
B.text(0.05, 0.05, "P = %.3f" % db["withinsite_slope_P"], transform=B.transAxes,
       ha="left", fontsize=7.5, c=INK)
B.set_title("b", loc="left", fontweight="bold", fontsize=11)
fig.tight_layout(pad=1.0)
save(fig, "FigS1_denominator")
print("=== FigS1 ===")
chk("winter mean block log2 ratio", lg[0], 0.5142, 1e-4)
chk("autumn difference of medians (points)", pts1[3], -2.905, 0.002)
chk("slope of the classification-rate contrast", db["withinsite_slope_dlog2rate_per_C"],
    -0.0488, 1e-4)

# ------------------------------------------------------------------- FigS2
rows = T("TableS11_env_correlation.tsv")
ORDER = ["Water temperature", "Precipitation", "Water velocity", "Turbidity",
         "Salinity", "pH", "Dissolved oxygen", "Chemical oxygen demand",
         "Total nitrogen", "Total phosphorus", "Sulfate",
         "Microplastic concentration", "Altitude", "River length",
         "Upstream area", "Upstream population"]
# Table S11 named two site characteristics "Catchment area" and "Catchment
# population" while the main text calls them upstream area and upstream
# population. The table is rewritten so the deposit and the text agree.
RENAME = {"Catchment area": "Upstream area",
          "Catchment population": "Upstream population"}
if any(RENAME.get(r["variable_i"]) or RENAME.get(r["variable_j"]) for r in rows):
    for r in rows:
        r["variable_i"] = RENAME.get(r["variable_i"], r["variable_i"])
        r["variable_j"] = RENAME.get(r["variable_j"], r["variable_j"])
    with open(os.path.join(D, "TableS11_env_correlation.tsv"), "w", newline="") as fh:
        wtr = csv.DictWriter(fh, delimiter="\t", lineterminator="\n",
                             fieldnames=list(rows[0].keys()))
        wtr.writeheader()
        for r in rows:
            wtr.writerow(r)
    print("Table S11 rewritten with upstream area and upstream population")
rho = {}
pv = {}
for r in rows:
    a = RENAME.get(r["variable_i"], r["variable_i"])
    b = RENAME.get(r["variable_j"], r["variable_j"])
    rho[(a, b)] = rho[(b, a)] = float(r["spearman_rho"])
    pv[(a, b)] = pv[(b, a)] = float(r["P"])
n = len(ORDER)
M = np.full((n, n), np.nan)
for i, a in enumerate(ORDER):
    for j, b in enumerate(ORDER):
        if i == j:
            M[i, j] = 1.0
        elif j < i:
            M[i, j] = rho[(a, b)]
fig, A = plt.subplots(figsize=(7.4, 6.4))
im = A.imshow(M, cmap="RdBu_r", vmin=-1, vmax=1)
for i, a in enumerate(ORDER):
    for j, b in enumerate(ORDER):
        if j < i:
            p = pv[(a, b)]
            s = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""
            if s:
                A.text(j, i, s, ha="center", va="center", fontsize=6.5,
                       c="white" if abs(M[i, j]) > 0.45 else INK)
A.set_xticks(range(n)); A.set_xticklabels(ORDER, rotation=55, ha="right", fontsize=7)
A.set_yticks(range(n)); A.set_yticklabels(ORDER, fontsize=7)
A.set_xticks(np.arange(-0.5, n, 1), minor=True)
A.set_yticks(np.arange(-0.5, n, 1), minor=True)
A.grid(which="minor", color="white", lw=1.2)
A.tick_params(which="minor", length=0)
for sp in A.spines.values():
    sp.set_visible(False)
cb = fig.colorbar(im, ax=A, fraction=0.030, pad=0.02)
cb.set_label("Spearman's rho")
cb.outline.set_visible(False)
fig.tight_layout(pad=0.8)
save(fig, "FigS2_env_correlation")
print("=== FigS2 ===")
chk("temperature vs precipitation rho", rho[("Water temperature", "Precipitation")],
    0.4168, 1e-4)
chk("temperature vs upstream population rho",
    rho[("Water temperature", "Upstream population")], 0.1580, 2e-3)
chk("no |rho| with temperature reaches 0.7",
    float(max(abs(rho[("Water temperature", v)]) for v in ORDER[1:])), 0.4168, 1e-4)
chk("banned label absent",
    "absent" if not any("Catchment" in v for v in ORDER) else "present", "absent")

# ------------------------------------------------------------------- FigS4
g = T("TableS15_fungal_genera.tsv")
g = [r for r in g if int(r["reads"]) > 0 and float(r["reads_per_distinct_minimizer"]) > 0]
X = np.array([math.log10(int(r["reads"])) for r in g])
Y = np.array([math.log10(float(r["reads_per_distinct_minimizer"])) for r in g])
lr = stats.linregress(X, Y)
resid = Y - (lr.intercept + lr.slope * X)
zz = resid / resid.std(ddof=1)
top = int(np.argmax([float(r["reads_per_distinct_minimizer"]) for r in g]))
fig, A = plt.subplots(figsize=(6.2, 4.3))
A.scatter(10 ** X, 10 ** Y, s=12, c=GREY, lw=0, alpha=0.65)
xs = np.linspace(X.min(), X.max(), 50)
A.plot(10 ** xs, 10 ** (lr.intercept + lr.slope * xs), c=INK, lw=1.6)
A.scatter([10 ** X[top]], [10 ** Y[top]], s=70, c=MP, edgecolor=INK, lw=0.9, zorder=5)
A.annotate("$\\it{Pichia}$ at the most\nenriched site\nresidual %s SD" % num("%.2f", zz[top]),
           xy=(10 ** X[top], 10 ** Y[top]), xytext=(0.62, 0.42),
           textcoords="axes fraction", fontsize=7, c=INK,
           arrowprops=dict(arrowstyle="-", lw=0.7, color=INK))
A.set_xscale("log"); A.set_yscale("log")
A.set_xlabel("Reads assigned to the genus")
A.set_ylabel("Reads per distinct minimizer")
A.text(0.03, 0.96, "slope %.3f, R$^2$ = %.3f, n = %d" % (lr.slope, lr.rvalue ** 2, len(g)),
       transform=A.transAxes, va="top", fontsize=7.5, c=INK)
fig.tight_layout(pad=0.8)
save(fig, "FigS4_complexity")
print("=== FigS4 ===")
chk("complexity slope", lr.slope, 0.472, 5e-4)
chk("complexity R2", lr.rvalue ** 2, 0.560, 1e-3)
chk("observations", float(len(g)), 686.0)
chk("top observation is Pichia", g[top]["genus"], "Pichia")
chk("top observation residual, SD", zz[top], -1.03, 0.02)
chk("top observation reads per distinct minimizer",
    float(g[top]["reads_per_distinct_minimizer"]), 0.442, 5e-4)

# ------------------------------------------------------------------- FigS6
st = [r for r in T("TableS20_structure.tsv") if r["arm"] in ("plastic", "surface_control")]
by = {(r["site"], r["arm"]): r for r in st}
sites = sorted({r["site"] for r in st})
fig, ax = plt.subplots(1, 2, figsize=(7.4, 3.4))
A = ax[0]
for s in sites:
    for a, c in (("plastic", MP), ("surface_control", NP)):
        A.plot([0, 1], [float(by[(s, a)]["genera_observed"]),
                        float(by[(s, a)]["expected_genera_rarefied"])],
               "-o", c=c, ms=4, lw=1.0, alpha=0.85)
A.set_xticks([0, 1]); A.set_xticklabels(["observed", "at matched\nfungal depth"])
A.set_xlim(-0.25, 1.25)
A.set_ylabel("Fungal genera")
drare = np.array([float(by[(s, "plastic")]["expected_genera_rarefied"])
                  - float(by[(s, "surface_control")]["expected_genera_rarefied"])
                  for s in sites])
A.plot([], [], "-o", c=MP, ms=4, lw=1.0, label="Microplastic")
A.plot([], [], "-o", c=NP, ms=4, lw=1.0, label="Natural particle")
A.legend(fontsize=6.6, loc="lower left", handletextpad=0.5, labelspacing=0.3,
         borderpad=0.2)
obs = [float(by[(s2, a)]["genera_observed"]) for s2 in sites for a in ("plastic", "surface_control")]
rar = [float(by[(s2, a)]["expected_genera_rarefied"]) for s2 in sites for a in ("plastic", "surface_control")]
A.set_ylim(min(obs + rar) - 1.5, max(obs + rar) + 6.0)
A.text(0.03, 0.97, "mean rarefied difference %s genera" % num("%.1f", drare.mean()),
       transform=A.transAxes, va="top", fontsize=7, c=INK)
A.set_title("a", loc="left", fontweight="bold", fontsize=11)
B = ax[1]
dsimp = np.array([float(by[(s, "plastic")]["inverse_simpson"])
                  - float(by[(s, "surface_control")]["inverse_simpson"]) for s in sites])
for s in sites:
    B.plot([0, 1], [float(by[(s, "plastic")]["inverse_simpson"]),
                    float(by[(s, "surface_control")]["inverse_simpson"])],
           "-o", c=GREY, ms=4, lw=1.0)
B.set_xticks([0, 1]); B.set_xticklabels(["microplastic", "natural particle"])
B.set_xlim(-0.25, 1.25)
B.set_ylabel("Inverse Simpson diversity")
wsimp = stats.wilcoxon(dsimp)
B.text(0.5, 0.96, "mean difference %s\nsigned-rank P = %.4f"
       % (num("%.2f", dsimp.mean()), wsimp.pvalue), transform=B.transAxes,
       ha="center", va="top", fontsize=7, c=INK)
B.set_title("b", loc="left", fontweight="bold", fontsize=11)
fig.tight_layout(pad=1.0)
save(fig, "FigS8_structure")
print("=== FigS6 ===")
chk("mean rarefied richness difference", drare.mean(), -3.6253, 1e-3)
chk("annotation rounds to one decimal", float("%.1f" % drare.mean()), -3.6, 1e-9)
chk("mean inverse Simpson difference", dsimp.mean(), -0.9379, 1e-3)
chk("inverse Simpson signed-rank P", wsimp.pvalue, 0.6788, 1e-3)

# ------------------------------------------------------------------- FigS10
ph = T("TableS23_phylum_concordance.tsv")
ph = [r for r in ph if r["phylum_GTDB"] and r["bao_log2_MP_over_NP"]
      and r["sign_agree"] in ("yes", "no")]
# Marker area and the stratification use the number of genomes recovered per
# phylum by the published reconstruction, read from the source study's
# Supporting Information and written back into Table S23.
NG = {r["phylum"]: float(r["n_rMAGs"])
      for r in T("bao2026_phylum_abundance_FROM_SI.tsv")}
missing = [r["phylum_GTDB"] for r in ph if r["phylum_GTDB"] not in NG]
if missing:
    raise SystemExit("no genome count for %s" % ", ".join(missing))
for r in ph:
    r["n_rMAGs"] = "%d" % NG[r["phylum_GTDB"]]
with open(os.path.join(D, "TableS23_phylum_concordance.tsv"), "w", newline="") as fh:
    wtr = csv.DictWriter(fh, delimiter="\t", lineterminator="\n",
                         fieldnames=["phylum_GTDB", "n_rMAGs", "bao_log2_MP_over_NP",
                                     "ours_log2_MP_over_NP", "sign_agree"])
    wtr.writeheader()
    for r in ph:
        wtr.writerow({k: r[k] for k in wtr.fieldnames})
print("Table S23 rewritten with the n_rMAGs column")
xs = np.array([float(r["bao_log2_MP_over_NP"]) for r in ph])
ys = np.array([float(r["ours_log2_MP_over_NP"]) for r in ph])
ng = np.array([NG[r["phylum_GTDB"]] for r in ph])
fig, A = plt.subplots(figsize=(6.4, 5.4))
big = ng >= 50
A.scatter(xs[~big], ys[~big], s=28 + 1.4 * ng[~big], c=GREY, edgecolor=INK, lw=0.7, zorder=3)
A.scatter(xs[big], ys[big], s=28 + 1.4 * ng[big], c=NP, edgecolor=INK, lw=0.7, zorder=3)
lo = min(xs.min(), ys.min()) - 0.4
hi = max(xs.max(), ys.max()) + 0.4
A.plot([lo, hi], [lo, hi], ls=":", c=INK, lw=1.0, zorder=1)
A.axhline(0, c="#c8cbd0", lw=0.8, zorder=0)
A.axvline(0, c="#c8cbd0", lw=0.8, zorder=0)
# Label placement. Candidate offsets are tried in display-normalised units
# and scored against every marker and every label already placed, so no label
# sits on a point or on another label.
xr = hi - lo
yr = (ys.max() + 0.6) - (ys.min() - 0.6)
cand = []
for ang in range(0, 360, 20):
    for rad in (0.030, 0.052, 0.080):
        cand.append((rad * math.cos(math.radians(ang)),
                     rad * math.sin(math.radians(ang))))
pts = [(xs[k] / xr, ys[k] / yr) for k in range(len(ph))]
placed = []
order = np.argsort(-ng)
for k in order:
    best, bestd = None, -1.0
    for dx, dy in cand:
        px, py = xs[k] / xr + dx, ys[k] / yr + dy
        dmark = min(math.hypot(px - qx, py - qy) for qx, qy in pts)
        dlab = min([math.hypot(px - qx, py - qy) for qx, qy in placed] or [9e9])
        d = min(dmark * 1.0, dlab * 0.75)
        if d > bestd:
            best, bestd, bdx, bdy = (px, py), d, dx, dy
    placed.append(best)
    A.annotate(ph[k]["phylum_GTDB"], (xs[k], ys[k]),
               xytext=(best[0] * xr, best[1] * yr),
               ha="left" if bdx > 0.012 else "right" if bdx < -0.012 else "center",
               va="bottom" if bdy > 0.012 else "top" if bdy < -0.012 else "center",
               fontsize=6.8, c=INK,
               arrowprops=dict(arrowstyle="-", lw=0.5, color="#9aa0a8",
                               shrinkA=0.5, shrinkB=3.0))
A.set_ylim(ys.min() - 0.6, ys.max() + 0.6)
A.set_xlabel("Genome-resolved log$_2$(MP/NP), source study")
A.set_ylabel("This study, read-based log$_2$(MP/NP)")
A.text(0.03, 0.97, "marker area scales with the number of recovered\n"
       "genomes. green marks at least 50", transform=A.transAxes, va="top",
       fontsize=7, c=INK)
A.text(0.03, 0.04, "%d of %d phyla agree in sign"
       % (sum(1 for r in ph if r["sign_agree"] == "yes"), len(ph)),
       transform=A.transAxes, fontsize=7, c=INK)
A.set_xlim(lo, hi)
fig.tight_layout(pad=0.8)
save(fig, "FigS10_phylum")
print("=== FigS10 ===")
chk("phyla plotted", float(len(ph)), 15.0)
chk("phyla agreeing in sign", float(sum(1 for r in ph if r["sign_agree"] == "yes")), 10.0)
chk("phyla with at least 50 genomes", float(big.sum()), 4.0)
chk("all with at least 50 genomes agree",
    float(sum(1 for r, b in zip(ph, big) if b and r["sign_agree"] == "yes")), 4.0)

# ------------------------------------------------------------------- FigS11
nz = T("TableS4_second_system.tsv")
nz = [r for r in nz if r.get("study", "").startswith("PRJ")]
blocks = T("TableS1_blocks.tsv")
hk = {s: {r["site"]: float(r["d"]) for r in blocks if r["season"] == s}
      for s in ("Summer", "Winter")}
com = sorted(set(hk["Summer"]) & set(hk["Winter"]))
hkd = np.array([hk["Summer"][s] - hk["Winter"][s] for s in com])
hkse = stats.sem(hkd); hkdf = len(hkd) - 1
hkm = (stats.t.ppf(0.975, hkdf) + stats.t.ppf(0.80, hkdf)) * hkse
SYS = [("Hong Kong rivers", hkd.mean(), stats.t.ppf(0.975, hkdf) * hkse, hkm),
       ("Auckland marina", 0.894, 1.931, 2.41),
       ("NZ polishing pond", 0.254, 2.708, 3.38)]
fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.1))
A = ax[0]
cols = {"Summer": MP, "Autumn": AMBER, "Winter": RW, "Spring": NP}
for i, r in enumerate(nz):
    A.barh(i, float(r["difference_log2"]), 0.6, color=cols[r["season_SH"]])
A.axvline(0, c=INK, lw=0.8)
A.set_yticks(range(len(nz)))
A.set_yticklabels(["%s  %s (%s)" % (r["system"].split()[0], r["date"], r["season_SH"])
                   for r in nz], fontsize=6.3)
A.invert_yaxis()
A.set_xlabel("Plastic minus glass control (log$_2$)")
A.set_title("a", loc="left", fontweight="bold", fontsize=11)
B = ax[1]
for i, (nm, est, ci, md) in enumerate(SYS):
    B.barh(i, md, 0.5, color=GREY, alpha=0.6,
           label="Detectable at 80% power" if i == 0 else None)
    B.errorbar([est], [i], xerr=[ci], fmt="o", ms=6, c=MP, lw=1.3, capsize=4,
               label="Observed, with 95% CI" if i == 0 else None)
B.axvline(0, c=INK, lw=0.8, ls="--")
B.set_yticks(range(len(SYS)))
B.set_yticklabels([s[0] for s in SYS], fontsize=7)
B.invert_yaxis()
B.set_xlabel("Summer minus winter contrast (log$_2$)")
B.set_title("b", loc="left", fontweight="bold", fontsize=11)
B.legend(fontsize=6.0, loc="upper center", bbox_to_anchor=(0.5, -0.20), ncol=2,
         handletextpad=0.5, labelspacing=0.3, borderpad=0.2, columnspacing=1.4)
fig.tight_layout(pad=1.0)
save(fig, "FigS11_second_system")
print("=== FigS11 ===")
chk("per-date contrasts", float(len(nz)), 6.0)
chk("Hong Kong summer minus winter", hkd.mean(), 4.2710, 0.02)
chk("Auckland detectable bound, fold", 2 ** 2.41, 5.31, 0.02)
chk("pond detectable bound, fold", 2 ** 3.38, 10.41, 0.02)

# ------------------------------------------------------------------- FigS12
th = T("TableS24_thermotolerance.tsv")
th = sorted(th, key=lambda r: -float(r["plastic_slope"]))
fig, ax = plt.subplots(1, 2, figsize=(7.4, 3.3))
for k, (axx, pre, title) in enumerate(((ax[0], "plastic", "Abundance on plastic"),
                                       (ax[1], "difference",
                                        "Plastic minus natural particle"))):
    v = [float(r[pre + "_slope"]) for r in th]
    p = [float(r[pre + "_P_holm"]) for r in th]
    axx.barh(range(len(th)), v, 0.6,
             color=[MP if q < 0.05 else GREY for q in p])
    for i, (vi, pi) in enumerate(zip(v, p)):
        axx.text(vi + max(v) * 0.03, i,
                 "P = %.4f" % pi if pi >= 0.0001 else "P < 0.0001",
                 va="center", fontsize=6.3, c=INK)
    axx.set_yticks(range(len(th)))
    axx.set_yticklabels(["$\\it{%s}$" % r["genus"] for r in th], fontsize=7.5)
    axx.invert_yaxis()
    axx.set_xlim(0, max(v) * 1.42)
    axx.set_xlabel("Temperature slope (log$_2$ units °C$^{-1}$)")
    axx.set_title("%s   %s" % ("ab"[k], title), loc="left", fontweight="bold",
                  fontsize=9)
ax[0].text(0.97, 0.05, "filled, Holm-adjusted P < 0.05", transform=ax[0].transAxes,
           ha="right", fontsize=6.3, c=INK)
fig.tight_layout(pad=1.0)
save(fig, "FigS12_thermotolerance")
print("=== FigS12 ===")
chk("row order starts with Pichia", th[0]["genus"], "Pichia")
chk("Pichia plastic slope", float(th[0]["plastic_slope"]), 0.8062, 1e-4)
chk("Pichia plastic Holm P", float(th[0]["plastic_P_holm"]), 0.0001, 1e-6)
chk("Lodderomyces plastic slope", float(th[1]["plastic_slope"]), 0.5359, 1e-4)
chk("Candida plastic slope", float(th[2]["plastic_slope"]), 0.3863, 1e-4)
chk("Candida plastic Holm P", float(th[2]["plastic_P_holm"]), 0.0132, 1e-6)
chk("plastic arm, genera surviving Holm",
    float(sum(1 for r in th if float(r["plastic_P_holm"]) < 0.05)), 3.0)
chk("difference, genera surviving Holm",
    float(sum(1 for r in th if float(r["difference_P_holm"]) < 0.05)), 2.0)
dif = {r["genus"]: float(r["difference_slope"]) for r in th}
chk("Yarrowia above Candida on the difference",
    float(dif["Yarrowia"] > dif["Candida"]), 1.0)


# ====================================================== S3, S5, S6
# ---------------------------------------------------------------- FigS5
yc = T("TableS17_yeast_carryover.tsv")
ARMS = [("plastic", "Microplastic", MP), ("surface_control", "Natural particle", NP),
        ("bulk_control", "River water", RW)]
drop = [r for r in yc if r["arm"] not in dict((a, 1) for a, _, _ in ARMS)]
print("FigS5: %d summary row(s) excluded from the arm series: %s"
      % (len(drop), [r["arm"] for r in drop]))
fig, ax = plt.subplots(1, 2, figsize=(6.6, 2.9))
for k, (col, lab) in enumerate((("median_yeast_RPM", "Median"), ("max_yeast_RPM", "Maximum"))):
    A = ax[k]
    w = 0.8 / len(ARMS)
    for j, (a, nm, c) in enumerate(ARMS):
        v = []
        for s in SE:
            g = [r for r in yc if r["season"] == s and r["arm"] == a]
            v.append(max(float(g[0][col]), 0.01) if g and g[0][col] else 0.01)
        A.bar(np.arange(4) + (j - len(ARMS) / 2.0 + 0.5) * w, v, w, color=c, label=nm)
    A.set_yscale("log"); A.set_xticks(range(4)); A.set_xticklabels(SE, fontsize=7)
    A.set_ylabel("%s yeast reads per million" % lab)
    A.set_title("ab"[k], loc="left", fontweight="bold", fontsize=11)
    if k == 0:
        A.legend(fontsize=6.2, loc="upper left", handletextpad=0.5, labelspacing=0.25)
fig.tight_layout(pad=1.0)
for e in ("png", "pdf"):
    fig.savefig(os.path.join(OUT, "FigS5_yeast_carryover." + e),
                **({"dpi": 400} if e == "png" else {}), bbox_inches="tight")
print("   wrote FigS5 (a medians, b maxima; 3 arms)")

# ---------------------------------------------------------------- FigS7
ic = T("TableS18_index_coverage.tsv")
fu = [r for r in ic if r["kingdom_group"] == "fungal"]
ot = [r for r in ic if r["kingdom_group"] != "fungal"]
print("FigS6: %d genera total = %d fungal + %d non-fungal" % (len(ic), len(fu), len(ot)))
fig, A = plt.subplots(figsize=(5.0, 3.2))
for rows, c, s, lab in ((ot, GREY, 12, "Non-fungal genera (n = %d)" % len(ot)),
                        (fu, MP, 16, "Fungal genera (n = %d)" % len(fu))):
    A.scatter([int(r["full_minimizers"]) for r in rows],
              [float(r["capped_over_full"]) for r in rows],
              s=s, c=c, lw=0, alpha=0.85, label=lab)
rr = [float(r["capped_over_full"]) for r in ic]
rf = [float(r["capped_over_full"]) for r in fu]
A.set_xscale("log"); A.set_ylim(0.125, 0.145)
A.set_xlabel("Minimizers in the full index")
A.set_ylabel("Capped / full retention")
A.legend(fontsize=6.4, loc="lower right")
A.text(0.03, 0.06, "all genera %.4f–%.4f\nfungal %.4f–%.4f"
       % (min(rr), max(rr), min(rf), max(rf)), transform=A.transAxes, fontsize=6.4, c=INK)
fig.tight_layout(pad=1.0)
for e in ("png", "pdf"):
    fig.savefig(os.path.join(OUT, "FigS6_index_equivalence." + e),
                **({"dpi": 400} if e == "png" else {}), bbox_inches="tight")
print("   wrote FigS6")

# ====================================================== S3, S5, S6
# ---------------------------------------------------------------- FigS5
yc = T("TableS17_yeast_carryover.tsv")
ARMS = [("plastic", "Microplastic", MP), ("surface_control", "Natural particle", NP),
        ("bulk_control", "River water", RW)]
drop = [r for r in yc if r["arm"] not in dict((a, 1) for a, _, _ in ARMS)]
print("FigS5: %d summary row(s) excluded from the arm series: %s"
      % (len(drop), [r["arm"] for r in drop]))
fig, ax = plt.subplots(1, 2, figsize=(6.6, 2.9))
for k, (col, lab) in enumerate((("median_yeast_RPM", "Median"), ("max_yeast_RPM", "Maximum"))):
    A = ax[k]
    w = 0.8 / len(ARMS)
    for j, (a, nm, c) in enumerate(ARMS):
        v = []
        for s in SE:
            g = [r for r in yc if r["season"] == s and r["arm"] == a]
            v.append(max(float(g[0][col]), 0.01) if g and g[0][col] else 0.01)
        A.bar(np.arange(4) + (j - len(ARMS) / 2.0 + 0.5) * w, v, w, color=c, label=nm)
    A.set_yscale("log"); A.set_xticks(range(4)); A.set_xticklabels(SE, fontsize=7)
    A.set_ylabel("%s yeast reads per million" % lab)
    A.set_title("ab"[k], loc="left", fontweight="bold", fontsize=11)
    if k == 0:
        A.legend(fontsize=6.2, loc="upper left", handletextpad=0.5, labelspacing=0.25)
fig.tight_layout(pad=1.0)
for e in ("png", "pdf"):
    fig.savefig(os.path.join(OUT, "FigS5_yeast_carryover." + e),
                **({"dpi": 400} if e == "png" else {}), bbox_inches="tight")
print("   wrote FigS5 (a medians, b maxima; 3 arms)")

# ---------------------------------------------------------------- FigS7
ic = T("TableS18_index_coverage.tsv")
fu = [r for r in ic if r["kingdom_group"] == "fungal"]
ot = [r for r in ic if r["kingdom_group"] != "fungal"]
print("FigS6: %d genera total = %d fungal + %d non-fungal" % (len(ic), len(fu), len(ot)))
fig, A = plt.subplots(figsize=(5.0, 3.2))
for rows, c, s, lab in ((ot, GREY, 12, "Non-fungal genera (n = %d)" % len(ot)),
                        (fu, MP, 16, "Fungal genera (n = %d)" % len(fu))):
    A.scatter([int(r["full_minimizers"]) for r in rows],
              [float(r["capped_over_full"]) for r in rows],
              s=s, c=c, lw=0, alpha=0.85, label=lab)
rr = [float(r["capped_over_full"]) for r in ic]
rf = [float(r["capped_over_full"]) for r in fu]
A.set_xscale("log"); A.set_ylim(0.125, 0.145)
A.set_xlabel("Minimizers in the full index")
A.set_ylabel("Capped / full retention")
A.legend(fontsize=6.4, loc="lower right")
A.text(0.03, 0.06, "all genera %.4f–%.4f\nfungal %.4f–%.4f"
       % (min(rr), max(rr), min(rf), max(rf)), transform=A.transAxes, fontsize=6.4, c=INK)
fig.tight_layout(pad=1.0)
for e in ("png", "pdf"):
    fig.savefig(os.path.join(OUT, "FigS6_index_equivalence." + e),
                **({"dpi": 400} if e == "png" else {}), bbox_inches="tight")
print("   wrote FigS6")

# ---------------------------------------------------------------- FigS3
# Recomputed from TableS1_blocks.tsv rather than read from the slope-sensitivity
# table, so each row carries its own slope, P value and count of positive sites.
from scipy import stats as _st
B = T("TableS1_blocks.tsv")
_d = np.array([float(r["d"]) for r in B]); _t = np.array([float(r["temp"]) for r in B])
_se = np.array([r["season"] for r in B]); _si = np.array([r["site"] for r in B])
SITES = sorted(set(_si))
def _slopes(mask):
    out = []
    for s in SITES:
        m = (_si == s) & mask
        if m.sum() >= 3 and _t[m].std() > 0:
            out.append(_st.linregress(_t[m], _d[m]).slope)
    return np.array(out)
VAR = [("All four seasons", np.ones(len(_d), bool)),
       ("Summer dropped",   _se != "Summer"),
       ("Winter dropped",   _se != "Winter"),
       ("Warm seasons only", _se != "Winter")]
VAR = VAR[:3]
fig, A = plt.subplots(figsize=(4.8, 2.4))
for i, (lab, m) in enumerate(VAR):
    a = _slopes(m); tt = _st.ttest_1samp(a, 0.0)
    ci = _st.t.interval(0.95, len(a) - 1, loc=a.mean(), scale=_st.sem(a))
    A.errorbar([a.mean()], [i], xerr=[[a.mean() - ci[0]], [ci[1] - a.mean()]],
               fmt="o", ms=6, c=MP, lw=1.3, capsize=4)
    A.text(ci[1] + 0.05, i, "%d/%d sites positive, P = %.4f"
           % (int((a > 0).sum()), len(a), tt.pvalue), va="center", fontsize=6.2, c=INK)
A.axvline(0, c=INK, lw=0.8, ls="--")
A.set_yticks(range(len(VAR))); A.set_yticklabels([v[0] for v in VAR], fontsize=7)
A.invert_yaxis()
A.set_xlabel("Within-site temperature slope (log$_2$ per °C)")
A.set_xlim(-0.1, 2.2)
fig.tight_layout(pad=1.0)
for e in ("png", "pdf"):
    fig.savefig(os.path.join(OUT, "FigS3_slope_sensitivity." + e),
                **({"dpi": 400} if e == "png" else {}), bbox_inches="tight")
print("FigS3: recomputed from TableS1; each row now carries its own P")
for lab, m in VAR:
    a = _slopes(m)
    print("   %-20s slope %+.4f  P = %.4f  %d/%d positive"
          % (lab, a.mean(), _st.ttest_1samp(a, 0.0).pvalue, int((a > 0).sum()), len(a)))

# ============================================================ S7, S9
blocks = T("TableS1_blocks.tsv")
rate = {r["run"]: float(r["class_rate"]) for r in T("TableS2_detection_envelope.tsv")}
bact = {(r["season"], r["site"]): r for r in T("TableS22_bacterial_blocks.tsv")}

temp = np.array([float(r["temp"]) for r in blocks])
dfun = np.array([float(r["d"]) for r in blocks])
# the classification-rate contrast, in percentage points, within each block
drate = np.array([rate[r["run_MP"]] - rate[r["run_NP"]] for r in blocks])
# and on the log2 scale, which is what the bacterial endpoint is compared with
dlograte = np.array([math.log2(rate[r["run_MP"]] / rate[r["run_NP"]]) for r in blocks])
dbact = np.array([float(bact[(r["season"], r["site"])]["d_log2"]) for r in blocks])
season = np.array([r["season"] for r in blocks])
site = np.array([r["site"] for r in blocks])
sites = sorted(set(site))
print("blocks %d   sites %d" % (len(blocks), len(sites)))

# ------------------------------------------------------------------ Figure S7
r_fun = stats.pearsonr(drate, dfun)
r_tmp = stats.pearsonr(drate, temp)
sl_tmp = stats.linregress(temp, drate)
within = []
for s in sites:
    m = site == s
    if np.std(temp[m]) > 0:
        within.append(stats.linregress(temp[m], drate[m]).slope)
within = np.array(within)
t_within = stats.ttest_1samp(within, 0.0)

fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.1))
A = ax[0]
A.scatter(drate, dfun, s=18, c=GREY, lw=0, alpha=0.85)
A.axhline(0, c=INK, lw=0.6, ls=":")
A.axvline(0, c=INK, lw=0.6, ls=":")
A.set_xlabel("Difference in classification rate (percentage points)")
A.set_ylabel("Difference in log$_2$ fungal fraction")
A.text(0.03, 0.97, "r = %s, R$^2$ = %.3f\nP = %.3f, n = %d"
       % (num("%+.3f", r_fun[0]), r_fun[0] ** 2, r_fun[1], len(drate)),
       transform=A.transAxes, va="top", fontsize=7, c=INK)
A.set_title("a", loc="left", fontweight="bold", fontsize=11)
B = ax[1]
B.scatter(temp, drate, s=18, c=GREY, lw=0, alpha=0.85)
B.axhline(0, c=INK, lw=0.9, ls="--")
B.set_xlabel("Water temperature (%sC)" % DEG)
B.set_ylabel("Difference in classification rate\n(percentage points)")
B.text(0.03, 0.97, "r = %s\nP = %.3f" % (num("%+.3f", r_tmp[0]), r_tmp[1]),
       transform=B.transAxes, va="top", fontsize=7, c=INK)
B.set_title("b", loc="left", fontweight="bold", fontsize=11)
fig.tight_layout(pad=1.0)
save(fig, "FigS7_classrate_confound")

print("=== FigS7, against the values reported in the caption and Table S19 ===")
chk("classification-rate contrast predicts the fungal contrast, r", r_fun[0], 0.1818, 5e-4)
chk("  R2", r_fun[0] ** 2, 0.0331, 5e-4)
chk("  P", r_fun[1], 0.1644, 5e-4)
chk("rate contrast against temperature, pooled r", r_tmp[0], -0.0663, 5e-4)
chk("  slope, points per degree", sl_tmp.slope, -0.1437, 5e-4)
chk("  P", r_tmp[1], 0.6148, 5e-4)
chk("rate contrast against temperature, within-site mean slope",
    float(within.mean()), -0.2670, 5e-4)
chk("  P", float(t_within.pvalue), 0.2771, 5e-4)

# ------------------------------------------------------------------ Figure S9
fslope = {r["site"]: float(r["slope_log2_per_C"]) for r in T("TableS14_site_slopes.tsv")}
fsl = np.array([fslope[s] for s in sites])
bsl = []
for s in sites:
    m = site == s
    bsl.append(stats.linregress(temp[m], dbact[m]).slope)
bsl = np.array(bsl)
tf, tb = stats.ttest_1samp(fsl, 0.0), stats.ttest_1samp(bsl, 0.0)

# Table S19 is a stacked file. csv.DictReader would read the second block with
# the first block's field names, so the per-season block is parsed explicitly.
rows, head = [], None
for line in open(os.path.join(D, "TableS19_classrate_confound.tsv")):
    q = line.rstrip("\n").split("\t")
    if q[0] == "season":
        head = q
        continue
    if head and len(q) == len(head) and q[0] in SE:
        rows.append(dict(zip(head, q)))
assert len(rows) == 4, rows
fmean = {r["season"]: float(r["mean_dFungal_log2"]) for r in rows}
cmean = {r["season"]: float(r["mean_dClassRate_points"]) for r in rows}
bmean = {s: float(np.mean(dbact[season == s])) for s in SE}
tmean = {s: float(np.mean(temp[season == s])) for s in SE}

r_coll = stats.pearsonr(dlograte, dbact)

fig, ax = plt.subplots(1, 3, figsize=(7.8, 2.9))
A = ax[0]
rng = np.random.default_rng(0)
for k, (v, c, lab, tt) in enumerate(((fsl, FUNGI, "Fungi", tf), (bsl, BACT, "Bacteria", tb))):
    A.scatter(k + (rng.random(len(v)) - 0.5) * 0.3, v, s=16, c=c, lw=0, alpha=0.9)
    ci = stats.t.ppf(0.975, len(v) - 1) * stats.sem(v)
    A.errorbar([k], [v.mean()], yerr=[ci], fmt="o", ms=7, c=INK, lw=1.2, capsize=4, zorder=4)
    # annotation text is black, the points carry the colour
    A.text(k, max(v) + 0.10, "%s\nP = %.3f" % (num("%+.3f", v.mean()), tt.pvalue),
           ha="center", va="bottom", fontsize=6.8, c=INK)
A.axhline(0, c=GREY, lw=0.9, ls="--")
A.set_xticks([0, 1]); A.set_xticklabels(["Fungi", "Bacteria"])
A.set_xlim(-0.6, 1.6)
A.set_ylim(min(bsl.min(), fsl.min()) - 0.12, max(bsl.max(), fsl.max()) + 0.62)
A.set_ylabel("Temperature slope of MP$-$NP\n(log$_2$ units %sC$^{-1}$)" % DEG)
A.set_title("a", loc="left", fontweight="bold", fontsize=11)
B = ax[1]
for mean, c, lab in ((fmean, FUNGI, "Fungi"), (bmean, BACT, "Bacteria")):
    o = sorted(SE, key=lambda s: tmean[s])
    B.plot([tmean[s] for s in o], [mean[s] for s in o], "-o", c=c, ms=6, lw=1.4, label=lab)
B.axhline(0, c=GREY, lw=0.9, ls="--")
B.set_xlabel("Mean water temperature (%sC)" % DEG)
B.set_ylabel("Mean MP$-$NP difference (log$_2$)")
B.legend(fontsize=6.6, loc="upper left")
B.set_title("b", loc="left", fontweight="bold", fontsize=11)
C = ax[2]
C.scatter(dlograte, dbact, s=16, c=BACT, lw=0, alpha=0.85)
lo = min(dlograte.min(), dbact.min()) - 0.3
hi = max(dlograte.max(), dbact.max()) + 0.3
C.plot([lo, hi], [lo, hi], ls=":", c=INK, lw=0.9)
C.set_xlabel("Difference in log$_2$ classification rate")
C.set_ylabel("Difference in log$_2$ bacterial fraction")
C.text(0.03, 0.97, "r = %s\nR$^2$ = %.4f" % (num("%+.4f", r_coll[0]), r_coll[0] ** 2),
       transform=C.transAxes, va="top", fontsize=7, c=INK)
C.set_title("c", loc="left", fontweight="bold", fontsize=11)
fig.tight_layout(pad=1.0)
save(fig, "FigS9_bacterial")

print("\n=== FigS9, against the values reported in the caption and Table S21 ===")
chk("fungal within-site temperature slope, mean", float(fsl.mean()), 0.4133, 5e-4)
chk("  P", float(tf.pvalue), 0.0079, 5e-4)
chk("bacterial within-site temperature slope, mean", float(bsl.mean()), -0.0500, 5e-3)
chk("bacterial endpoint against classification rate, r", r_coll[0], 0.9986, 5e-4)
chk("  R2", r_coll[0] ** 2, 0.9971, 5e-4)
for s, want in zip(SE, (-1.2215, -0.2577, 3.0495, 0.9473)):
    chk("fungal mean contrast, %-7s" % s, fmean[s], want, 5e-4)
for s, want in zip(SE, (3.3216, -1.3854, 2.0695, -0.7110)):
    chk("classification-rate mean contrast, %-7s" % s, cmean[s], want, 5e-4)
chk("per-block rates reproduce the per-season fungal means",
    float(np.mean(dfun[season == "Summer"])), 3.0495, 5e-4)
chk("per-block rates reproduce the per-season rate means",
    float(np.mean(drate[season == "Winter"])), 3.3216, 5e-4)
chk("winter bacterial mean contrast", bmean["Winter"], 0.5200, 5e-3)
chk("summer fungal mean contrast", fmean["Summer"], 3.0495, 5e-4)
chk("winter fungal mean contrast", fmean["Winter"], -1.2215, 5e-4)
chk("bacterial winter is the highest of the four seasons",
    float(bmean["Winter"] == max(bmean.values())), 1.0)
chk("fungal winter is the lowest of the four seasons",
    float(fmean["Winter"] == min(fmean.values())), 1.0)

print("")
print("\n%d checks, %d mismatches" % (nchk[0], len(fails)))
if fails:
    raise SystemExit("mismatched: " + ", ".join(fails))
