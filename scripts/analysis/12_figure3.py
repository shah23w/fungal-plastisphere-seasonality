"""Figure 3, fungal yield across the archived runs and false fungal assignment in
the fungus-free mock community.

Panel d fits the 2 x 2 factorial on log10 false positives per million with
library as a blocking factor, which is the model the Methods describe. Two cells
return no false reads and are undefined on the log scale, leaving 34 of 36. The
script also refits all 36 cells with a half-read continuity correction to show
that dropping them is conservative.

    python3 12_figure3.py [table_dir] [out_dir]

Exit status is 0 only if the fitted terms match the reported values.
"""
import csv, math, os, sys
import numpy as np
from scipy import stats
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
D = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "..", "tables")
OUT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, "..", "..", "figures")


def T(n):
    p = os.path.join(D, n)
    if not os.path.exists(p):
        raise SystemExit("missing %s" % p)
    return list(csv.DictReader(open(p), delimiter="\t"))


FLOOR = 0.16
env = T("TableS2_detection_envelope.tsv")
mk = T("TableS7_mock2x2.tsv")
print("envelope runs: %d   mock rows: %d" % (len(env), len(mk)))

INK, GREY, MP = "#22252a", "#b9bcc2", "#c1443a"
ARM = {"plastic": "#c1443a", "natural": "#2f7d6f", "water": "#4a6fa5"}
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.linewidth": 0.7,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "legend.frameon": False})
fig, ax = plt.subplots(2, 2, figsize=(7.6, 6.6))
fig.subplots_adjust(wspace=0.30, hspace=0.42, left=0.085,
                    right=0.985, top=0.955, bottom=0.085)

# ------------------------------------------------------------------ panel a
A = ax[0][0]
co = [r for r in env if not r["hk_arm"]]
hk = [r for r in env if r["hk_arm"]]
A.scatter([int(r["total_pairs"]) for r in co], [int(r["fungal_reads"]) + 0.5 for r in co],
          s=5, c=GREY, lw=0, alpha=0.75,
          label="Other studies (n=%d)" % len(co))
for a in ("water", "natural", "plastic"):
    g = [r for r in hk if r["hk_arm"] == a]
    A.scatter([int(r["total_pairs"]) for r in g], [int(r["fungal_reads"]) + 0.5 for r in g],
              s=11, c=ARM[a], lw=0, alpha=0.9,
              label="%s (n=%d)" % (a.capitalize(), len(g)))
xs = [10 ** 5.5, 10 ** 9.0]
A.plot(xs, [FLOOR * x / 1e6 for x in xs], ls="--", c=INK, lw=0.9)
for thr, lb in ((1000, "T1-T2"), (10000, "T1-T3")):
    A.axhline(thr, c="#8a8f98", lw=0.6, ls=":")
    A.text(10 ** 8.6, thr * 1.25, lb, fontsize=6.0, c="#6b7078", ha="right")
A.set_xscale("log"); A.set_yscale("log")
A.set_xlabel("Total read pairs"); A.set_ylabel("Fungal reads")
A.set_title("a", loc="left", fontweight="bold", fontsize=11)
A.legend(fontsize=5.8, loc="lower right", handletextpad=0.4, frameon=False,
         labelspacing=0.25, borderpad=0.2)

# ------------------------------------------------------------------ panel b
B = ax[0][1]
tiers = ["T1-T3", "T1-T2", "STOP"]
nm = {"T1-T3": "At least 10,000\n(full analysis)", "T1-T2": "1,000-9,999\n(limited)",
      "STOP": "Under 1,000\n(insufficient)"}
n = [sum(1 for r in env if r["gate1a"] == t) for t in tiers]
B.bar(range(3), n, color=["#2f7d6f", "#d9a43b", "#c1443a"], width=0.62)
for i, v in enumerate(n):
    B.text(i, v + 12, "%d\n(%.1f%%)" % (v, 100.0 * v / len(env)), ha="center",
           fontsize=6.8, c=INK)
B.set_xticks(range(3)); B.set_xticklabels([nm[t] for t in tiers], fontsize=6.6)
B.set_ylabel("Runs"); B.set_ylim(0, max(n) * 1.22); B.set_xlabel("Fungal reads per run")
B.set_title("b", loc="left", fontweight="bold", fontsize=11)

# ------------------------------------------------------------------ panel c
C = ax[1][0]
CELLS = [1, 2, 3, 4, 5]
LAB = {1: "Fungal-only\nconf. 0", 2: "Fungal-only\nconf. 0.1", 3: "PlusPF\nconf. 0",
       4: "PlusPF\nconf. 0.1", 5: "PlusPF\nconf. 0.1\nMHG 3"}
SC = {"A": "#4a6fa5", "B": "#2f7d6f", "C": "#c1443a"}
SN = {"A": "Bacteria", "B": "Out-of-database eukaryotes", "C": "Human"}
w = 0.26
for j, s in enumerate("ABC"):
    vals = []
    for c in CELLS:
        v = [float(r["false_per_million"]) for r in mk
             if int(r["cell"]) == c and r["stratum"] == s]
        vals.append(max(np.mean(v), 1e-3) if v else 1e-3)
    C.bar(np.arange(len(CELLS)) + (j - 1) * w, vals, w, color=SC[s], label=SN[s])
C.axhline(FLOOR, c=INK, lw=0.9, ls="--")
C.set_yscale("log"); C.set_xticks(range(len(CELLS)))
C.set_xticklabels([LAB[c] for c in CELLS], fontsize=6.0)
C.set_ylabel("False fungal reads per million")
C.set_title("c", loc="left", fontweight="bold", fontsize=11)
C.legend(fontsize=5.8, loc="upper right", frameon=False, labelspacing=0.25,
         handletextpad=0.5, borderpad=0.2)

# ------------------------------------------------------------------ panel d
# Library-blocked factorial on log10 false positives per million, cells 1-4.
# Two cells returned zero false fungal reads and are undefined on the log
# scale, so 34 of the 36 cells enter, exactly as the Methods state.
P = ax[1][1]
sub = [r for r in mk if int(r["cell"]) in (1, 2, 3, 4) and float(r["false_per_million"]) > 0]
y = np.array([math.log10(float(r["false_per_million"])) for r in sub])
sc = np.array([0.0 if r["index"] == "fungal-only" else 1.0 for r in sub])
cf = np.array([0.0 if float(r["confidence"]) == 0.0 else 1.0 for r in sub])
runs = sorted({r["run"] for r in sub})
blocks = [np.array([1.0 if r["run"] == rn else 0.0 for r in sub]) for rn in runs[1:]]
X = np.column_stack([np.ones(len(y)), sc, cf, sc * cf] + blocks)
beta, *_ = np.linalg.lstsq(X, y, rcond=None)
res = y - X @ beta
df = len(y) - np.linalg.matrix_rank(X)
se = np.sqrt(np.diag(np.linalg.pinv(X.T @ X)) * (res @ res) / df)
names = ["Capped PlusPF\nindex", "Confidence\n0.1", "Interaction"]
vals = [-beta[1], -beta[2], -beta[3]]
errs = [stats.t.ppf(0.975, df) * se[i] for i in (1, 2, 3)]
ps = [2 * stats.t.sf(abs(beta[i] / se[i]), df) for i in (1, 2, 3)]
cols = [MP if p < 0.05 else GREY for p in ps]
P.barh(range(3), vals, xerr=errs, color=cols, height=0.55,
       error_kw=dict(lw=1.0, capsize=4, ecolor=INK))
for i in range(3):
    P.text(vals[i] + errs[i] + 0.12, i, "%.1f×" % (10 ** vals[i]),
           va="center", fontsize=6.4, c=INK)
P.axvline(0, c=INK, lw=0.7); P.set_yticks(range(3))
P.set_yticklabels(names, fontsize=6.6); P.invert_yaxis()
P.set_xlabel("Reduction in false positives (log$_{10}$)")
P.set_title("d", loc="left", fontweight="bold", fontsize=11)

fig.tight_layout(pad=1.1)
os.makedirs(OUT, exist_ok=True)
for ext, kw in (("png", {"dpi": 400}), ("pdf", {}),
                ("tiff", {"dpi": 400, "pil_kwargs": {"compression": "tiff_lzw"}})):
    fig.savefig(os.path.join(OUT, "Figure3.%s" % ext), bbox_inches="tight", **kw)
print("wrote Figure3 to %s" % OUT)

# ------------------------------------------------------------------ checks
fails = []
nchk = [0]


def chk(name, got, want, tol):
    nchk[0] += 1
    ok = abs(got - want) <= tol
    print("   %-44s %14.6g  want %-14.6g %s" % (name, got, want, "OK" if ok else "MISMATCH"))
    if not ok:
        fails.append(name)


print("\n=== VERIFICATION, panel d against the reported values ===")
print("   cells entering the model: %d of 36   residual df: %d" % (len(sub), df))
for i, (lab, wf, wp) in enumerate(
        (("reference scope", 664.09, 3.726e-10),
         ("confidence 0.1", 4.4267, 0.02319),
         ("interaction", 11.293, 0.01347)), start=1):
    chk("%s, fold" % lab, 10 ** -beta[i], wf, max(0.02, wf * 0.002))
    chk("%s, P" % lab, ps[i - 1], wp, max(1e-11, wp * 0.02))
chk("all three bars significant", float(sum(1 for p in ps if p < 0.05)), 3.0, 0.0)

print("\n=== VERIFICATION, panels a and b ===")
chk("runs from other studies", float(len(co)), 679.0, 0.0)
chk("Hong Kong libraries", float(len(hk)), 180.0, 0.0)
chk("tier at least 10,000", float(n[0]), 21.0, 0.0)
chk("tier 1,000-9,999", float(n[1]), 132.0, 0.0)
chk("tier under 1,000", float(n[2]), 706.0, 0.0)
chk("percent under 1,000", 100.0 * n[2] / len(env), 82.2, 0.05)

print("\n=== CONTROL, the zero-preserving specification ===")
# Every one of the 36 cells, with a half-read continuity correction, to show
# that dropping the two structural zeros was conservative for the interaction.
allc = [r for r in mk if int(r["cell"]) in (1, 2, 3, 4)]
y2 = np.array([math.log10((float(r["false_fungal_reads"]) + 0.5)
                          / float(r["total_reads"]) * 1e6) for r in allc])
sc2 = np.array([0.0 if r["index"] == "fungal-only" else 1.0 for r in allc])
cf2 = np.array([0.0 if float(r["confidence"]) == 0.0 else 1.0 for r in allc])
runs2 = sorted({r["run"] for r in allc})
X2 = np.column_stack([np.ones(len(y2)), sc2, cf2, sc2 * cf2]
                     + [np.array([1.0 if r["run"] == rn else 0.0 for r in allc])
                        for rn in runs2[1:]])
b2, *_ = np.linalg.lstsq(X2, y2, rcond=None)
r2 = y2 - X2 @ b2
df2 = len(y2) - np.linalg.matrix_rank(X2)
se2 = np.sqrt(np.diag(np.linalg.pinv(X2.T @ X2)) * (r2 @ r2) / df2)
for i, lab in enumerate(("reference scope", "confidence 0.1", "interaction"), start=1):
    print("   %-20s %9.2f-fold   P = %.4g" % (lab, 10 ** -b2[i],
                                              2 * stats.t.sf(abs(b2[i] / se2[i]), df2)))
chk("zero-preserving interaction, fold", 10 ** -b2[3], 11.792, 0.02)
chk("zero-preserving interaction exceeds the reported one",
    float(10 ** -b2[3] > 10 ** -beta[3]), 1.0, 0.0)

print("\n%d checks, %d mismatches" % (nchk[0], len(fails)))
if fails:
    raise SystemExit("mismatched: " + ", ".join(fails))
