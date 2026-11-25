"""Recompute the headline values of the paper from the deposited tables and compare
each with the value in the text. Covers the endpoint identity, the season by
substrate interaction, the leave-one-site-out range, the within-deposit refit,
the Pichia-removed interaction, the temperature slope, the per-season
Hodges-Lehmann estimates and Holm-adjusted P values, the exact Walsh interval,
the false-positive floor series, the cross-study meta-analysis and the mock
community factorial.

    python3 01_check_main_values.py [table_dir]

Exit status is 0 only if every value matches.
"""
import csv, math, os, sys, collections
import numpy as np
from scipy import stats

D = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "tables")
def T(n):
    p = os.path.join(D, n)
    if not os.path.exists(p): raise SystemExit("missing %s" % p)
    return list(csv.DictReader(open(p), delimiter="\t"))

FAIL = []
NCHK = [0]
def chk(label, got, want, tol=5e-4):
    NCHK[0] += 1
    ok = (abs(got - want) <= tol) if isinstance(want, float) else (got == want)
    print("  %-52s %14s  vs %-14s %s"
          % (label, ("%.5g" % got) if isinstance(got, float) else got,
             ("%.5g" % want) if isinstance(want, float) else want,
             "MATCH" if ok else "*** MISMATCH ***"))
    if not ok: FAIL.append(label)

B = T("TableS1_blocks.tsv")
d = np.array([float(r["d"]) for r in B]); sea = np.array([r["season"] for r in B])
site = np.array([r["site"] for r in B]); tp = np.array([float(r["temp"]) for r in B])
SE = ["Winter", "Spring", "Summer", "Autumn"]; SI = sorted(set(site))

print("\n== endpoint identity ==")
dev = 0.0
for r in B:
    for arm in ("MP", "NP"):
        v = float(r["fungal_" + arm]) - (float(r["total_" + arm]) / 1e6) * 0.16
        dev = max(dev, abs(math.log2(v / float(r["total_" + arm])) - float(r["log2_" + arm])))
chk("max deviation, recomputed vs deposited log2", dev, 0.0, 1e-4)

def interaction(y, s, si, k=3):
    X = np.column_stack([np.ones(len(y))] + [(s == x).astype(float) for x in sorted(set(s))[1:]]
                        + [(si == x).astype(float) for x in sorted(set(si))[1:]])
    b, *_ = np.linalg.lstsq(X, y, rcond=None); r = y - X @ b; rss = r @ r
    dfr = len(y) - X.shape[1]
    X0 = np.column_stack([np.ones(len(y))] + [(si == x).astype(float) for x in sorted(set(si))[1:]])
    b0, *_ = np.linalg.lstsq(X0, y, rcond=None); r0 = y - X0 @ b0
    F = ((r0 @ r0 - rss) / k) / (rss / dfr)
    return F, dfr, stats.f.sf(F, k, dfr)

print("\n== primary interaction ==")
F, dfr, p = interaction(d, sea, site)
chk("season x substrate F(3,42)", F, 6.064, 0.002)
chk("  P", p, 0.0016, 2e-4)
fr = stats.friedmanchisquare(*[d[sea == s] for s in SE])
chk("Friedman Q", fr.statistic, 13.96, 0.01)
chk("  P", fr.pvalue, 0.0030, 2e-4)

print("\n== leave-one-site-out on the interaction ==")
ps = [interaction(d[site != s], sea[site != s], site[site != s])[2] for s in SI]
chk("min P", min(ps), 0.00009, 2e-5); chk("max P", max(ps), 0.0079, 2e-4)

print("\n== within-deposit (spring/summer/autumn) ==")
m = sea != "Winter"
F2, dfr2, p2 = interaction(d[m], sea[m], site[m], k=2)
chk("F(2,28)", F2, 5.234, 0.01); chk("  P", p2, 0.0117, 5e-4)

print("\n== Pichia-removed interaction, from the genus matrix ==")
G = T("TableS25_genus_matrix_180.tsv")
pich = collections.defaultdict(int); ftot = {}
for r in G:
    ftot[r["run"]] = int(r["fungal_clade_reads"])
    if r["genus"] == "Pichia": pich[r["run"]] = int(r["reads"])
tot = {}
for r in B:
    tot[r["run_MP"]] = int(r["total_MP"]); tot[r["run_NP"]] = int(r["total_NP"])
def ep(run, drop=False):
    Tp = tot[run]; f = ftot.get(run, 0) - (pich.get(run, 0) if drop else 0)
    v = f - (Tp / 1e6) * 0.16
    return math.log2(v / Tp)
dd = np.array([ep(r["run_MP"], True) - ep(r["run_NP"], True) for r in B])
F3, _, p3 = interaction(dd, sea, site)
chk("Pichia-removed F(3,42)", F3, 4.790, 0.002); chk("  P", p3, 0.0058, 2e-4)
chk("summer mean, Pichia removed", float(dd[sea == "Summer"].mean()), 2.398, 0.002)

print("\n== temperature slope ==")
sl = [stats.linregress(tp[site == s], d[site == s]).slope for s in SI]
t = stats.ttest_1samp(sl, 0.0)
chk("mean within-site slope", float(np.mean(sl)), 0.4133, 2e-4)
chk("  t(14)", float(t.statistic), 3.093, 0.002)
chk("  P", float(t.pvalue), 0.0079, 2e-4)
chk("  sites positive", int(sum(1 for x in sl if x > 0)), 13)

print("\n== per-season Hodges-Lehmann and Holm ==")
def hl(v):
    w = [(v[i] + v[j]) / 2 for i in range(len(v)) for j in range(i, len(v))]
    return float(np.median(w))
raw = {s: stats.wilcoxon([x for x in d[sea == s] if x != 0], method="exact").pvalue for s in SE}
adj, run = {}, 0.0
for i, s in enumerate(sorted(SE, key=lambda z: raw[z])):
    run = max(run, (4 - i) * raw[s]); adj[s] = min(1.0, run)
chk("summer HL (log2)", hl(list(d[sea == "Summer"])), 2.4514, 2e-4)
chk("winter HL (log2)", hl(list(d[sea == "Winter"])), -1.1731, 2e-4)
chk("summer Holm P", adj["Summer"], 0.0105, 2e-4)
chk("winter Holm P", adj["Winter"], 0.0905, 2e-4)

print("\n== exact Walsh interval for the summer estimate ==")
def walsh(v):
    return sorted((v[i] + v[j]) / 2 for i in range(len(v)) for j in range(i, len(v)))
def srank_cdf(n):
    N = n * (n + 1) // 2
    c = np.zeros(N + 1); c[0] = 1.0
    for r in range(1, n + 1):
        nx = c.copy(); nx[r:] += c[:N + 1 - r]; c = nx
    return np.cumsum(c) / c.sum()
sv = list(d[sea == "Summer"]); w = walsh(sv); N = len(w); cdf = srank_cdf(len(sv))
ci_c = max(k for k in range(0, N // 2) if 1.0 - 2.0 * cdf[k] >= 0.95)
chk("summer CI lower (fold)", 2 ** w[ci_c], 1.620, 0.002)
chk("summer CI upper (fold)", 2 ** w[N - 1 - ci_c], 39.035, 0.02)
chk("  exact coverage", 1.0 - 2.0 * cdf[ci_c], 0.9521, 2e-4)
_pout = stats.wilcoxon([x - (w[ci_c] - 1e-9) for x in sv], method="exact").pvalue
chk("  coverage equals 1 minus P just outside", 1.0 - _pout, 1.0 - 2.0 * cdf[ci_c], 1e-9)

print("\n== floor sensitivity ==")
for fl, want in ((0.0, 6.209), (0.30, 5.892), (0.45, 5.601)):
    y = []
    for r in B:
        v = [float(r["fungal_" + a]) - (float(r["total_" + a]) / 1e6) * fl for a in ("MP", "NP")]
        y.append(math.log2(v[0] / float(r["total_MP"])) - math.log2(v[1] / float(r["total_NP"])))
    chk("interaction F at floor %.2f" % fl, interaction(np.array(y), sea, site)[0], want, 0.002)

print("\n== meta-analysis (Figure 4) ==")
C = [r for r in T("TableS3_crossstudy.tsv") if int(r["n_surface_control"]) >= 3]
def dl(rows):
    y = np.array([float(r["diff_log2"]) for r in rows]); se = np.array([float(r["se"]) for r in rows])
    w = 1 / se ** 2; fe = (w * y).sum() / w.sum(); Q = (w * (y - fe) ** 2).sum(); dq = len(y) - 1
    Cc = w.sum() - (w ** 2).sum() / w.sum(); t2 = max(0.0, (Q - dq) / Cc)
    wr = 1 / (se ** 2 + t2)
    return float((wr * y).sum() / wr.sum()), (max(0.0, (Q - dq) / Q) * 100 if Q > 0 else 0.0), wr
ext = [r for r in C if "Hong Kong" not in r["system"]]
e1, i1, _ = dl(ext); e2, i2, wr = dl(C)
chk("external-only pooled (fold)", 2 ** e1, 0.834, 0.002); chk("  I2 (%)", i1, 77.6, 0.3)
chk("all-nine pooled (fold)", 2 ** e2, 0.893, 0.002); chk("  I2 (%)", i2, 78.7, 0.3)
hk = [i for i, r in enumerate(C) if "Hong Kong" in r["system"]]
chk("own deposits' random-effects weight (%)", 100 * wr[hk].sum() / wr.sum(), 18.14, 0.05)

print("\n== mock 2x2, blocked on library ==")
Mk = [r for r in T("TableS7_mock2x2.tsv")
      if float(r["min_hit_groups"]) == 2 and float(r["false_per_million"]) > 0]
y = np.array([math.log10(float(r["false_per_million"])) for r in Mk])
sc = np.array([0.0 if r["index"] == "fungal-only" else 1.0 for r in Mk])
cf = np.array([0.0 if float(r["confidence"]) == 0.0 else 1.0 for r in Mk])
lib = np.array([r["run"] for r in Mk]); U = sorted(set(lib))
X = np.column_stack([np.ones(len(y)), sc, cf, sc * cf]
                    + [(lib == u).astype(float) for u in U[1:]])
b, *_ = np.linalg.lstsq(X, y, rcond=None); res = y - X @ b; df = len(y) - X.shape[1]
se = np.sqrt(np.diag(np.linalg.pinv(X.T @ X)) * (res @ res) / df)
for i, (nm, wf, wp) in enumerate([("index scope", 664.1, 3.7e-10), ("confidence", 4.43, 0.023),
                                  ("interaction", 11.29, 0.0135)], start=1):
    chk("%s fold" % nm, 10 ** -b[i], wf, max(0.02, wf * 0.002))
    chk("  P", 2 * stats.t.sf(abs(b[i] / se[i]), df), wp, max(1e-10, wp * 0.05))

print("\n%d checks, %d mismatches" % (NCHK[0], len(FAIL)))
sys.exit(0 if not FAIL else 1)
