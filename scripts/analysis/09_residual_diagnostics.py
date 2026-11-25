"""Residual normality and equality of variance for the season by substrate model,
and refits of the same contrast on ranks and on site-demeaned differences, which
assume neither.

    python3 09_residual_diagnostics.py [table_dir]

Appends the diagnostics to Table S10. Exit status is 0 only if every check passes.
"""
import csv, os, sys
import numpy as np
from scipy import stats

D = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "..", "tables")
P = lambda n: os.path.join(D, n)
SE = ["Winter", "Spring", "Summer", "Autumn"]
FAIL = []
NCHK = [0]


def chk(label, got, want, tol=5e-4):
    NCHK[0] += 1
    ok = abs(got - want) <= tol
    print("  %-56s %10.5g  want %-10.5g %s" % (label, got, want, "ok" if ok else "FAILED"))
    if not ok:
        FAIL.append(label)


B = list(csv.DictReader(open(P("TableS1_blocks.tsv")), delimiter="\t"))
d = np.array([float(r["d"]) for r in B])
sea = np.array([r["season"] for r in B])
site = np.array([r["site"] for r in B])
sites = sorted(set(site))
assert len(B) == 60 and len(sites) == 15

X = np.column_stack([np.ones(60)]
                    + [(sea == s).astype(float) for s in SE[1:]]
                    + [(site == s).astype(float) for s in sites[1:]])
b, *_ = np.linalg.lstsq(X, d, rcond=None)
res = d - X @ b
dfres = 60 - np.linalg.matrix_rank(X)

sw = stats.shapiro(res)
dp = stats.normaltest(res)
bf = stats.levene(*[res[sea == s] for s in SE], center="median")
ba = stats.bartlett(*[res[sea == s] for s in SE])
sd_raw = {s: float(d[sea == s].std(ddof=1)) for s in SE}
sd_res = {s: float(res[sea == s].std(ddof=1)) for s in SE}

print("=== residuals of the site-blocked season model, n = 60, df = %d ===" % dfres)
print("  Shapiro-Wilk            W = %.4f   P = %.4f" % (sw.statistic, sw.pvalue))
print("  D'Agostino-Pearson      K2 = %.3f   P = %.4f" % (dp.statistic, dp.pvalue))
print("  skewness %+.3f   excess kurtosis %+.3f" % (stats.skew(res), stats.kurtosis(res)))
print("  Brown-Forsythe          W = %.4f   P = %.4f" % (bf.statistic, bf.pvalue))
print("  Bartlett                T = %.4f   P = %.4f" % (ba.statistic, ba.pvalue))
print("\n  per-season spread, raw differences against residuals")
for s in SE:
    print("    %-7s SD(d) = %.3f   SD(residual) = %.3f   skew(d) = %+.3f"
          % (s, sd_raw[s], sd_res[s], stats.skew(d[sea == s])))
vr_raw = (max(sd_raw.values()) / min(sd_raw.values())) ** 2
vr_res = (max(sd_res.values()) / min(sd_res.values())) ** 2
print("  variance ratio, raw %.2f   residual %.2f" % (vr_raw, vr_res))

# distribution-free refits of the same contrast
rk = stats.rankdata(d)
br, *_ = np.linalg.lstsq(X, rk, rcond=None)
rr = rk - X @ br
X0 = np.column_stack([np.ones(60)] + [(site == s).astype(float) for s in sites[1:]])
b0, *_ = np.linalg.lstsq(X0, rk, rcond=None)
r0 = rk - X0 @ b0
Frank = ((r0 @ r0 - rr @ rr) / 3) / ((rr @ rr) / dfres)
Prank = stats.f.sf(Frank, 3, dfres)
dm = d.copy()
for s in sites:
    m = site == s
    dm[m] = d[m] - d[m].mean()
kw = stats.kruskal(*[dm[sea == s] for s in SE])
fr = stats.friedmanchisquare(*[np.array([d[(sea == s) & (site == st)][0] for st in sites])
                               for s in SE])
print("\n=== the same contrast without the normality assumption ===")
print("  rank-transformed F        F(3,%d) = %.3f   P = %.4f" % (dfres, Frank, Prank))
print("  Kruskal-Wallis, site-demeaned   H = %.3f   P = %.5f" % (kw.statistic, kw.pvalue))
print("  Friedman, as reported     chi2 = %.3f   P = %.4f" % (fr.statistic, fr.pvalue))

print("\n=== checks ===")
chk("Shapiro-Wilk W", float(sw.statistic), 0.9852, 5e-4)
chk("Shapiro-Wilk P", float(sw.pvalue), 0.681, 2e-3)
chk("Brown-Forsythe W", float(bf.statistic), 0.5054, 5e-4)
chk("Brown-Forsythe P", float(bf.pvalue), 0.680, 2e-3)
chk("raw variance ratio across seasons", vr_raw, 3.84, 0.01)
chk("residual variance ratio across seasons", vr_res, 2.10, 0.01)
chk("rank-transformed F", Frank, 6.530, 5e-3)
chk("rank-transformed P", Prank, 0.0010, 5e-4)
chk("Friedman P reproduces the reported value", float(fr.pvalue), 0.0030, 5e-4)
chk("residuals pass normality at 0.05", float(sw.pvalue > 0.05), 1.0, 0.0)
chk("residuals pass equal variance at 0.05", float(bf.pvalue > 0.05), 1.0, 0.0)
chk("every distribution-free refit stays below 0.05",
    float(max(Prank, kw.pvalue, fr.pvalue) < 0.05), 1.0, 0.0)

rows = [
    ("residual_normality", "Shapiro-Wilk on the residuals of the site-blocked season model",
     "W=%.4f" % sw.statistic, "", str(dfres), "%.4f" % sw.pvalue,
     "skewness %+.3f, excess kurtosis %+.3f" % (stats.skew(res), stats.kurtosis(res))),
    ("residual_homoscedasticity", "Brown-Forsythe across the four seasons",
     "W=%.4f" % bf.statistic, "", "3,56", "%.4f" % bf.pvalue,
     "per-season residual SD %s; variance ratio %.2f against %.2f for the raw differences"
     % ("/".join("%.2f" % sd_res[s] for s in SE), vr_res, vr_raw)),
    ("interaction_rank_transformed", "season x substrate on ranks of the block differences",
     "F=%.3f" % Frank, "", "3,%d" % dfres, "%.4f" % Prank,
     "distribution-free refit of the primary test"),
    ("interaction_kruskal_sitedemeaned", "Kruskal-Wallis on site-demeaned block differences",
     "H=%.3f" % kw.statistic, "", "3", "%.5f" % kw.pvalue,
     "distribution-free refit of the primary test"),
]
old = list(csv.reader(open(P("TableS10_slope_sensitivity.tsv")), delimiter="\t"))
head, body = old[0], [r for r in old[1:] if r and r[0]]
keep = [r for r in body if r[0] not in {x[0] for x in rows}]
with open(P("TableS10_slope_sensitivity.tsv"), "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(head)
    for r in keep:
        w.writerow(r)
    for r in rows:
        w.writerow(list(r))
print("\n  Table S10 now holds %d analyses" % (len(keep) + len(rows)))

print("\n%d checks, %d mismatches" % (NCHK[0], len(FAIL)))
if FAIL:
    raise SystemExit("FAILED: %d check(s): %s" % (len(FAIL), ", ".join(FAIL)))

