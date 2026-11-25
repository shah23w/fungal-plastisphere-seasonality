"""Build Table S13, the per-season summary of the microplastic minus natural
particle contrast, and recompute the leave-one-site-out F range recorded in
Table S10.

Effect sizes are Hodges-Lehmann estimates with exact distribution-free intervals.
For n paired differences there are N = n(n+1)/2 Walsh averages, and on the sorted
array the interval (w[c], w[N-1-c]) excludes the true shift with probability
2*P(W+ <= c), so its exact coverage is 1 - 2*P(W+ <= c). c is the largest index
whose coverage still reaches 1 - alpha. Coverage is reported alongside each
interval because at n = 15 no choice of c gives exactly 95%. Each interval is
checked by inverting the exact signed-rank test.

Per-arm medians are recomputed from the read counts in Table S2 rather than
differenced from the rounded cells of Table S1.

    python3 03_season_summary.py [table_dir]

Rewrites TableS13_season_summary.tsv. Exit status is 0 only if every check passes.
"""
import csv, math, os, sys
import numpy as np
from scipy import stats

D = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "..", "tables")
ALPHA = 0.05
SE = ["Winter", "Spring", "Summer", "Autumn"]
FAIL = []
NCHK = [0]


def chk(label, cond):
    NCHK[0] += 1
    print("  %-58s %s" % (label, "ok" if cond else "*** FAILED ***"))
    if not cond:
        FAIL.append(label)


def walsh(v):
    return sorted((v[i] + v[j]) / 2 for i in range(len(v)) for j in range(i, len(v)))


def hl(v):
    return float(np.median(walsh(v)))


def srank_cdf(n):
    """Exact null distribution of W+ for n non-zero differences."""
    N = n * (n + 1) // 2
    cnt = np.zeros(N + 1)
    cnt[0] = 1.0
    for r in range(1, n + 1):
        nxt = cnt.copy()
        nxt[r:] += cnt[:N + 1 - r]
        cnt = nxt
    return np.cumsum(cnt) / cnt.sum()


def hl_ci(v, alpha=ALPHA):
    w = walsh(v)
    n = len(v)
    N = n * (n + 1) // 2
    cdf = srank_cdf(n)
    c = max(k for k in range(0, N // 2) if 1.0 - 2.0 * cdf[k] >= 1.0 - alpha)
    return w[c], w[N - 1 - c], 1.0 - 2.0 * cdf[c], c


def exact_p(v, theta=0.0):
    x = [a - theta for a in v if a != theta]
    return stats.wilcoxon(x, method="exact").pvalue if len(x) > 1 else 1.0


B = list(csv.DictReader(open(os.path.join(D, "TableS1_blocks.tsv")), delimiter="\t"))
d = {s: [float(r["d"]) for r in B if r["season"] == s] for s in SE}
tp = {s: [float(r["temp"]) for r in B if r["season"] == s] for s in SE}
L2 = {a: {s: [float(r["log2_" + a]) for r in B if r["season"] == s] for s in SE}
      for a in ("MP", "NP", "RW")}

raw = {s: exact_p(d[s]) for s in SE}
adj, run = {}, 0.0
for i, s in enumerate(sorted(SE, key=lambda z: raw[z])):
    run = max(run, (len(SE) - i) * raw[s])
    adj[s] = min(1.0, run)

print("== per-season estimates ==")
rows = []
for s in SE:
    h = hl(d[s])
    lo, hi, cov, c = hl_ci(d[s])
    rows.append({
        "season": s, "n": len(d[s]),
        "mean_temp_C": "%.2f" % np.mean(tp[s]),
        "mean_d": "%.4f" % np.mean(d[s]),
        "fold": "%.3f" % 2 ** np.mean(d[s]),
        "hl_log2": "%.4f" % h,
        "hl_fold": "%.3f" % 2 ** h,
        "hl_ci_lo_fold": "%.3f" % 2 ** lo,
        "hl_ci_hi_fold": "%.3f" % 2 ** hi,
        "hl_ci_coverage": "%.4f" % cov,
        "walsh_c": c,
        "P_signed_rank_exact": "%.5f" % raw[s],
        "P_holm": "%.5f" % adj[s],
        "median_log2_MP": "%.4f" % np.median(L2["MP"][s]),
        "median_log2_NP": "%.4f" % np.median(L2["NP"][s]),
        "median_log2_RW": "%.4f" % np.median(L2["RW"][s]),
        "n_MP_gt_NP": sum(1 for x in d[s] if x > 0),
    })
    print("  %-7s HL %+.4f log2 = %.3f-fold  CI %.3f-%.3f (%.2f%% coverage, c=%d)"
          "  P %.5f  Holm %.5f  %d/%d positive"
          % (s, h, 2 ** h, 2 ** lo, 2 ** hi, cov * 100, c,
             raw[s], adj[s], rows[-1]["n_MP_gt_NP"], len(d[s])))

print("\n== interval verified by inverting the exact signed-rank test ==")
for s in SE:
    lo, hi, cov, c = hl_ci(d[s])
    eps = 1e-9
    chk("%s: inside both endpoints not rejected" % s,
        exact_p(d[s], lo + eps) > ALPHA and exact_p(d[s], hi - eps) > ALPHA)
    chk("%s: just outside both endpoints rejected" % s,
        exact_p(d[s], lo - eps) <= ALPHA and exact_p(d[s], hi + eps) <= ALPHA)
    chk("%s: stated coverage equals 1 minus the P just outside" % s,
        abs(cov - (1.0 - exact_p(d[s], lo - eps))) < 1e-9)

print("\n== values the manuscript states ==")
want = {"Winter": (0.443, 0.09045), "Spring": (0.729, 0.45862),
        "Summer": (5.470, 0.01050), "Autumn": (2.043, 0.45862)}
for s, r in zip(SE, rows):
    chk("%s HL fold" % s, abs(float(r["hl_fold"]) - want[s][0]) < 0.002)
    chk("%s Holm P" % s, abs(float(r["P_holm"]) - want[s][1]) < 2e-5)
summer = rows[SE.index("Summer")]
chk("summer interval is 1.62 to 39.0-fold",
    abs(float(summer["hl_ci_lo_fold"]) - 1.620) < 0.005
    and abs(float(summer["hl_ci_hi_fold"]) - 39.035) < 0.05)
chk("summer coverage is at least 95%", float(summer["hl_ci_coverage"]) >= 0.95)
chk("the next narrower interval would drop below 95%",
    hl_ci(d["Summer"], ALPHA)[3] + 1 not in
    [k for k in range(0, 60) if 1.0 - 2.0 * srank_cdf(15)[k] >= 0.95])

print("\n== winter-to-summer arm shifts the manuscript reports ==")
# The medians are recomputed from the raw read counts rather than differenced
# from the four-decimal cells of Table S1. The river-water shift sits at
# -0.8850 on those rounded cells, right on the two-decimal boundary, but its
# value from the counts is -0.884958, which rounds to -0.88. Rounding from the
# stored cell would print -0.89.
E0 = {r["run"]: r for r in csv.DictReader(open(os.path.join(D, "TableS2_detection_envelope.tsv")),
                                          delimiter="\t")}


FLOOR = 0.16


def endpoint(run):
    """The paper's endpoint, recomputed from the raw counts."""
    t, f = int(E0[run]["total_pairs"]), int(E0[run]["fungal_reads"])
    return math.log2((f - t / 1e6 * FLOOR) / t)


for a, wantv in (("MP", 2.78), ("NP", 0.10), ("RW", -0.88)):
    med = {s2: float(np.median([endpoint(r["run_" + a]) for r in B if r["season"] == s2]))
           for s2 in ("Winter", "Summer")}
    delta = med["Summer"] - med["Winter"]
    stored = float(np.median(L2[a]["Summer"]) - np.median(L2[a]["Winter"]))
    print("  %s %+.6f log2 = %.3f-fold   (from the stored cells %+.4f)"
          % (a, delta, 2 ** delta, stored))
    chk("%s shift is %+.2f to the reported precision" % (a, wantv),
        abs(delta - wantv) < 0.005)
chk("the river-water shift does not round to -0.89",
    round(float(np.median([endpoint(r["run_RW"]) for r in B if r["season"] == "Summer"])
                - np.median([endpoint(r["run_RW"]) for r in B if r["season"] == "Winter"])), 2)
    == -0.88)

print("\n== leave-one-site-out F range for Table S10 ==")
y = np.array([float(r["d"]) for r in B])
sea = np.array([r["season"] for r in B])
site = np.array([r["site"] for r in B])


def interaction(yy, ss, ii, k=3):
    X = np.column_stack([np.ones(len(yy))]
                        + [(ss == x).astype(float) for x in sorted(set(ss))[1:]]
                        + [(ii == x).astype(float) for x in sorted(set(ii))[1:]])
    b, *_ = np.linalg.lstsq(X, yy, rcond=None)
    r = yy - X @ b
    rss = r @ r
    dfr = len(yy) - X.shape[1]
    X0 = np.column_stack([np.ones(len(yy))]
                         + [(ii == x).astype(float) for x in sorted(set(ii))[1:]])
    b0, *_ = np.linalg.lstsq(X0, yy, rcond=None)
    r0 = yy - X0 @ b0
    F = ((r0 @ r0 - rss) / k) / (rss / dfr)
    return F, dfr, stats.f.sf(F, k, dfr)


Fs, Ps, FRs = [], [], []
for s in sorted(set(site)):
    m = site != s
    F, _, p = interaction(y[m], sea[m], site[m])
    Fs.append(F)
    Ps.append(p)
    sub = [r for r in B if r["site"] != s]
    FRs.append(stats.friedmanchisquare(
        *[[float(r["d"]) for r in sub if r["season"] == z] for z in SE]).pvalue)
print("  F %.4f to %.4f, P %.6f to %.6f, Friedman P %.5f to %.5f"
      % (min(Fs), max(Fs), min(Ps), max(Ps), min(FRs), max(FRs)))
chk("all 15 deletions keep P < 0.05", all(p < 0.05 for p in Ps))
chk("all 15 deletions keep Friedman P < 0.05", all(p < 0.05 for p in FRs))
chk("F range is 4.55 to 9.27",
    abs(min(Fs) - 4.5523) < 2e-3 and abs(max(Fs) - 9.2684) < 2e-3)

if FAIL:
    raise SystemExit("FAILED: %s" % FAIL)

out = os.path.join(D, "TableS13_season_summary.tsv")
cols = list(rows[0].keys())
with open(out, "w", newline="") as fh:
    wr = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
    wr.writeheader()
    wr.writerows(rows)
print("\nwrote %s with %d columns" % (out, len(cols)))
# Table S10 records this range, so check it rather than print advice about it.
want = "F=%.2f-%.2f" % (min(Fs), max(Fs))
row = [r for r in csv.reader(open(os.path.join(D, "TableS10_slope_sensitivity.tsv")),
                             delimiter="\t") if r and r[0] == "interaction_leave_one_site_out"]
chk("Table S10 records the leave-one-site-out F range as %s" % want,
    bool(row) and row[0][2] == want)

print("\n%d checks, %d mismatches" % (NCHK[0], len(FAIL)))
if FAIL:
    raise SystemExit("mismatched: " + ", ".join(FAIL))
