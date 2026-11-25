"""Recompute every value in the Results subsection on water temperature, including
the within-site slope and its distribution-free checks, the leave-one-site-out
and river-clustered refits, the covariate screen against per-variable detectable
bounds, the multivariable fit, and the within-season and site-demeaned models.

Table S9 reports two standardizations that are not interchangeable. Univariate
rows divide by the raw pooled standard deviation of the predictor across all 60
blocks, multivariable rows by the within-site standard deviation after demeaning.
For water temperature these are 1.6047 and 1.3902 log2 units per standard
deviation. The adjusted coefficient of 1.4055 is quoted against the demeaned
value, and the script prints both.

    python3 04_temperature_slopes.py [table_dir]

Exit status is 0 only if every value matches.
"""
import csv, math, os, sys
import numpy as np
from scipy import stats

D = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "tables")
def path(n):
    for c in (os.path.join(D, n), n):
        if os.path.exists(c): return c
    raise SystemExit("cannot find %s" % n)

B = list(csv.DictReader(open(path("TableS1_blocks.tsv")), delimiter="\t"))
C = {(r["season"], r["site"]): r for r in csv.DictReader(open(path("bao2026_covariates.tsv")), delimiter="\t")}
assert len(B) == 60 and len(C) == 60

d = np.array([float(r["d"]) for r in B])
site = np.array([r["site"] for r in B])
sea = np.array([r["season"] for r in B])
key = [(r["season"], r["site"]) for r in B]
sites = sorted(set(site))
SEASONS = ["Winter", "Spring", "Summer", "Autumn"]
NS = len(sites)
assert NS == 15
def col(n):
    """Missing entries come back as NaN. River length and catchment area are not
    recorded for every site, which is why Table S11 reports n = 48 and n = 56 for
    those two rows and n = 60 for the rest."""
    out = np.empty(len(key))
    for i, k in enumerate(key):
        v = C[k][n].strip()
        out[i] = float(v) if v else np.nan
    return out
tp = col("water_temp_C")

FAIL = []
NCHK = [0]
def chk(label, got, want, tol):
    NCHK[0] += 1
    ok = abs(float(got) - float(want)) <= tol
    if not ok: FAIL.append(label)
    print("  %-56s %12.5f  manuscript %-11s %s" % (label, got, "%g" % want, "ok" if ok else "MISMATCH"))

def demean(v):
    o = v.astype(float).copy()
    for s in sites:
        m = site == s
        o[m] = v[m] - v[m].mean()
    return o

def slopes(x, y=d, mask=None):
    """One slope per site. Sites where the predictor does not vary carry no
    information about a slope and are dropped, which is why Table S9 reports
    n = 13 for salinity and n = 15 for every other variable."""
    m0 = np.ones(len(x), bool) if mask is None else mask
    out = []
    for s in sites:
        m = m0 & (site == s)
        if x[m].std() == 0.0:
            continue
        out.append(stats.linregress(x[m], y[m]).slope)
    return np.array(out)

def ols(X, y, df=None):
    """df is given explicitly where the data were site-demeaned before fitting,
    because the 15 absorbed site means cost degrees of freedom that the design
    matrix no longer shows."""
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    r = y - X @ b
    if df is None:
        df = len(y) - np.linalg.matrix_rank(X)
    se = np.sqrt((r @ r / df) * np.diag(np.linalg.pinv(X.T @ X)))
    return b, se, df

S4 = np.column_stack([(sea == s).astype(float) for s in SEASONS])
SITE = np.column_stack([(site == s).astype(float) for s in sites])

print("\n=== 1. within-site temperature slope, per degree C ===")
sl = slopes(tp)
t = stats.ttest_1samp(sl, 0.0)
chk("mean slope (log2 units per C)", sl.mean(), 0.413, 5e-4)
chk("  fold per C", 2 ** sl.mean(), 1.33, 5e-3)
lo, hi = stats.t.interval(0.95, NS - 1, loc=sl.mean(), scale=stats.sem(sl))
chk("  95% CI lower", lo, 0.127, 1e-3)
chk("  95% CI upper", hi, 0.700, 1e-3)
chk("  t statistic", t.statistic, 3.09, 5e-3)
chk("  P", t.pvalue, 0.0079, 5e-5)
chk("  sites with a positive slope", (sl > 0).sum(), 13, 0)
chk("  exact signed-rank P", stats.wilcoxon(sl, method="exact").pvalue, 0.0034, 5e-5)
chk("  exact sign test P", stats.binomtest(int((sl > 0).sum()), NS, 0.5).pvalue, 0.0074, 5e-5)
chk("  largest leave-one-site-out P", max(stats.ttest_1samp(np.delete(sl, i), 0.0).pvalue for i in range(NS)), 0.017, 0.0005)
chk("  water temperature minimum (C)", tp.min(), 18.0, 1e-9)
chk("  water temperature maximum (C)", tp.max(), 31.9, 1e-9)

rivers = {}
for s in sites: rivers.setdefault(s[0], []).append(sites.index(s))
rs = np.array([sl[v].mean() for v in map(np.array, rivers.values())])
tr = stats.ttest_1samp(rs, 0.0)
chk("river-averaged slope (10 rivers)", rs.mean(), 0.426, 5e-4)
chk("  t statistic", tr.statistic, 2.75, 5e-3)
chk("  df", len(rs) - 1, 9, 0)
chk("  P", tr.pvalue, 0.0225, 5e-5)

for drop, wb, wt, wp in (("Summer", 0.311, None, 0.020), ("Winter", 0.879, 3.31, None)):
    s2 = slopes(tp, mask=(sea != drop))
    tt = stats.ttest_1samp(s2, 0.0)
    chk("slope with %s excluded" % drop, s2.mean(), wb, 5e-4)
    if wt is not None: chk("  t statistic", tt.statistic, wt, 5e-3)
    if wp is not None: chk("  P", tt.pvalue, wp, 5e-5)
    if drop == "Winter":
        a, b2 = stats.t.interval(0.95, NS - 1, loc=s2.mean(), scale=stats.sem(s2))
        chk("  95% CI lower (three warm seasons)", a, 0.310, 1e-3)
        chk("  95% CI upper (three warm seasons)", b2, 1.447, 1e-3)

print("\n=== 2. screen across the 12 measured variables ===")
TWELVE = ["water_temp_C", "precip_mm_day", "velocity_m_s", "turbidity_NTU", "salinity_psu", "pH",
          "DO_mg_L", "COD_mg_L", "TN_ug_L", "TP_ug_L", "SO4_ug_L", "total_MPs_per_L"]
chk("number of measured variables screened", len(TWELVE), 12, 0)
K = (stats.t.ppf(0.975, NS - 1) + stats.t.ppf(0.80, NS - 1)) / math.sqrt(NS)
chk("design MDES, in SD of the site slope distribution", K, 0.78, 5e-3)
chk("Holm threshold across 12 tests", 0.05 / 12, 0.0042, 5e-5)
res = {}
for n in TWELVE:
    v = col(n)
    z = slopes(v) * v.std(ddof=1)
    res[n] = (z.mean(), K * float(np.std(z, ddof=1)), stats.ttest_1samp(z, 0.0).pvalue)
order = sorted(res, key=lambda n: res[n][2])
chk("smallest P belongs to water temperature", order[0] == "water_temp_C", True, 0)
chk("  its P exceeds the Holm threshold", res["water_temp_C"][2] > 0.05 / 12, True, 0)
chk("second smallest P (average daily precipitation)", res[order[1]][2], 0.137, 5e-4)
chk("  is precip_mm_day", order[1] == "precip_mm_day", True, 0)
chk("  its observed slope per SD", res["precip_mm_day"][0], 2.48, 5e-3)
chk("  its detectable bound per SD", res["precip_mm_day"][1], 4.73, 5e-3)
chk("ambient MP concentration, slope per SD", res["total_MPs_per_L"][0], -1.90, 5e-3)
chk("  P", res["total_MPs_per_L"][2], 0.25, 5e-3)
chk("  detectable bound per SD", res["total_MPs_per_L"][1], 4.80, 5e-3)

print("\n=== 3. the two standardizations, and adjustment for sewage loading ===")
chk("water temperature per RAW POOLED SD (Table S9 univariate)", res["water_temp_C"][0], 1.6047, 5e-4)
NAMES = ["water_temp_C", "COD_mg_L", "TN_ug_L", "TP_ug_L"]
Z = np.column_stack([demean(col(n)) / demean(col(n)).std(ddof=1) for n in NAMES])
y = demean(d)
DFR = 60 - 4 - 15
b1, se1, _ = ols(np.column_stack([np.ones(60), Z[:, :1]]), y, df=60 - 1 - 15)
chk("water temperature per DEMEANED SD, unadjusted", b1[1], 1.390, 5e-4)
b, se, _ = ols(np.column_stack([np.ones(60), Z]), y, df=DFR)
chk("  adjusted", b[1], 1.405, 5e-4)
chk("  percent change on adjustment", 100 * (b[1] / b1[1] - 1), 1.1, 0.1)
chk("  t statistic", b[1] / se[1], 2.94, 5e-3)
chk("  residual df", DFR, 41, 0)
chk("  P", 2 * stats.t.sf(abs(b[1] / se[1]), DFR), 0.0054, 5e-5)
for i, (n, wp) in enumerate(zip(NAMES[1:], (0.21, 0.59, 0.32)), start=2):
    chk("  %s P" % n, 2 * stats.t.sf(abs(b[i] / se[i]), DFR), wp, 5e-3)
vifs = []
for i in range(4):
    o = [j for j in range(4) if j != i]
    bb, *_ = np.linalg.lstsq(np.column_stack([np.ones(60), Z[:, o]]), Z[:, i], rcond=None)
    r = Z[:, i] - np.column_stack([np.ones(60), Z[:, o]]) @ bb
    cen = Z[:, i] - Z[:, i].mean()
    vifs.append(1.0 / (1.0 - (1.0 - r @ r / (cen @ cen))))
chk("  smallest variance inflation factor", min(vifs), 1.39, 5e-3)
chk("  largest variance inflation factor", max(vifs), 2.08, 5e-3)
for n, wr in (("COD_mg_L", -0.446), ("TP_ug_L", -0.424)):
    r, p = stats.pearsonr(demean(tp), demean(col(n)))
    chk("within-site r, temperature vs %s" % n, r, wr, 5e-4)
    chk("  P below 0.001", p < 0.001, True, 0)

print("\n=== 4. collinearity with water temperature, across the 60 blocks ===")
ALL = TWELVE[1:] + ["altitude_m", "river_length_km", "catchment_km2", "population"]
rh = {}
for n in ALL:
    v = col(n)
    m = np.isfinite(v)
    rh[n] = stats.spearmanr(tp[m], v[m])
top = max(rh, key=lambda n: abs(rh[n].statistic))
chk("largest |rho| with water temperature", abs(rh[top].statistic), 0.42, 5e-3)
chk("  it is average daily precipitation", top == "precip_mm_day", True, 0)
chk("  its P", rh[top].pvalue, 0.0009, 5e-5)
chk("  none reaches 0.7", max(abs(rh[n].statistic) for n in ALL) < 0.7, True, 0)
chk("upstream population rho", rh["population"].statistic, 0.158, 5e-4)
chk("  P", rh["population"].pvalue, 0.229, 5e-4)

print("\n=== 5. antecedent rainfall (see 05_antecedent_rainfall.py for the index itself) ===")
import datetime as dt
rain = {r["date"]: float(r["territory_mean_mm"]) for r in csv.DictReader(open(path("HKO_daily_rainfall_2023.tsv")), delimiter="\t")}
dates = {(r["season"], r["site"]): r["date"] for r in csv.DictReader(open(path("block_sampling_dates.tsv")), delimiter="\t")}
def ante(k, win):
    day = dt.date.fromisoformat(dates[k])
    return sum(rain.get((day - dt.timedelta(days=j)).isoformat(), 0.0) for j in range(1, win + 1))
for win, wobs, wmd, wp in ((3, 1.175, 1.128, 0.0073), (7, 1.584, 1.558, 0.0085)):
    v = np.array([ante(k, win) for k in key])
    z = slopes(v) * v.std(ddof=1)
    chk("%d-day rainfall, slope per SD" % win, z.mean(), wobs, 5e-4)
    chk("  its detectable bound per SD", K * float(np.std(z, ddof=1)), wmd, 5e-4)
    chk("  P", stats.ttest_1samp(z, 0.0).pvalue, wp, 5e-5)
for win, wb, wp in ((1, -0.071, 0.435), (3, 0.032, 0.309), (7, 0.017, 0.480)):
    v = np.array([ante(k, win) for k in key])
    bb, ss, df = ols(np.column_stack([v, S4]), d)
    chk("%d-day rainfall given season, per mm" % win, bb[0], wb, 5e-4)
    chk("  P", 2 * stats.t.sf(abs(bb[0] / ss[0]), df), wp, 5e-4)
def interaction(y):
    X = np.column_stack([np.ones(60), S4[:, 1:], SITE[:, 1:]])
    b2, *_ = np.linalg.lstsq(X, y, rcond=None); r2 = y - X @ b2
    X0 = np.column_stack([np.ones(60), SITE[:, 1:]])
    b0, *_ = np.linalg.lstsq(X0, y, rcond=None); r0 = y - X0 @ b0
    dfres = 60 - np.linalg.matrix_rank(X)
    Fv = ((r0 @ r0 - r2 @ r2) / 3) / (r2 @ r2 / dfres)
    return Fv, dfres, stats.f.sf(Fv, 3, dfres)
for win, wF, wp in ((1, 6.18, 0.0014), (3, 4.29, 0.0099), (7, 2.92, 0.045)):
    v = np.array([ante(k, win) for k in key])
    bb, *_ = np.linalg.lstsq(np.column_stack([np.ones(60), v]), d, rcond=None)
    Fv, dfres, pv = interaction(d - np.column_stack([np.ones(60), v]) @ bb)
    chk("interaction with %d-day rainfall removed, F" % win, Fv, wF, 5e-3)
    chk("  residual df", dfres, 42, 0)
    chk("  P", pv, wp, 5e-4)
r7 = np.array([ante(k, 7) for k in key])
chk("7-day rainfall vs water temperature, r", stats.pearsonr(r7, tp).statistic, 0.66, 5e-3)
for s, w in (("Spring", 42.4), ("Summer", 55.4)):
    chk("mean 7-day antecedent rainfall, %s (mm)" % s, np.mean([ante(k, 7) for k in key if k[0] == s]), w, 0.05)
chk("spring to summer gap in mean enrichment (log2)", d[sea == "Summer"].mean() - d[sea == "Spring"].mean(), 3.31, 5e-3)
su = [k for k in key if k[0] == "Summer"]
lag = np.array([(dt.date.fromisoformat(dates[k]) - dt.date(2023, 7, 29)).days for k in su])
ds = np.array([d[key.index(k)] for k in su])
chk("summer contrast vs days since the heaviest preceding rainfall, r", stats.pearsonr(lag, ds).statistic, -0.21, 5e-3)
chk("  P", stats.pearsonr(lag, ds).pvalue, 0.45, 5e-3)

print("\n=== 6. temperature cannot be separated from season ===")
for s, wr, wp in (("Winter", 0.16, 0.58), ("Spring", 0.17, 0.55), ("Summer", 0.03, 0.93), ("Autumn", 0.09, 0.74)):
    m = sea == s
    r, p = stats.pearsonr(tp[m], d[m])
    chk("%s within-season r" % s, r, wr, 5e-3)
    chk("  P", p, wp, 5e-3)
bb, ss, df = ols(np.column_stack([tp, S4]), d)
tv = bb[0] / ss[0]
chk("pooled within-season slope", bb[0], 0.179, 5e-4)
chk("  t statistic", tv, 0.73, 5e-3)
chk("  residual df (60 - 1 slope - 4 season levels)", df, 55, 0)
chk("  P", 2 * stats.t.sf(abs(tv), df), 0.47, 5e-4)
crit = stats.t.ppf(0.975, df)
chk("  95% CI lower", bb[0] - crit * ss[0], -0.314, 1e-3)
chk("  95% CI upper", bb[0] + crit * ss[0], 0.671, 1e-3)
for lab, X, wb, wt, wdf, wp in (
        ("site-demeaned, temperature alone", np.column_stack([tp, SITE]), 0.376, 3.58, 44, 0.00084),
        ("site-demeaned, temperature given season", np.column_stack([tp, SITE, S4[:, 1:]]), 0.033, 0.09, 41, 0.93)):
    bb, ss, df = ols(X, d)
    tv = bb[0] / ss[0]
    chk(lab, bb[0], wb, 5e-4)
    chk("  t statistic", tv, wt, 5e-3)
    chk("  residual df", df, wdf, 0)
    chk("  P", 2 * stats.t.sf(abs(tv), df), wp, max(5e-6, wp * 0.02))
chk("winter maximum temperature (C)", tp[sea == "Winter"].max(), 23.0, 1e-9)
chk("coldest observation outside winter (C)", tp[sea != "Winter"].min(), 23.4, 1e-9)
sp, au, sm = (np.array([d[key.index((s, k))] for k in sites]) for s in ("Spring", "Autumn", "Summer"))
tsp, tau = (np.array([tp[key.index((s, k))] for k in sites]) for s in ("Spring", "Autumn"))
chk("spring to autumn gap in mean temperature (C)", abs(tsp.mean() - tau.mean()), 0.31, 5e-3)
chk("spring to autumn gap in mean enrichment (log2)", abs(sp.mean() - au.mean()), 1.20, 5e-3)
chk("  as a percentage of the summer mean", 100 * abs(sp.mean() - au.mean()) / sm.mean(), 40.0, 0.6)
chk("  paired exact signed-rank P", stats.wilcoxon(au - sp, method="exact").pvalue, 0.28, 5e-3)

print("\n%d checks, %d mismatches" % (NCHK[0], len(FAIL)))
if FAIL:
    for f in FAIL: print("   MISMATCH:", f)
    sys.exit(1)
