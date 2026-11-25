"""Derive the antecedent-rainfall variables and every rainfall statistic in the
Results and in Tables S9 and S10, from the daily Hong Kong Observatory index and
the per-block sampling dates.

Inputs: HKO_daily_rainfall_2023.tsv, block_sampling_dates.tsv, TableS1_blocks.tsv.

    python3 05_antecedent_rainfall.py [table_dir]
"""
import csv, math, os, sys, datetime as dt
import numpy as np
from scipy import stats

D = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
def path(n):
    for c in (os.path.join(D, n), os.path.join(D, "..", "..", "tables", n), n):
        if os.path.exists(c): return c
    raise SystemExit("cannot find %s" % n)

rain = {r["date"]: float(r["territory_mean_mm"])
        for r in csv.DictReader(open(path("HKO_daily_rainfall_2023.tsv")), delimiter="\t")}
dates = {(r["season"], r["site"]): r["date"]
         for r in csv.DictReader(open(path("block_sampling_dates.tsv")), delimiter="\t")}
B = list(csv.DictReader(open(path("TableS1_blocks.tsv")), delimiter="\t"))
print("rainfall days %d | blocks with dates %d | blocks %d" % (len(rain), len(dates), len(B)))

def antecedent(season, site, win):
    d = dt.date.fromisoformat(dates[(season, site)])
    return sum(rain.get((d - dt.timedelta(days=k)).isoformat(), 0.0) for k in range(1, win + 1))

d = np.array([float(r["d"]) for r in B])
sea = np.array([r["season"] for r in B])
site = np.array([r["site"] for r in B])
SEASONS = ["Winter", "Spring", "Summer", "Autumn"]
sites = sorted(set(site))
NS = len(sites)
K = (stats.t.ppf(0.975, NS - 1) + stats.t.ppf(0.80, NS - 1)) / math.sqrt(NS)

print("\n=== univariate bounds (Table S9 rows; same method as the 12 measured covariates) ===")
print("%-10s %12s %10s %14s %10s" % ("window", "obs_perSD", "MDES", "obs_over_MDES", "P"))
for win in (1, 3, 7):
    vals = np.array([antecedent(r["season"], r["site"], win) for r in B])
    sd = float(np.std(vals, ddof=1))
    sl = []
    for s in sites:
        m = site == s
        if vals[m].std() == 0: continue
        sl.append(stats.linregress(vals[m], d[m]).slope)
    z = np.array(sl) * sd
    mdes = K * float(np.std(z, ddof=1)); obs = float(z.mean())
    p = stats.ttest_1samp(z, 0.0).pvalue
    print("%-10s %+12.4f %10.4f %14.3f %10.4f" % ("%d-day" % win, obs, mdes, abs(obs)/mdes, p))

print("\n=== season-adjusted slopes (Table S10 rain_*_pooled rows) ===")
S4 = np.column_stack([(sea == s).astype(float) for s in SEASONS])
for win in (1, 3, 7):
    v = np.array([antecedent(r["season"], r["site"], win) for r in B])
    X = np.column_stack([v, S4])
    b, *_ = np.linalg.lstsq(X, d, rcond=None); res = d - X @ b
    df = len(d) - X.shape[1]
    se = np.sqrt((res @ res / df) * np.diag(np.linalg.pinv(X.T @ X)))
    t = b[0] / se[0]
    print("   %d-day: b=%+.4f t=%+.3f df=%d P=%.3f" % (win, b[0], t, df, 2 * stats.t.sf(abs(t), df)))

def interaction(y):
    X = np.column_stack([np.ones(len(y))] + [(sea == s).astype(float) for s in SEASONS[1:]]
                        + [(site == s).astype(float) for s in sites[1:]])
    b, *_ = np.linalg.lstsq(X, y, rcond=None); r = y - X @ b; rss = r @ r
    dfres = len(y) - X.shape[1]
    X0 = np.column_stack([np.ones(len(y))] + [(site == s).astype(float) for s in sites[1:]])
    b0, *_ = np.linalg.lstsq(X0, y, rcond=None); r0 = y - X0 @ b0
    F = ((r0 @ r0 - rss) / 3) / (rss / dfres)
    return F, dfres, stats.f.sf(F, 3, dfres)

print("\n=== interaction after removing each rainfall term (Table S10) ===")
for win in (1, 3, 7):
    v = np.array([antecedent(r["season"], r["site"], win) for r in B])
    X = np.column_stack([np.ones(len(v)), v])
    b, *_ = np.linalg.lstsq(X, d, rcond=None)
    F, dfres, p = interaction(d - X @ b)
    print("   %d-day removed: F(3,%d)=%.3f  P=%.4f" % (win, dfres, F, p))

print("\n=== per-campaign mean antecedent rainfall (mm) ===")
for win in (1, 3, 7):
    out = []
    for s in SEASONS:
        vv = [antecedent(r["season"], r["site"], win) for r in B if r["season"] == s]
        out.append("%s %.1f" % (s[:2], float(np.mean(vv))))
    print("   %d-day: %s" % (win, "  ".join(out)))

print("\n=== cyclone lag within the summer campaign ===")
REF = dt.date(2023, 7, 29)   # heaviest-rainfall day of the Doksuri period
lag, con = [], []
for r in B:
    if r["season"] != "Summer": continue
    lag.append((dt.date.fromisoformat(dates[("Summer", r["site"])]) - REF).days)
    con.append(float(r["d"]))
lag = np.array(lag); con = np.array(con)
rr, pp = stats.pearsonr(lag, con)
U, pm = stats.mannwhitneyu(con[lag <= 6], con[lag > 6], alternative="two-sided")
print("   r = %+.3f  df = %d  P = %.3f" % (rr, len(lag) - 2, pp))
print("   within 6 d (n=%d) vs later (n=%d): U = %.1f  P = %.3f"
      % ((lag <= 6).sum(), (lag > 6).sum(), U, pm))
print("   territory-mean rainfall 29/30/31 Jul: %.1f / %.1f / %.1f mm"
      % (rain["2023-07-29"], rain["2023-07-30"], rain["2023-07-31"]))
print("   7-8 Sep: %.1f / %.1f mm" % (rain["2023-09-07"], rain["2023-09-08"]))
print("\n=== 7-day rainfall vs water temperature ===")
v7 = np.array([antecedent(r["season"], r["site"], 7) for r in B])
tp = np.array([float(r["temp"]) for r in B])
print("   r = %+.3f" % stats.pearsonr(v7, tp)[0])
