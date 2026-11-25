"""Within-site temperature slopes for the six yeast genera, on microplastics and on
the microplastic minus natural particle contrast, with Holm adjustment across the
12 tests. Uses the paper's endpoint, floor-corrected reads over total read pairs,
with the pseudocount at one read, which is the single-read detection limit of
these libraries. Also reports the association between the slope ordering and
documented genus-level maximum growth temperatures.

    python3 06_genus_slopes.py [table_dir]

Writes TableS24_thermotolerance.tsv.
"""
import csv, math, os, sys, collections
import numpy as np
from scipy import stats

D = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
def path(n):
    for c in (os.path.join(D, n), os.path.join(D, "..", "..", "tables", n), n):
        if os.path.exists(c):
            return c
    raise SystemExit("cannot find %s (pass the tables directory as argv[1])" % n)

FLOOR = 0.16
GEN = ["Pichia", "Lodderomyces", "Yarrowia", "Candida", "Torulaspora", "Debaryomyces"]
# Documented genus-level maximum growth temperatures (deg C). Approximations:
# a genus-level value is not meaningful for a polyphyletic genus such as Candida.
TMAX = {"Pichia": 45.0, "Yarrowia": 32.0, "Candida": 45.0,
        "Lodderomyces": 37.0, "Torulaspora": 35.0, "Debaryomyces": 35.0}

B = list(csv.DictReader(open(path("TableS1_blocks.tsv")), delimiter="\t"))
G = list(csv.DictReader(open(path("TableS25_genus_matrix_180.tsv")), delimiter="\t"))
TEMP = {(r["season"], r["site"]): float(r["temp"]) for r in B}
tot = {}
for r in B:
    tot[r["run_MP"]] = int(r["total_MP"])
    tot[r["run_NP"]] = int(r["total_NP"])
reads = collections.defaultdict(dict)
meta = {}
for r in G:
    reads[r["run"]][r["genus"]] = int(r["reads"])
    meta[r["run"]] = (r["season"], r["site"], r["arm"])
run_of = {v: k for k, v in meta.items()}
sites = sorted({m[1] for m in meta.values()})
SEASONS = ["Winter", "Spring", "Summer", "Autumn"]
print("blocks %d | genus rows %d | sites %d" % (len(B), len(G), len(sites)))


def endpoint(run, g):
    """Floor-corrected log2 read share for one genus in one library.
    Pseudocount = 1 read, the single-read detection limit."""
    T = tot[run]
    f = reads[run].get(g, 0) - (T / 1e6) * FLOOR
    return math.log2(max(f, 1.0) / T)


def slopes(g, arms):
    out = []
    for site in sites:
        xs, ys = [], []
        for se in SEASONS:
            try:
                v = {a: endpoint(run_of[(se, site, a)], g) for a in arms}
            except KeyError:
                continue
            ys.append(v["MP"] - v["NP"] if len(arms) == 2 else v["MP"])
            xs.append(TEMP[(se, site)])
        if len(xs) >= 3:
            out.append(stats.linregress(xs, ys).slope)
    return np.array(out)


rows = []
for g in GEN:
    rec = {"genus": g, "documented_Tmax_C": TMAX[g]}
    for lab, arms in (("plastic", ["MP"]), ("difference", ["MP", "NP"])):
        a = slopes(g, arms)
        t = stats.ttest_1samp(a, 0.0)
        rec[lab + "_slope"] = round(float(a.mean()), 4)
        rec[lab + "_t"] = round(float(t.statistic), 3)
        rec[lab + "_P"] = float(t.pvalue)
        rec["n_sites"] = len(a)
    rows.append(rec)

# Holm across the 12 tests
keys = [(r["genus"], lab) for r in rows for lab in ("plastic", "difference")]
pv = sorted(((k, dict((r["genus"], r) for r in rows)[k[0]][k[1] + "_P"]) for k in keys),
            key=lambda x: x[1])
adj, run = {}, 0.0
for i, (k, p) in enumerate(pv):
    run = max(run, (len(pv) - i) * p)
    adj[k] = min(1.0, run)
for r in rows:
    for lab in ("plastic", "difference"):
        r[lab + "_P_holm"] = round(adj[(r["genus"], lab)], 4)
        r[lab + "_P"] = float("%.3g" % r[lab + "_P"])

cols = ["genus", "documented_Tmax_C", "n_sites",
        "plastic_slope", "plastic_t", "plastic_P", "plastic_P_holm",
        "difference_slope", "difference_t", "difference_P", "difference_P_holm"]
OUT = path("TableS1_blocks.tsv").replace("TableS1_blocks.tsv",
                                         "TableS24_thermotolerance.tsv")
with open(OUT, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
    w.writeheader()
    for r in rows:
        w.writerow({c: r[c] for c in cols})
print("wrote", OUT)
print("")
print("%-14s %6s %10s %10s %10s   %12s %10s %10s"
      % ("genus", "Tmax", "slope(MP)", "P", "Holm", "slope(MP-NP)", "P", "Holm"))
for r in rows:
    print("%-14s %6.0f %10.4f %10.2e %10.4f   %12.4f %10.2e %10.4f"
          % (r["genus"], r["documented_Tmax_C"], r["plastic_slope"], r["plastic_P"],
             r["plastic_P_holm"], r["difference_slope"], r["difference_P"],
             r["difference_P_holm"]))

# the thermotolerance test the manuscript now reports
print("")
print("=== ordering vs documented maximum growth temperature ===")
for lab in ("plastic", "difference"):
    x = [r["documented_Tmax_C"] for r in rows]
    y = [r[lab + "_slope"] for r in rows]
    rr, pp = stats.pearsonr(x, y)
    print("   %-11s r = %+.3f   P = %.3f   n = %d genera" % (lab, rr, pp, len(x)))
print("   Six genera resolve only very large effects; a genus-level Tmax is an")
print("   approximation and is not meaningful for polyphyletic Candida.")
