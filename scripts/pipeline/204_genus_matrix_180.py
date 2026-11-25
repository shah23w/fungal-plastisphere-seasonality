"""Build Table S25, the per-library fungal genus read counts for all 180 libraries,
from the Kraken2 reports, and recompute the season by substrate interaction from
the file it writes, with and without the Pichia clade.
"""
import csv, os
import numpy as np
from scipy import stats

R1 = os.path.expanduser("~/fpar/results/134_reconcile_PRJNA1159715.tsv")
R2 = os.path.expanduser("~/fpar/results/134_reconcile_PRJNA1330252.tsv")
KD = os.path.expanduser("~/fpar/results/kreports")
F1 = os.path.expanduser("~/fpar/figures/Figure1/Figure1_source_data.tsv")
OUT = os.path.expanduser("~/fpar/FINAL/tables")
os.makedirs(OUT, exist_ok=True)
FLOOR = 0.16

ARMNAME = {"plastic": "MP", "surface_control": "NP", "bulk_control": "RW"}
meta = {}
for r in csv.DictReader(open(R1), delimiter="\t"):
    meta[r["run"]] = ("Winter", r["site"], r["archive_arm"])
for r in csv.DictReader(open(R2), delimiter="\t"):
    s = r["site"]
    se, si = (s.split("-", 1) + [""])[:2]
    meta[r["run"]] = (se, si, r["archive_arm"])
print("runs in reconcile files: %d" % len(meta))


def fungal_genera(path):
    """(fungal_clade_reads, {genus: clade_reads}) from one kreport."""
    if not os.path.exists(path):
        return None, None
    rows = []
    for ln in open(path):
        p = ln.rstrip("\n").split("\t")
        if len(p) < 8:
            continue
        nm = p[7]
        ind = len(nm) - len(nm.lstrip())
        try:
            rows.append((int(p[1]), p[5], nm.strip(), ind))
        except ValueError:
            continue
    inside = False
    base = None
    ftot = 0
    gen = {}
    for cl, rk, nm, ind in rows:
        if nm == "Fungi" and rk == "K" and not inside:
            inside = True
            base = ind
            ftot = cl
            continue
        if inside:
            if ind <= base:
                break
            if rk == "G" and cl > 0:
                gen[nm] = gen.get(nm, 0) + cl
    return (ftot, gen) if inside else (None, None)


# ---------------------------------------------------------------- write matrix
tot = {}
for r in csv.DictReader(open(F1), delimiter="\t"):
    for a, rk, tk in (("MP", "run_MP", "total_MP"), ("NP", "run_NP", "total_NP")):
        tot[r[rk]] = int(r[tk])

rows = []
missing = []
for run, (se, si, arm) in sorted(meta.items()):
    ftot, gen = fungal_genera(os.path.join(KD, "%s.kreport" % run))
    if ftot is None:
        missing.append(run)
        continue
    for g, c in sorted(gen.items(), key=lambda kv: -kv[1]):
        rows.append(dict(season=se, site=si, arm=ARMNAME.get(arm, arm), run=run,
                         total_read_pairs=tot.get(run, ""), fungal_clade_reads=ftot,
                         genus=g, reads=c))
if missing:
    print("   kreports MISSING for %d runs: %s" % (len(missing), missing[:8]))
P = os.path.join(OUT, "TableS25_genus_matrix_180.tsv")
cols = ["season", "site", "arm", "run", "total_read_pairs",
        "fungal_clade_reads", "genus", "reads"]
with open(P, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
libs = {r["run"] for r in rows}
print("wrote %s" % P)
print("   %d rows, %d libraries, %d distinct genera"
      % (len(rows), len(libs), len({r["genus"] for r in rows})))

# ------------------------------------------- reproduce the published interaction
src = list(csv.DictReader(open(F1), delimiter="\t"))
byrun = {}
for r in rows:
    byrun.setdefault(r["run"], {})[r["genus"]] = int(r["reads"])
ftot_by = {r["run"]: int(r["fungal_clade_reads"]) for r in rows}


def endpoint(run, drop=None):
    T = tot[run]
    F = ftot_by.get(run, 0)
    if drop:
        F -= byrun.get(run, {}).get(drop, 0)
    v = F - (T / 1e6) * FLOOR
    return np.log2(v / T) if v > 0 else np.nan


def interaction(drop):
    y, sea, site = [], [], []
    for r in src:
        a, b = endpoint(r["run_MP"], drop), endpoint(r["run_NP"], drop)
        if np.isnan(a) or np.isnan(b):
            continue
        y.append(a - b); sea.append(r["season"]); site.append(r["site"])
    y = np.array(y); sea = np.array(sea); site = np.array(site)
    SE = ["Winter", "Spring", "Summer", "Autumn"]; SI = sorted(set(site))
    X = np.column_stack([np.ones(len(y))]
                        + [(sea == s).astype(float) for s in SE[1:]]
                        + [(site == s).astype(float) for s in SI[1:]])
    b_, *_ = np.linalg.lstsq(X, y, rcond=None); r_ = y - X @ b_; rss = r_ @ r_
    dfres = len(y) - X.shape[1]
    X0 = np.column_stack([np.ones(len(y))]
                         + [(site == s).astype(float) for s in SI[1:]])
    b0, *_ = np.linalg.lstsq(X0, y, rcond=None); r0 = y - X0 @ b0; rss0 = r0 @ r0
    F = ((rss0 - rss) / 3) / (rss / dfres)
    mu = {s: float(np.mean(y[sea == s])) for s in SE}
    return F, dfres, stats.f.sf(F, 3, dfres), len(y), mu


print("")
print("=== VERIFICATION: recomputed from the file just written ===")
for lab, drop, eF, eP in (("all fungal reads", None, 6.064, 0.0016),
                          ("Pichia reads removed", "Pichia", 4.790, 0.0058)):
    F, dfres, p, n, mu = interaction(drop)
    ok = "MATCH" if abs(F - eF) < 0.02 else "*** MISMATCH ***"
    print("  %-22s F(3,%d) = %.3f  P = %.4f   n = %d blocks   [published %.3f / %.4f]  %s"
          % (lab, dfres, F, p, n, eF, eP, ok))
    print("     season means: " + "  ".join("%s %+.3f" % (k, v) for k, v in mu.items()))
print("")
print("  summer mean should fall from 3.05 to 2.40 log2 units on Pichia removal")
