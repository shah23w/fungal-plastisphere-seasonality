"""Bound the technical contribution to the winter result using the bacterial
contrast computed on the same 60 blocks.

All 15 winter blocks come from PRJNA1159715 and all 45 others from PRJNA1330252,
so for winter alone season is confounded with submission, and therefore with
library preparation and sequencing batch. Batch and season are perfectly
collinear here, so ComBat, RUV and related methods cannot separate them.

Anything technical that displaced the microplastic libraries relative to the
natural-particle libraries within a submission displaces both kingdoms. This
script quantifies that shared component and refits the interaction with the
bacterial contrast entered as a covariate. It bounds artifacts shared by the two
kingdoms, not a fungus-specific contaminant, which these data cannot exclude.

    python3 08_bacterial_contrast.py [table_dir]

Exit status is 0 only if every check passes.
"""
import csv, math, os, sys
import numpy as np
from scipy import stats

D = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "tables")
P = lambda n: os.path.join(D, n)
rd = lambda n: list(csv.DictReader(open(P(n)), delimiter="\t"))
FAIL, NCHK = [], [0]
def chk(label, got, want, tol):
    NCHK[0] += 1
    ok = abs(float(got) - float(want)) <= tol
    if not ok: FAIL.append(label)
    print("  %-56s %12.5f  expected %-11s %s" % (label[:56], float(got), "%g" % want, "ok" if ok else "MISMATCH"))

S1, BK, S13 = rd("TableS1_blocks.tsv"), rd("TableS22_bacterial_blocks.tsv"), rd("TableS12_interaction_tests.tsv")
S13 = {r["kingdom"]: r for r in S13}
key = lambda r: (r["season"], r["site"])
bd = {key(r): float(r["d_log2"]) for r in BK}
tb = {key(r): float(r["temp_C"]) for r in BK}
SE = ["Winter", "Spring", "Summer", "Autumn"]
f = np.array([float(r["d"]) for r in S1])
b = np.array([bd[key(r)] for r in S1])
sea = np.array([r["season"] for r in S1]); site = np.array([r["site"] for r in S1])
sites = sorted(set(site))

print("=== the bacterial table is the same 60 blocks ===")
chk("blocks", len(BK), 60, 0)
chk("every block matched", float({key(r) for r in BK} == {key(r) for r in S1}), 1, 0)
chk("temperatures agree with Table S1", float(all(abs(tb[key(r)] - float(r["temp"])) < 1e-9 for r in S1)), 1, 0)
for s in SE:
    chk("%s bacterial mean reproduces Table S12" % s, b[sea == s].mean(), float(S13["Bacteria"]["mean_d_%s" % s]), 5e-4)

def interaction(y, cov=None):
    S4 = np.column_stack([(sea == s).astype(float) for s in SE])
    SI = np.column_stack([(site == s).astype(float) for s in sites])
    base = [np.ones(60)] + ([cov] if cov is not None else [])
    X = np.column_stack(base + [S4[:, 1:], SI[:, 1:]])
    bb, *_ = np.linalg.lstsq(X, y, rcond=None); r = y - X @ bb
    X0 = np.column_stack(base + [SI[:, 1:]])
    b0, *_ = np.linalg.lstsq(X0, y, rcond=None); r0 = y - X0 @ b0
    df = 60 - np.linalg.matrix_rank(X)
    F = ((r0 @ r0 - r @ r) / 3) / ((r @ r) / df)
    return F, df, stats.f.sf(F, 3, df)

print("\n=== controls ===")
F, df, p = interaction(b)
chk("bacterial interaction F", F, float(S13["Bacteria"]["F"]), 5e-3)
chk("  df", df, 42, 0)
chk("  P", p, float(S13["Bacteria"]["P"]), 5e-5)
F, df, p = interaction(f)
chk("fungal interaction F", F, float(S13["Fungi"]["F"]), 5e-3)
chk("  P", p, float(S13["Fungi"]["P"]), 5e-5)

print("\n=== direction: the two kingdoms move opposite ways in winter ===")
chk("fungal winter mean", f[sea == "Winter"].mean(), -1.2215, 5e-4)
chk("  it is the lowest of the four seasons",
    float(f[sea == "Winter"].mean() == min(f[sea == s].mean() for s in SE)), 1, 0)
chk("bacterial winter mean", b[sea == "Winter"].mean(), 0.5237, 5e-4)
chk("  it is the highest of the four seasons",
    float(b[sea == "Winter"].mean() == max(b[sea == s].mean() for s in SE)), 1, 0)
tt = stats.ttest_ind(b[sea == "Winter"], b[sea != "Winter"])
chk("  winter vs the rest, bacterial P", tt.pvalue, 0.163, 5e-3)

print("\n=== what a shared winter offset would have had to be to escape notice ===")
K = (stats.t.ppf(0.975, 14) + stats.t.ppf(0.80, 14)) / math.sqrt(15)
chk("detectable bacterial winter offset at 80% power", K * float(np.std(b[sea == "Winter"], ddof=1)), 0.98, 5e-3)
chk("  observed bacterial winter offset", b[sea == "Winter"].mean(), 0.52, 5e-3)
chk("  fungal winter contrast to be explained", f[sea == "Winter"].mean(), -1.22, 5e-3)

print("\n=== adjusting the fungal contrast for the bacterial contrast ===")
r = stats.pearsonr(f[sea == "Winter"], b[sea == "Winter"])
chk("within winter, fungal vs bacterial r", r.statistic, 0.714, 5e-3)
chk("  P", r.pvalue, 0.003, 5e-4)
r = stats.pearsonr(f, b)
chk("across all 60 blocks, r", r.statistic, 0.240, 5e-3)
chk("  P", r.pvalue, 0.065, 5e-3)
F, df, p = interaction(f, cov=b)
chk("interaction with the bacterial contrast entered, F", F, 7.474, 5e-3)
chk("  df", df, 41, 0)
chk("  P", p, 0.0004, 5e-5)
chk("  it does not weaken the interaction", float(p < float(S13["Fungi"]["P"])), 1, 0)
X = np.column_stack([np.ones(60), b] + [(site == s).astype(float) for s in sites[1:]])
bb, *_ = np.linalg.lstsq(X, f, rcond=None)
chk("slope of fungal on bacterial contrast", bb[1], 0.712, 5e-3)
adj = f - bb[1] * (b - b.mean())
chk("winter fungal mean after removing that component", adj[sea == "Winter"].mean(), -1.476, 5e-3)
chk("  it becomes more negative, not less", float(adj[sea == "Winter"].mean() < f[sea == "Winter"].mean()), 1, 0)

print("\n%d checks, %d mismatches" % (NCHK[0], len(FAIL)))
if FAIL:
    for x in FAIL: print("   MISMATCH:", x)
    sys.exit(1)
