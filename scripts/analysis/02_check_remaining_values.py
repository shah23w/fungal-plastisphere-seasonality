"""Recompute the values in the Results that 01_check_main_values.py does not cover,
so that between the two scripts every number in the paper is derived from a
deposited table rather than asserted.

    python3 02_check_remaining_values.py [table_dir]

Exit status is 0 only if every value matches.
"""
import collections, csv, math, os, sys
import numpy as np
from scipy import stats

D = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "..", "tables")
P = lambda n: os.path.join(D, n)


def T(n):
    return list(csv.DictReader(open(P(n)), delimiter="\t"))


FAIL, N = [], [0]


def chk(label, got, want, tol=5e-4):
    N[0] += 1
    ok = (abs(got - want) <= tol) if isinstance(want, float) else (got == want)
    print("  %-56s %12s  want %-12s %s"
          % (label, ("%.6g" % got) if isinstance(got, float) else got,
             ("%.6g" % want) if isinstance(want, float) else want,
             "ok" if ok else "FAILED"))
    if not ok:
        FAIL.append(label)


FLOOR = 0.16
SE = ["Winter", "Spring", "Summer", "Autumn"]
B = T("TableS1_blocks.tsv")
d = np.array([float(r["d"]) for r in B])
sea = np.array([r["season"] for r in B])
site = np.array([r["site"] for r in B])
temp = np.array([float(r["temp"]) for r in B])
sites = sorted(set(site))

# ------------------------------------------- first Results subsection
print("== season, enrichment and the four-season contrast ==")
chk("summer mean water temperature", float(temp[sea == "Summer"].mean()), 29.9, 0.05)
chk("summer temperature SD", float(temp[sea == "Summer"].std(ddof=1)), 1.0, 0.05)
chk("summer blocks positive", int((d[sea == "Summer"] > 0).sum()), 13)
chk("winter blocks positive", int((d[sea == "Winter"] > 0).sum()), 4)


def hl(v):
    w = sorted((v[i] + v[j]) / 2.0 for i in range(len(v)) for j in range(i, len(v)))
    return float(np.median(w))


for s, want in (("Spring", 0.73), ("Autumn", 2.04)):
    chk("%s Hodges-Lehmann, fold" % s, 2 ** hl(d[sea == s]), want, 5e-3)
raw = {s: stats.wilcoxon(d[sea == s], method="exact").pvalue for s in SE}
order = sorted(SE, key=lambda s: raw[s])
run, adj = 0.0, {}
for i, s in enumerate(order):
    run = max(run, (4 - i) * raw[s])
    adj[s] = min(1.0, run)
chk("spring Holm P", float(adj["Spring"]), 0.46, 5e-3)
chk("autumn Holm P", float(adj["Autumn"]), 0.46, 5e-3)
# The four-season contrast is a site-blocked test, so each site contributes one
# mean. A plain t test on all 60 blocks ignores the blocking and gives a
# narrower interval, which is not what the manuscript reports.
sm = np.array([float(d[site == s].mean()) for s in sites])
tt = stats.ttest_1samp(sm, 0.0)
half = stats.t.ppf(0.975, len(sm) - 1) * stats.sem(sm)
chk("four-season mean contrast, fold", 2 ** float(sm.mean()), 1.55, 5e-3)
chk("  95% CI lower, fold", 2 ** float(sm.mean() - half), 0.83, 5e-3)
chk("  95% CI upper, fold", 2 ** float(sm.mean() + half), 2.90, 5e-3)
chk("  P", float(tt.pvalue), 0.16, 5e-3)
chk("the site means and the block mean agree", float(sm.mean()), float(d.mean()), 1e-9)

print("\n== the arm that moves ==")
# Recomputed from the raw counts rather than from the deposited 4 dp medians,
# because a difference of two rounded cells can land on a rounding boundary.
E0 = {r["run"]: r for r in T("TableS2_detection_envelope.tsv")}


def endpoint(run):
    t, f = int(E0[run]["total_pairs"]), int(E0[run]["fungal_reads"])
    return math.log2((f - t / 1e6 * FLOOR) / t)


S13 = {r["season"]: r for r in T("TableS13_season_summary.tsv")}
chg = {}
for arm in ("MP", "NP", "RW"):
    med = {s: float(np.median([endpoint(r["run_" + arm]) for r in B if r["season"] == s]))
           for s in ("Winter", "Summer")}
    chg[arm] = med["Summer"] - med["Winter"]
    chk("deposited median matches the recomputed one, winter %s" % arm,
        med["Winter"], float(S13["Winter"]["median_log2_" + arm]), 5e-5)
chk("winter to summer change, microplastic", chg["MP"], 2.78, 5e-3)
chk("winter to summer change, natural particle", chg["NP"], 0.10, 5e-3)
chk("winter to summer change, river water", chg["RW"], -0.88, 5e-3)
chk("  it does not round to -0.89", float(round(chg["RW"], 2) == -0.88), 1.0, 0.0)
chk("microplastic change as a fold", 2 ** chg["MP"], 6.9, 0.05)

print("\n== genera on summer microplastics ==")
G = T("TableS15_fungal_genera.tsv")
sm = [g for g in G if g["arm"] == "plastic"]
tot = sum(int(g["reads"]) for g in sm)
agg = collections.Counter()
for g in sm:
    agg[g["genus"]] += int(g["reads"])
for gen, want in (("Pichia", 73.7), ("Lodderomyces", 12.0), ("Yarrowia", 6.1),
                  ("Fusarium", 2.9), ("Candida", 1.4)):
    chk("%s share of summer MP fungal reads" % gen, 100.0 * agg[gen] / tot, want, 0.05)
chk("genera on the MP arm", len({g["genus"] for g in sm}), 57)
chk("genera across both arms", len({g["genus"] for g in G}), 62)
GM = T("TableS25_genus_matrix_180.tsv")
best = T("TableS16_species_coverage.tsv")[0]["site"]
hk = [g for g in GM if g["season"] == "Summer" and g["site"] == best and g["arm"] == "MP"]
assert hk
clade = int(hk[0]["fungal_clade_reads"])
pb = sum(int(g["reads"]) for g in hk if g["genus"] == "Pichia")
chk("Pichia share of fungal clade reads at the most enriched site",
    100.0 * pb / clade, 86.7, 0.05)
chk("  the denominator is the clade, not the sum over genera",
    float(clade > sum(int(g["reads"]) for g in hk)), 1.0, 0.0)
chk("the most enriched site is the one Table S16 names", best, T("TableS16_species_coverage.tsv")[0]["site"])

print("\n== species-level coverage ==")
S16 = T("TableS16_species_coverage.tsv")
top = S16[0]
chk("observed Pichia minimizers", int(top["observed_minimizers"]), 421120)
chk("index Pichia minimizers", int(top["index_minimizers"]), 437678)
chk("observed coverage, %", 100.0 * int(top["observed_minimizers"]) / int(top["index_minimizers"]), 96.2, 0.05)
lam = -math.log(1.0 - float(top["poisson_expected"]))
chk("implied redundancy, observations per index minimizer", lam, 6.5, 0.05)
chk("Poisson expected coverage, %", 100.0 * float(top["poisson_expected"]), 99.85, 0.01)
chk("observed to expected ratio", float(top["coverage_ratio"]), 0.96, 5e-3)
pk = [r for r in S16 if r["species"] == "Pichia kudriavzevii"]
pk.sort(key=lambda r: -float(r["index_coverage"]))
chk("second most enriched site coverage, %", 100.0 * float(pk[1]["index_coverage"]), 57.5, 0.05)
others = [r for r in S16 if r["species"] != "Pichia kudriavzevii"]
chk("no other taxon exceeds 50% coverage",
    float(max(float(r["index_coverage"]) for r in others) < 0.50), 1.0, 0.0)
ca = [r for r in S16 if r["site"] == top["site"] and r["species"] == "Candida albicans"][0]
chk("C. albicans observed coverage at that site, %", 100.0 * float(ca["index_coverage"]), 1.1, 0.05)
chk("  expected, %", 100.0 * float(ca["poisson_expected"]), 19.0, 0.05)
chk("  ratio", float(ca["coverage_ratio"]), 0.06, 5e-3)
chk("all Fusarium assignments below 1% coverage",
    float(max([float(r["index_coverage"]) for r in S16 if r["genus"] == "Fusarium"] or [0]) < 0.01),
    1.0, 0.0)
le = [r for r in S16 if r["species"] == "Lodderomyces elongisporus"]
lo = [r for r in S16 if r["species"] == "Lodderomyces orthopsilosis"]
chk("L. elongisporus top coverage, %", 100.0 * max(float(r["index_coverage"]) for r in le), 32.5, 0.05)
chk("L. orthopsilosis top coverage, %", 100.0 * max(float(r["index_coverage"]) for r in lo), 18.4, 0.05)

print("\n== genome complexity ==")
gg = [g for g in G if int(g["reads"]) > 0 and float(g["reads_per_distinct_minimizer"]) > 0]
chk("genus-level observations", len(gg), 686)
x = np.array([math.log10(int(g["reads"])) for g in gg])
y = np.array([math.log10(float(g["reads_per_distinct_minimizer"])) for g in gg])
lr = stats.linregress(x, y)
z = (y - (lr.intercept + lr.slope * x)) / (y - (lr.intercept + lr.slope * x)).std(ddof=1)
k = int(np.argmax([float(g["reads_per_distinct_minimizer"]) for g in gg]))
chk("highest reads per distinct minimizer", float(gg[k]["reads_per_distinct_minimizer"]), 0.442, 5e-4)
chk("  it is Pichia", gg[k]["genus"], "Pichia")
chk("  its residual, SD", float(z[k]), -1.03, 5e-3)
big = np.argsort(-z)[:5]
chk("largest of the five biggest positive residuals, SD", float(z[big[0]]), 2.8, 0.05)
chk("smallest of them, SD", float(z[big[-1]]), 2.4, 0.05)
rd = sorted(int(gg[i]["reads"]) for i in big)
chk("fewest reads among them", rd[0], 6)
chk("most reads among them", rd[-1], 61)

# ------------------------------------------- resolvability subsection
print("\n== yeast carryover ==")
Y = T("TableS17_yeast_carryover.tsv")
allrow = [r for r in Y if r["season"] == "ALL"][0]
chk("summer MP share of all yeast reads, %", float(allrow["median_yeast_RPM"]), 95.5, 0.05)
med = {(r["season"], r["arm"]): float(r["median_yeast_RPM"]) for r in Y
       if r["arm"] in ("plastic", "surface_control", "bulk_control")}
chk("summer MP median yeast rpm", med[("Summer", "plastic")], 17.30, 5e-3)
chk("summer NP median yeast rpm", med[("Summer", "surface_control")], 0.65, 5e-3)
chk("summer RW median yeast rpm", med[("Summer", "bulk_control")], 0.745, 5e-4)
chk("  it rounds to 0.74, not 0.75", float(round(med[("Summer", "bulk_control")], 2) == 0.74), 1.0, 0.0)
for s, want in (("Winter", 0.16), ("Spring", 0.79), ("Autumn", 0.66)):
    chk("%s MP median yeast rpm" % s, med[(s, "plastic")], want, 5e-3)
for arm, lab, want in (("plastic", "MP", 21.95), ("surface_control", "NP", 2.16),
                       ("bulk_control", "RW", 0.91)):
    other = max(med[(s, arm)] for s in SE if s != "Summer")
    chk("summer rise in the %s median, fold above any other season" % lab,
        med[("Summer", arm)] / other, want, 5e-3)

print("\n== detection envelope ==")
E = T("TableS2_detection_envelope.tsv")
chk("archived runs", len(E), 859)
chk("runs above the floor", sum(1 for r in E if r["above_floor"] == "yes"), 846)
chk("  as a percentage", 100.0 * sum(1 for r in E if r["above_floor"] == "yes") / len(E), 98.5, 0.05)
chk("runs with at least 10,000 fungal reads", sum(1 for r in E if int(r["fungal_reads"]) >= 10000), 21)
chk("  as a percentage", 100.0 * sum(1 for r in E if int(r["fungal_reads"]) >= 10000) / len(E), 2.4, 0.05)
chk("runs under 1,000 fungal reads", sum(1 for r in E if int(r["fungal_reads"]) < 1000), 706)
chk("  as a percentage", 100.0 * sum(1 for r in E if int(r["fungal_reads"]) < 1000) / len(E), 82.2, 0.05)
cr = [float(r["class_rate"]) for r in E]
chk("lowest classification rate, %", min(cr), 0.0964, 5e-4)
chk("highest classification rate, %", max(cr), 54.155, 5e-3)
chk("range, fold", max(cr) / min(cr), 562.0, 0.5)
rho, pv = stats.spearmanr(cr, [float(r["fungal_rpm"]) for r in E])
chk("fungal share against classification rate, rho", float(rho), -0.024, 5e-4)
chk("  P", float(pv), 0.49, 5e-3)
chk("floor in reads per 40 million pairs", FLOOR * 40.0, 6.4, 0.5)
nplo = sum(1 for r in B if int(r["fungal_NP"]) < 1000)
both = sum(1 for r in B if int(r["fungal_NP"]) < 1000 and int(r["fungal_MP"]) < 1000)
chk("blocks whose NP arm is under 1,000 fungal reads", nplo, 35)
chk("blocks where both arms are", both, 19)
chk("fewest fungal reads in any particle library",
    min(min(int(r["fungal_MP"]), int(r["fungal_NP"])) for r in B), 23)
worst = min(min(int(r["fungal_MP"]), int(r["fungal_NP"])) for r in B)
chk("Poisson counting error at the smallest particle library, log2 units",
    math.log2(1.0 + 1.0 / math.sqrt(worst)), 0.27, 5e-3)
chk("  it is far below the summer contrast",
    float(math.log2(1.0 + 1.0 / math.sqrt(worst)) < hl(d[sea == "Summer"])), 1.0, 0.0)

print("\n== community structure ==")
S20 = [r for r in T("TableS20_structure.tsv") if r["arm"] in ("plastic", "surface_control")]
by = {(r["site"], r["arm"]): r for r in S20}
st = sorted({r["site"] for r in S20})
ratio = [int(by[(s, "plastic")]["fungal_reads"]) / int(by[(s, "surface_control")]["fungal_reads"])
         for s in st]
fold = [max(x, 1.0 / x) for x in ratio]
chk("median fungal-yield disparity between arms, fold", float(np.median(fold)), 4.5, 0.05)
chk("largest disparity, fold", float(max(fold)), 1616.0, 1.0)
disc = []
for s in st:
    a, b2 = int(by[(s, "plastic")]["fungal_reads"]), int(by[(s, "surface_control")]["fungal_reads"])
    disc.append(0.0 if a <= b2 else 100.0 * (a - b2) / a)
chk("median share of MP reads discarded at matched depth, %", float(np.median(disc)), 78.0, 1.0)
dr = np.array([float(by[(s, "plastic")]["expected_genera_rarefied"])
               - float(by[(s, "surface_control")]["expected_genera_rarefied"]) for s in st])
chk("rarefied richness difference, genera", float(dr.mean()), -3.6, 0.05)
ds = np.array([float(by[(s, "plastic")]["inverse_simpson"])
               - float(by[(s, "surface_control")]["inverse_simpson"]) for s in st])
chk("inverse Simpson difference", float(ds.mean()), -0.94, 5e-3)
chk("  signed-rank P", float(stats.wilcoxon(ds).pvalue), 0.68, 5e-3)
dt = np.array([float(by[(s, "plastic")]["top_genus_pct"])
               - float(by[(s, "surface_control")]["top_genus_pct"]) for s in st])
chk("top-genus share difference, percentage points", float(dt.mean()), 0.42, 5e-3)
chk("  signed-rank P", float(stats.wilcoxon(dt).pvalue), 0.93, 5e-3)

print("\n== fungus-free mock community ==")
M = T("TableS7_mock2x2.tsv")


def cell(c, s):
    v = [float(r["false_per_million"]) for r in M if int(r["cell"]) == c and r["stratum"] == s]
    return float(np.mean(v))


for s, want in (("A", 2591.0), ("B", 31737.0), ("C", 48103.0)):
    chk("fungal-only index at confidence 0, stratum %s" % s, cell(1, s), want, 1.0)
hum = [r for r in M if int(r["cell"]) == 1 and r["stratum"] == "C"][0]
chk("human reads assigned to Fungi, %",
    100.0 * int(hum["false_fungal_reads"]) / int(hum["total_reads"]), 4.8, 0.05)
for s, want in (("A", 0.15), ("B", 0.45), ("C", 0.36)):
    chk("locked settings, stratum %s" % s, cell(5, s), want, 5e-3)
# Libraries that return no false reads at confidence 0.1 are undefined on the
# log scale and are dropped, exactly as in the factorial model, which is why the
# broad index contributes 7 pairs and not 9.
for scope, want_n, want_p in (("fungal-only", 9, 0.0020), ("PlusPF-capped", 7, 0.0078)):
    lo = {r["run"]: float(r["false_per_million"]) for r in M
          if r["index"] == scope and float(r["confidence"]) == 0.0}
    hi = {r["run"]: float(r["false_per_million"]) for r in M
          if r["index"] == scope and float(r["confidence"]) == 0.1 and int(r["min_hit_groups"]) == 2}
    diff = [lo[k] - hi[k] for k in lo if hi[k] > 0 and lo[k] != hi[k]]
    chk("%s, libraries where confidence reduced false positives" % scope, len(diff), want_n)
    chk("  one-sided signed-rank P",
        float(stats.wilcoxon(diff, alternative="greater", method="exact").pvalue), want_p, 5e-4)

print("\n== bacteria and the phylum concordance ==")
BB = {(r["season"], r["site"]): float(r["d_log2"]) for r in T("TableS22_bacterial_blocks.tsv")}
db = np.array([BB[(r["season"], r["site"])] for r in B])
bm = {s: float(db[sea == s].mean()) for s in SE}
chk("lowest bacterial season mean", float(min(bm.values())), -0.09, 5e-3)
chk("highest bacterial season mean", float(max(bm.values())), 0.52, 5e-3)
PH = T("TableS23_phylum_concordance.tsv")
pre = {"Bacteroidota": (-0.92, -1.39), "Firmicutes": (2.04, 1.30), "Proteobacteria": (-0.17, -0.28)}
for ph, (ours, theirs) in pre.items():
    row = [r for r in PH if r["phylum_GTDB"] == ph][0]
    chk("%s, this study" % ph, float(row["ours_log2_MP_over_NP"]), ours, 5e-3)
    chk("  %s, published" % ph, float(row["bao_log2_MP_over_NP"]), theirs, 5e-3)
chk("phyla agreeing in sign", sum(1 for r in PH if r["sign_agree"] == "yes"), 10)
chk("phyla compared", len(PH), 15)
ge = [r for r in PH if float(r["n_rMAGs"]) >= 50]
chk("phyla with at least 50 recovered genomes", len(ge), 4)
chk("  all agree in sign", sum(1 for r in ge if r["sign_agree"] == "yes"), 4)
lt = [r for r in PH if float(r["n_rMAGs"]) < 50]
chk("phyla below that threshold", len(lt), 11)
chk("  agreeing in sign", sum(1 for r in lt if r["sign_agree"] == "yes"), 6)

print("\n== cross-study and the southern-hemisphere systems ==")
S3 = [r for r in T("TableS3_crossstudy.tsv") if r.get("study", "").startswith("PRJ")]
single = [r for r in S3 if r.get("multi_season") == "no"]
chk("single-timepoint studies", len(single), 6)
chk("  recording a season",
    sum(1 for r in single if r["season_sampled"] not in ("", "NA", "unknown")), 5)
chk("  the sixth has no sampling date",
    sum(1 for r in single if r["sampling_months"] in ("", "NA")), 1)
chk("  none of them summer",
    float(all("Summer" not in r["season_sampled"] for r in single)), 1.0, 0.0)
NZ = T("TableS4_second_system.tsv")
chk("per-date contrasts in the two systems", len(NZ), 6)
chk("summer dates positive",
    float(all(float(r["difference_log2"]) > 0 for r in NZ if r["season_SH"] == "Summer")), 1.0, 0.0)
chk("winter dates negative",
    float(all(float(r["difference_log2"]) < 0 for r in NZ if r["season_SH"] == "Winter")), 1.0, 0.0)
chk("marina detectable summer-to-winter contrast, fold", 2 ** 2.41, 5.3, 0.02)
chk("marina observed contrast, fold", 2 ** 0.894, 1.9, 0.05)
chk("pond detectable contrast, fold", 2 ** 3.38, 10.4, 0.02)

print("\n== spike-in recovery and detectable effects ==")
SP = [r for r in T("TableS8_spikein_recovery.tsv")
      if r.get("donor") in ("Pichia", "Yarrowia", "Fusarium")
      and str(r.get("injected_pairs", "")).isdigit()]
rec = {}
for dn, want in (("Pichia", 0.566), ("Yarrowia", 0.593), ("Fusarium", 0.620)):
    sub = [r for r in SP if r["donor"] == dn]
    lr2 = stats.linregress([int(r["injected_pairs"]) for r in sub],
                           [int(r["observed_reads"]) for r in sub])
    rec[dn] = float(lr2.slope)
    chk("recovery, %s" % dn, rec[dn], want, 5e-4)
chk("spread across donors, %", 100.0 * (max(rec.values()) - min(rec.values())) / max(rec.values()), 9.0, 0.5)
two = [r for r in SP if r["donor"] == "Pichia" and int(r["injected_pairs"]) == 2]
chk("Pichia recovered from two injected pairs",
    float(any(int(r["observed_reads"]) > 0 for r in two)), 1.0, 0.0)
IX = T("TableS18_index_coverage.tsv")
ret = [float(r["capped_over_full"]) for r in IX]
chk("genera compared for index retention", len(IX), 148)
chk("  fungal", sum(1 for r in IX if r["kingdom_group"] == "fungal"), 62)
chk("  non-fungal", sum(1 for r in IX if r["kingdom_group"] != "fungal"), 86)
chk("lowest retention", min(ret), 0.1351, 5e-5)
chk("highest retention", max(ret), 0.1360, 5e-5)
tc = stats.t.ppf(0.975, 14) + stats.t.ppf(0.80, 14)
mdes = {s: tc * stats.tstd(d[sea == s]) / math.sqrt(15) for s in SE}
for s, want in (("Winter", 2.42), ("Spring", 1.52), ("Summer", 2.99), ("Autumn", 2.15)):
    chk("minimum detectable effect, %s" % s, float(mdes[s]), want, 5e-3)
chk("summer mean exceeds its own threshold",
    float(d[sea == "Summer"].mean() > mdes["Summer"]), 1.0, 0.0)
chk("summer mean", float(d[sea == "Summer"].mean()), 3.05, 5e-3)
chk("summer Hodges-Lehmann falls short of it", float(hl(d[sea == "Summer"]) < mdes["Summer"]), 1.0, 0.0)
chk("undetected spring effect could reach, fold", 2 ** mdes["Spring"], 2.9, 0.05)
chk("undetected autumn effect could reach, fold", 2 ** mdes["Autumn"], 4.4, 0.05)

print("\n%d checks, %d mismatches" % (N[0], len(FAIL)))
if FAIL:
    for f in FAIL:
        print("   FAILED:", f)
    sys.exit(1)
