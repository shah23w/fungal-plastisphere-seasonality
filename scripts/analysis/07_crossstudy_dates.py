"""Join sampling months and the inferred season into Table S3 and check the
statement that five of the six single-timepoint studies sampled outside summer,
the sixth recording no sampling date.

    python3 07_crossstudy_dates.py [table_dir]

Rewrites TableS3_crossstudy.tsv. Re-running is safe. Exit status is 0 only if
every check passes.
"""
import csv, os, sys

D = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "tables")
P = lambda n: os.path.join(D, n)
FAIL, NCHK = [], [0]
def chk(label, ok):
    NCHK[0] += 1
    if not ok: FAIL.append(label)
    print("  %-64s %s" % (label[:64], "ok" if ok else "FAILED"))

with open(P("TableS3_crossstudy.tsv")) as fh:
    rd = csv.DictReader(fh, delimiter="\t"); COLS, S3 = rd.fieldnames, list(rd)
CS = {r["study"]: r for r in csv.DictReader(open(P("crossstudy_season.tsv")), delimiter="\t")}

chk("Table S3 holds 10 studies", len(S3) == 10)
chk("every study has a season record", all(r["study"] in CS for r in S3))
chk("no extra studies in the season record", set(CS) == {r["study"] for r in S3})
already = "sampling_months" in COLS
if already:
    print("  columns already present, verifying only")
else:
    COLS = COLS + ["sampling_months", "season_sampled", "multi_season"]
for r in S3:
    c = CS[r["study"]]
    r["sampling_months"] = c["months"]
    r["season_sampled"] = c["season"]
    r["multi_season"] = c["multi_season"]

single = [r for r in S3 if r["multi_season"] == "no"]
chk("six studies sampled at a single timepoint", len(single) == 6)
known = [r for r in single if r["season_sampled"] != "unknown"]
chk("  five of them record a season", len(known) == 5)
chk("  none of those five is summer", not any("Summer" in r["season_sampled"] for r in known))
chk("  the sixth records no month", [r["sampling_months"] for r in single if r["season_sampled"] == "unknown"] == ["NA"])
chk("no multi-season study is miscounted as single",
    all(len(r["sampling_months"].split(",")) > 1 or r["sampling_months"] == "NA"
        for r in S3 if r["multi_season"] == "yes"))
chk("every month field is a month list or NA",
    all(r["sampling_months"] == "NA" or
        all(m.isdigit() and 1 <= int(m) <= 12 for m in r["sampling_months"].split(","))
        for r in S3))

if not FAIL and not already:
    with open(P("TableS3_crossstudy.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, COLS, delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(S3)
    print("  wrote TableS3_crossstudy.tsv with three added columns")

print("\n  single-timepoint studies and the season each sampled:")
for r in sorted(single, key=lambda x: x["study"]):
    print("    %-14s %-22s months %-6s %s" % (r["study"], r["system"], r["sampling_months"], r["season_sampled"]))
print("\n%d checks, %d mismatches" % (NCHK[0], len(FAIL)))
if FAIL:
    for f in FAIL: print("   FAILED:", f)
    sys.exit(1)
