"""Reconstruct the arm assignment of PRJNA1159715 from the archive metadata rather
than inheriting it from derived files, and compare the three counts that can
disagree: how many runs the archive holds, how many were downloaded, and how many
carry an arm label. Also checks whether any sample maps to more than one run.
"""
import os
import re
import sys
import urllib.request
from collections import Counter, defaultdict

STUDY = "PRJNA1159715"
LABELS = os.path.expanduser("~/fpar/data/units_v3/UNITS_v3_labels.tsv")
CORPUS = os.path.expanduser("~/fpar/corpus")
OUTDIR = os.path.expanduser("~/fpar/results")

FIELDS = ",".join([
    "run_accession", "sample_accession", "experiment_accession",
    "sample_alias", "sample_title", "library_name",
    "library_strategy", "base_count",
])
URL = ("https://www.ebi.ac.uk/ena/portal/api/filereport"
       "?accession=%s&result=read_run&fields=%s&format=tsv&limit=0" % (STUDY, FIELDS))

print("=" * 72)
print("RECONCILIATION: %s rebuilt from ENA" % STUDY)
print("=" * 72)

try:
    with urllib.request.urlopen(URL, timeout=120) as r:
        raw = r.read().decode("utf-8", "replace")
except Exception as e:
    print("ENA FETCH FAILED: %s" % e)
    sys.exit(1)

lines = [l for l in raw.rstrip("\n").split("\n") if l.strip()]
hdr = lines[0].split("\t")
rows = [dict(zip(hdr, l.split("\t"))) for l in lines[1:]]
print("ENA returned %d run rows" % len(rows))
print()

print("--- first 3 raw rows (to confirm the alias format, not assume it) ---")
for r in rows[:3]:
    print("  run=%s  sample=%s" % (r.get("run_accession"), r.get("sample_accession")))
    print("     alias=%r" % r.get("sample_alias"))
    print("     title=%r" % r.get("sample_title"))
    print("     libname=%r" % r.get("library_name"))
print()

# ---- Q0: does any SAMPLE map to >1 RUN? -------------------------------
per_sample = defaultdict(list)
for r in rows:
    per_sample[r.get("sample_accession", "")].append(r.get("run_accession", ""))
multi = {s: v for s, v in per_sample.items() if len(v) > 1}
print("--- Q0: sample-to-run mapping ---")
print("unique samples: %d   unique runs: %d" % (len(per_sample), len(rows)))
if multi:
    print("SAMPLES WITH >1 RUN (%d) -- 'samples' and 'runs' are NOT interchangeable here:" % len(multi))
    for s, v in sorted(multi.items())[:10]:
        print("   %s -> %s" % (s, ",".join(v)))
else:
    print("every sample maps to exactly one run")
print()

# ---- decode substrate + site from alias/title/library_name -------------
SUB = [("MP", "plastic"), ("NP", "surface_control"), ("RW", "bulk_control")]


def decode(r):
    """Return (site, substrate, source_field, raw) or (None, None, None, raw)."""
    for f in ("sample_alias", "library_name", "sample_title"):
        v = (r.get(f) or "").strip()
        if not v:
            continue
        for code, arm in SUB:
            m = re.search(r"[-_.]%s(?:[-_.]|$)" % code, v, re.IGNORECASE)
            if m:
                site = v[:m.start()]
                site = re.sub(r"^HKR[.\-_]?W?[-_.]?", "", site, flags=re.IGNORECASE)
                return site.strip("-_. ") or "UNKNOWN", arm, f, v
    return None, None, None, (r.get("sample_alias") or r.get("sample_title") or "")


dec = []
undecoded = []
for r in rows:
    site, arm, fld, raw_v = decode(r)
    if arm is None:
        undecoded.append((r.get("run_accession"), raw_v))
    dec.append({"run": r.get("run_accession"), "site": site, "arm": arm,
                "field": fld, "raw": raw_v,
                "sample": r.get("sample_accession"),
                "strategy": r.get("library_strategy"),
                "bases": r.get("base_count")})

print("--- Q1: ARCHIVE structure, decoded from scratch ---")
c = Counter(d["arm"] for d in dec)
for code, arm in SUB:
    print("  %-16s (%s): %d" % (arm, code, c.get(arm, 0)))
print("  undecodable            : %d" % len(undecoded))
if undecoded:
    print("  undecodable aliases (first 10):")
    for run, v in undecoded[:10]:
        print("     %s  %r" % (run, v))
print("  TOTAL in archive       : %d" % len(dec))
print()

fields_used = Counter(d["field"] for d in dec if d["field"])
print("  decoded from field(s): %s" % dict(fields_used))
print()

print("--- site x substrate matrix (the diagnostic) ---")
sites = sorted({d["site"] for d in dec if d["site"]})
arms = [a for _, a in SUB]
print("  %-22s %9s %17s %15s" % ("site", "plastic", "surface_control", "bulk_control"))
incomplete = []
for s in sites:
    cells = []
    for a in arms:
        n = sum(1 for d in dec if d["site"] == s and d["arm"] == a)
        cells.append(n)
    flag = "  <== MISSING" if 0 in cells else ""
    if 0 in cells:
        incomplete.append((s, dict(zip(arms, cells))))
    print("  %-22s %9d %17d %15d%s" % (s, cells[0], cells[1], cells[2], flag))
print("  sites: %d   sites with a missing arm: %d" % (len(sites), len(incomplete)))
if incomplete:
    print("  INCOMPLETE SITES:")
    for s, d in incomplete:
        print("     %s -> %s" % (s, d))
print()

# ---- Q2: on disk ------------------------------------------------------
print("--- Q2: which archive runs are ON DISK ---")
on_disk, missing_disk = [], []
for d in dec:
    f1 = os.path.join(CORPUS, d["run"], "%s_1.fastq.gz" % d["run"])
    (on_disk if os.path.exists(f1) and os.path.getsize(f1) > 0 else missing_disk).append(d)
print("  on disk: %d / %d" % (len(on_disk), len(dec)))
if missing_disk:
    print("  NOT DOWNLOADED (%d):" % len(missing_disk))
    for d in missing_disk:
        print("     %s  arm=%s  site=%s  bases=%s" % (d["run"], d["arm"], d["site"], d["bases"]))
print()

# ---- Q3: in UNITS_v3 --------------------------------------------------
print("--- Q3: which archive runs are IN UNITS_v3_labels.tsv ---")
lab = {}
if os.path.exists(LABELS):
    with open(LABELS) as f:
        h = f.readline().rstrip("\n").split("\t")
        for line in f:
            p = line.rstrip("\n").split("\t")
            row = dict(zip(h, p))
            if row.get("study_accession") == STUDY:
                lab[row.get("run_accession")] = row.get("arm")
else:
    print("  LABELS FILE NOT FOUND: %s" % LABELS)

print("  runs labelled for %s: %d" % (STUDY, len(lab)))
arch_runs = {d["run"] for d in dec}
missing_lab = [d for d in dec if d["run"] not in lab]
extra_lab = [r for r in lab if r not in arch_runs]
print("  in archive but NOT labelled: %d" % len(missing_lab))
for d in missing_lab:
    print("     %s  archive_arm=%s  site=%s  on_disk=%s"
          % (d["run"], d["arm"], d["site"],
             os.path.exists(os.path.join(CORPUS, d["run"], "%s_1.fastq.gz" % d["run"]))))
print("  labelled but NOT in archive: %d" % len(extra_lab))
for r in extra_lab:
    print("     %s  label_arm=%s" % (r, lab[r]))

disagree = [(d["run"], d["arm"], lab[d["run"]])
            for d in dec if d["run"] in lab and lab[d["run"]] != d["arm"]]
print("  ARM DISAGREEMENTS archive vs labels: %d" % len(disagree))
for run, a, b in disagree:
    print("     %s  archive=%s  labels=%s" % (run, a, b))
print()

# ---- summary -----------------------------------------------------------
print("=" * 72)
n_arch, n_disk, n_lab = len(dec), len(on_disk), len(lab)
print("  archive=%d   downloaded=%d   labelled=%d" % (n_arch, n_disk, n_lab))
if n_disk < n_arch:
    print("  %d run(s) are present at ENA but were not downloaded." % (n_arch - n_disk))
if n_lab < n_disk:
    print("  %d run(s) are on disk but carry no arm label." % (n_disk - n_lab))
if n_arch == n_disk == n_lab:
    print("  complete and consistent at %d runs." % n_arch)
print("=" * 72)
