"""Figure 2, the in silico recovery curve, the recovery coefficient per donor, and
the minimum detectable effect per season against the effects observed.

    python3 11_figure2.py [table_dir] [out_dir]
"""
import csv, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import sys
HERE = os.path.dirname(os.path.abspath(__file__))
D1 = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "..", "tables")
D = sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, "..", "..", "figures")
os.makedirs(D, exist_ok=True)
REC = os.path.join(D1, "TableS8_spikein_recovery.tsv")
SRC = os.path.join(D1, "TableS1_blocks.tsv")
FLOOR = 0.16
SE = ["Winter", "Spring", "Summer", "Autumn"]

sec, cur = [], []
for ln in open(REC):
    if ln.strip() == "":
        if cur:
            sec.append(cur)
        cur = []
    else:
        cur.append(ln.rstrip("\n").split("\t"))
if cur:
    sec.append(cur)
print("sections in recovery table: %d" % len(sec))


def tab(block):
    h = block[0]
    return [dict(zip(h, r)) for r in block[1:]]


pts = tab(sec[0])
fit = tab(sec[1]) if len(sec) > 1 else []
mde = tab(sec[2]) if len(sec) > 2 else []
DON = ["Pichia", "Yarrowia", "Fusarium"]
COL = {"Pichia": "#9c27a8", "Yarrowia": "#1b7f5c", "Fusarium": "#d9a43b"}
INK, GREY, MP_C = "#22252a", "#b9bcc2", "#c1443a"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8,
                     "axes.linewidth": 0.7, "axes.spines.top": False,
                     "axes.spines.right": False, "legend.frameon": False})
fig, ax = plt.subplots(1, 3, figsize=(7.6, 3.1))
fig.subplots_adjust(wspace=0.42, left=0.085, right=0.985,
                    top=0.93, bottom=0.30)

A = ax[0]
lo, hi = 0.7, 6e4
A.plot([lo, hi], [lo, hi], ls=":", c=INK, lw=0.9, zorder=1)
for d in DON:
    x = [float(r["injected_pairs"]) for r in pts if r["donor"] == d
         and float(r["injected_pairs"]) > 0 and float(r["observed_reads"]) > 0]
    y = [float(r["observed_reads"]) for r in pts if r["donor"] == d
         and float(r["injected_pairs"]) > 0 and float(r["observed_reads"]) > 0]
    A.scatter(x, y, s=16, c=COL[d], lw=0, alpha=0.85, label="$%s$" % d, zorder=3)
k = {r["donor"]: float(r["recovery_k"]) for r in fit} if fit else {}
if k:
    xs = np.array([lo, hi])
    A.plot(xs, xs * np.mean(list(k.values())), "-", c=INK, lw=1.2, zorder=2)
A.set_xscale("log")
A.tick_params(which="minor", length=0)
A.set_yscale("log")
A.set_xlim(lo, hi)
A.set_ylim(lo, hi)
A.set_xlabel("Reads injected")
A.set_ylabel("Reads recovered")
A.set_title("a", loc="left",
            fontweight="bold", fontsize=11)
_h, _l = A.get_legend_handles_labels()
_h += [Line2D([], [], ls=":", c=INK, lw=0.9, label="Perfect recovery (1:1)"),
       Line2D([], [], ls="-", c=INK, lw=1.2, label="Mean fitted recovery")]
A.legend(handles=_h, fontsize=6.0, loc="upper left", handletextpad=0.5,
         labelspacing=0.28, borderpad=0.2)

B = ax[1]
for i, d in enumerate(DON):
    f = [r for r in fit if r["donor"] == d]
    if not f:
        continue
    kk = float(f[0]["recovery_k"])
    se = float(f[0]["se"])
    B.errorbar([i], [kk], yerr=[1.96 * se], fmt="o", ms=7, c=COL[d],
               lw=1.4, capsize=5)
    B.text(i, kk + 0.045, "%.3f" % kk, ha="center", fontsize=7, c=INK)
B.axhline(1.0, c=GREY, lw=0.8, ls="--")
B.set_xticks(range(len(DON)))
B.set_xticklabels(["$%s$" % d for d in DON], fontsize=7.5)
B.set_ylim(0, 1.12)
B.set_ylabel("Fraction of reads recovered")
B.set_title("b", loc="left",
            fontweight="bold", fontsize=11)

C = ax[2]
src = list(csv.DictReader(open(SRC), delimiter="\t"))
obs = {}
for s in SE:
    v = [float(r["d"]) for r in src if r["season"] == s]
    if v:
        obs[s] = float(np.mean(v))
md = {r["season"]: float(r["mdes_log2"]) for r in mde} if mde else {}
yv = np.arange(len(SE))
for i, s in enumerate(SE):
    if s not in md:
        continue
    C.barh(i, 2 * md[s], 0.5, left=-md[s], color=GREY, alpha=0.75)
    C.plot([obs.get(s, 0)], [i], "o", ms=7,
           c=(MP_C if abs(obs.get(s, 0)) > md[s] else INK), zorder=3)
C.axvline(0, c=INK, lw=0.7)
C.set_yticks(yv)
C.set_yticklabels(SE, fontsize=7.5)
C.invert_yaxis()
C.legend(handles=[
    Line2D([], [], ls="-", lw=6, c=GREY, alpha=0.5,
           label="Detectable at 80% power ($\\pm$)"),
    Line2D([], [], marker="o", ls="", ms=5, c=INK,
           label="Observed, below threshold"),
    Line2D([], [], marker="o", ls="", ms=5, c=MP_C,
           label="Observed, above threshold")],
    fontsize=5.8, loc="lower center", bbox_to_anchor=(0.5, -0.46),
    frameon=False, ncol=1, handletextpad=0.5, borderpad=0.1,
    labelspacing=0.28)
_lim = max(max(md.values()), max(abs(v) for v in obs.values())) * 1.18
C.set_xlim(-_lim, _lim)
C.set_xlabel("Effect size (log$_2$)")
C.set_title("c", loc="left",
            fontweight="bold", fontsize=11)

fig.tight_layout(pad=1.0)
for ext, kw in (("png", {"dpi": 400}), ("pdf", {}),
                ("tiff", {"dpi": 400, "pil_kwargs": {"compression": "tiff_lzw"}})):
    fig.savefig(os.path.join(D, "Figure2.%s" % ext), bbox_inches="tight", **kw)
print("wrote Figure2.png/.pdf/.tiff in %s" % D)
print("recovery k: %s" % ", ".join("%s %.4f" % (d, k[d]) for d in DON if d in k))
print("")
print("=== VERIFICATION ===")
print("a  legend names both lines: perfect recovery (1:1) and mean fitted recovery")
print("c  detectable bound is mirrored about zero (two-sided):")
for _s in SE:
    if _s in md:
        _o = obs.get(_s, 0.0)
        print("      %-7s observed %+6.3f   bound +-%.3f   |obs|<bound: %s"
              % (_s, _o, md[_s], "yes" if abs(_o) < md[_s] else "NO"))
print("      expected: summer alone exceeds its bound (+3.050 vs 2.987);")
print("      winter, spring and autumn fall inside theirs, so those are")
print("      bounded nulls rather than demonstrated absences")
print("observed: %s" % ", ".join("%s %+.2f" % (s, obs[s]) for s in SE if s in obs))
print("seasons where observed exceeds MDES: %s"
      % ", ".join(s for s in SE if s in md and obs.get(s, 0) > md[s]))

print("Figure 2 written to %s" % D)
