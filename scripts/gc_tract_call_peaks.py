"""
How many times each position of the V gene is covered by a called tract.

Two counts, deliberately kept apart
-----------------------------------
CALLS counts every transcript-level tract call, and DISTINCT counts each
(parent, donor, start, end) once.  They differ by a lot -- the busiest IGL
position is covered by 91 calls but only 13 distinct tracts -- because many
transcripts are clonally related and carry the same tract.  Plotting calls alone
would let one expanded clone look like a hotspot; plotting distinct tracts alone
would throw away the fact that a tract really was seen 91 times.  Both are drawn.

One panel per locus
-------------------
IGH pools seven parents, and they cannot be pooled in absolute coordinates: the
genes fall into two length classes whose CDR2 boundaries sit at 106 and ~144, so
an absolute trace is bimodal purely because of that.  The IGH panel is therefore
drawn relative to each gene's own CDR2 start.  That is exact rather than
approximate here -- FR2 is 39 bp in all seven parents, so CDR1 and CDR2 occupy
the SAME relative interval in every gene (CDR2 at 0..+30, CDR1 at -63..-39) and
the shaded blocks are real, not an average.

IGL has a single parent, so it is drawn in its own real coordinates; there is
nothing to align.
"""
import argparse
import csv
import sys
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from gc_conversion_peaks import landmarks, merge
from gc_tract_sequence_map import regions_of
from gc_lib import read_fasta
from gc_palette import save_figure, LOCUS, SEGMENT, GREY, GREY_DARK, INK, NO, YES


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def short(g):
    p = g.split(".")
    return p[1] if len(p) > 1 else g


def profiles(seqs, tracts):
    """per parent: (calls per position, distinct tracts per position, CDR2 start).

    calls    = one per transcript-level event (primary donor only)
    distinct = one per distinct tract POSITION, so clonally related transcripts
               carrying the same tract collapse to one
    """
    out = {}
    byp = defaultdict(list)
    for t in tracts:
        byp[t["parent"]].append(t)
    for p, ts in byp.items():
        n = len(seqs[p])
        calls = np.zeros(n)
        distinct = np.zeros(n)
        seen = set()
        for t in ts:
            s, e = int(t["start"]), min(int(t["end"]), n - 1)
            calls[s:e + 1] += 1
            # keyed on position, not donor: two candidate donors for one tract
            # are one event, and the point of this trace is where events land
            key = (t["start"], t["end"])
            if key not in seen:
                seen.add(key)
                distinct[s:e + 1] += 1
        lm = dict((lab, x) for x, lab in landmarks(seqs[p])[0])
        out[p] = (calls, distinct, lm.get("FR2|CDR2"))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--locus", action="append", required=True,
                    metavar="NAME=fasta,tracts")
    ap.add_argument("--out-figure", required=True)
    ap.add_argument("--out-table", required=True)
    args = ap.parse_args()

    data = []
    for spec in args.locus:
        name, rest = spec.split("=", 1)
        fa, tr = rest.split(",")
        seqs = read_fasta(fa)
        # primary_donor filter matters: a tract whose donor is ambiguous is
        # listed once per candidate donor, so 42 significant IGH rows are only
        # 28 events. Counting the rows would inflate IGH by ~50% and would also
        # disagree with panel C of the main figure, which filters here.
        tracts = [t for t in read_tsv(tr) if t.get("significant") == "True"
                  and t.get("primary_donor", "True") == "True"]
        data.append((name, seqs, profiles(seqs, tracts)))

    with open(args.out_table, "w") as fh:
        fh.write("locus\tparent\tposition\tposition_rel_cdr2\tn_calls\tn_distinct_tracts\n")
        for name, seqs, prof in data:
            for p, (c, d, cdr2) in prof.items():
                for i in range(len(c)):
                    if c[i]:
                        rel = i - cdr2 if cdr2 is not None else ""
                        fh.write(f"{name}\t{short(p)}\t{i}\t{rel}\t"
                                 f"{int(c[i])}\t{int(d[i])}\n")

    fig, axes = plt.subplots(len(data), 1, figsize=(11.0, 3.1 * len(data) + 0.9),
                             gridspec_kw={"hspace": 0.48}, squeeze=False)
    axes = axes[:, 0]
    FR2_LEN, CDR1_LEN, CDR2_LEN = 39, 24, 30

    for ax, (name, seqs, prof) in zip(axes, data):
        col = LOCUS.get(name, INK)
        multi = len(prof) > 1
        if multi:
            # anchored on CDR2, because absolute coordinates cannot pool genes
            span = 150
            calls = np.zeros(2 * span)
            distinct = np.zeros(2 * span)
            for p, (c, d, cdr2) in prof.items():
                if cdr2 is None:
                    continue
                for i in range(len(c)):
                    j = i - cdr2 + span
                    if 0 <= j < 2 * span:
                        calls[j] += c[i]
                        distinct[j] += d[i]
            x = np.arange(-span, span)
            cdrs = [(0, CDR2_LEN, "CDR2"), (-FR2_LEN - CDR1_LEN, -FR2_LEN, "CDR1")]
            xlab = "position relative to the start of CDR2 (bp)"
            sub = f"{len(prof)} parent genes pooled"
        else:
            p = next(iter(prof))
            calls, distinct, cdr2 = prof[p]
            x = np.arange(len(calls))
            cdrs = ([(cdr2, cdr2 + CDR2_LEN, "CDR2"),
                     (cdr2 - FR2_LEN - CDR1_LEN, cdr2 - FR2_LEN, "CDR1")]
                    if cdr2 is not None else [])
            xlab = f"position in the V gene (bp)"
            sub = f"parent gene {short(p)}"

        top = max(calls.max(), 1)
        for s0, e0, nm in cdrs:
            ax.axvspan(s0, e0, color=SEGMENT["J"], alpha=0.22, lw=0, zorder=0)
            ax.annotate(nm, ((s0 + e0) / 2, top * 1.10), ha="center", va="bottom",
                        fontsize=8, color=SEGMENT["J"], fontweight="bold")
        ax.fill_between(x, calls, color=col, alpha=0.50, lw=0, zorder=2,
                        label="tract calls")
        ax.plot(x, distinct, color=INK, lw=1.6, zorder=3, label="distinct tracts")
        ax.set_xlim(x[0], x[-1])
        ax.set_ylim(0, top * 1.26)
        ax.set_xlabel(xlab, fontsize=9)
        ax.set_ylabel("tracts covering", fontsize=9)
        ax.tick_params(labelsize=8)
        ax.legend(fontsize=8.5, frameon=False, loc="upper right")
        ax.set_title(f"$\\bf{{{name}}}$   {sub}", fontsize=10.5, loc="left",
                     x=0.0, color=col)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)

    fig.tight_layout()
    save_figure(fig, args.out_figure)
    print(f"{len(data)} panels", file=sys.stderr)


if __name__ == "__main__":
    main()
