"""
The sequence of each recurrent conversion window, and where it sits in the gene.

Two things at once
------------------
The top panel is the parent V gene drawn to scale and partitioned into its
framework and CDR regions, with the conversion windows marked on it -- that is
the "where".  Every panel below is one window at single-base resolution: the
parent sequence, its translation, and every donor that converts into it, with
the bases that differ from the parent picked out.

Why the donors are worth seeing base by base
--------------------------------------------
The recurrent IGL window takes sequence from six different donor genes, and the
question that follows is whether they deliver the same change or different ones.
Rows of aligned letters answer that immediately, and no summary statistic does:
the donors fall into two families that differ at the 3' end of the window, which
is visible at a glance and invisible in any composition measure.

Region boundaries come from protein-level framework motifs (see
gc_conversion_peaks.py) rather than nucleotide regexes, because synonymous
variation between paralogues breaks the nucleotide patterns.
"""
import argparse
import csv
import sys
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio import Align
from Bio.Seq import Seq
from matplotlib.patches import Patch, Rectangle

from gc_conversion_peaks import landmarks, merge
from gc_lib import read_fasta
from gc_palette import save_figure, LOCUS, SEGMENT, GREY, GREY_DARK, INK, NO, YES

REGION_C = {"FR1": GREY, "CDR1": SEGMENT["J"], "FR2": GREY,
            "CDR2": SEGMENT["J"], "FR3": GREY, "CDR3": SEGMENT["J"]}


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def short(g):
    p = g.split(".")
    return p[1] if len(p) > 1 else g


def regions_of(seq):
    """[(start, end, name)] partitioning the gene into FR/CDR blocks."""
    lm = dict((lab, x) for x, lab in landmarks(seq)[0])
    a = lm.get("CDR1|FR2")
    b = lm.get("FR2|CDR2")
    c = lm.get("FR3|CDR3")
    out = []
    if a is not None:
        # CDR1 is the ~24 bp immediately before FR2; the FR1/CDR1 boundary is the
        # conserved cysteine, which is not reliably found in every pseudogene, so
        # it is taken as a fixed offset rather than guessed per gene.
        out += [(0, max(0, a - 24), "FR1"), (max(0, a - 24), a, "CDR1")]
    if a is not None and b is not None:
        out.append((a, b, "FR2"))
    if b is not None and c is not None:
        # CDR2 runs from the FR2 end to roughly a third of the way to the FR3 end
        out += [(b, min(c, b + 30), "CDR2"), (min(c, b + 30), c, "FR3")]
    if c is not None:
        out.append((c, len(seq), "CDR3"))
    return [(s, e, n) for s, e, n in out if e > s]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vgene-fasta", required=True)
    ap.add_argument("--tracts", required=True)
    ap.add_argument("--locus", required=True)
    ap.add_argument("--parent", help="restrict to one parent (default: the busiest)")
    ap.add_argument("--min-calls", type=int, default=1)
    ap.add_argument("--out-figure", required=True)
    ap.add_argument("--out-table", required=True)
    args = ap.parse_args()

    v = read_fasta(args.vgene_fasta)
    tracts = [t for t in read_tsv(args.tracts) if t.get("significant") == "True"]
    byp = defaultdict(list)
    for t in tracts:
        byp[t["parent"]].append((int(t["start"]), int(t["end"])))
    if args.parent:
        P = next(g for g in byp if short(g) == args.parent)
    else:
        P = max(byp, key=lambda g: len(byp[g]))
    seq = v[P]
    wins = merge(byp[P])

    calls = defaultdict(int)
    donors = defaultdict(set)
    for t in tracts:
        if t["parent"] != P:
            continue
        for s, e in wins:
            if s <= int(t["start"]) <= e:
                calls[(s, e)] += 1
                donors[(s, e)].add(t["donor"])
    wins = [(s, e) for s, e in wins if calls[(s, e)] >= args.min_calls]

    # donors projected onto parent coordinates
    al = Align.PairwiseAligner(mode="global", open_gap_score=-10,
                              extend_gap_score=-0.5, match_score=2, mismatch_score=-1)
    proj = {}
    for g in {d for w in wins for d in donors[w]}:
        aln = al.align(seq, v[g])[0]
        dp = [None] * len(seq)
        for (ps, pe), (ds, de) in zip(aln.aligned[0], aln.aligned[1]):
            for o in range(pe - ps):
                dp[ps + o] = v[g][ds + o]
        proj[g] = dp

    regs = regions_of(seq)

    def region_at(i):
        for s, e, n in regs:
            if s <= i < e:
                return n
        return "?"

    with open(args.out_table, "w") as fh:
        fh.write("locus\tparent\tstart\tend\tregion\tn_calls\tn_donors\trow\tsequence\n")
        for s, e in wins:
            fh.write(f"{args.locus}\t{short(P)}\t{s}\t{e}\t{region_at(s)}\t"
                     f"{calls[(s,e)]}\t{len(donors[(s,e)])}\tparent\t{seq[s:e+1]}\n")
            for g in sorted(donors[(s, e)], key=short):
                fh.write(f"{args.locus}\t{short(P)}\t{s}\t{e}\t{region_at(s)}\t"
                         f"{calls[(s,e)]}\t{len(donors[(s,e)])}\t{short(g)}\t"
                         + "".join(proj[g][i] or "-" for i in range(s, e + 1)) + "\n")

    rows_needed = sum(2 + len(donors[w]) for w in wins)
    fig = plt.figure(figsize=(13.0, 1.5 + 0.30 * rows_needed + 0.42 * len(wins)))
    gs = fig.add_gridspec(len(wins) + 1, 1,
                          height_ratios=[1.15] + [(2 + len(donors[w])) * 0.42
                                                  for w in wins],
                          hspace=0.62, left=0.155, right=0.985,
                          top=0.965, bottom=0.075)

    # ── where the windows land ────────────────────────────────────────────────
    ax = fig.add_subplot(gs[0])
    for s, e, n in regs:
        ax.add_patch(Rectangle((s, 0.18), e - s, 0.30,
                               facecolor=REGION_C.get(n, GREY),
                               alpha=0.85 if n.startswith("CDR") else 0.45,
                               edgecolor="white", lw=0.8))
        if e - s > 14:
            ax.text((s + e) / 2, 0.33, n, ha="center", va="center", fontsize=7,
                    color="white" if n.startswith("CDR") else INK,
                    fontweight="bold")
    for s, e in wins:
        ax.add_patch(Rectangle((s, 0.56), e - s + 1, 0.20, facecolor=YES,
                               edgecolor="none"))
        ax.annotate(f"{calls[(s,e)]} calls · {len(donors[(s,e)])} donors",
                    ((s + e) / 2, 0.80), ha="center", va="bottom", fontsize=7,
                    color=YES, fontweight="bold")
    ax.set_xlim(-4, len(seq) + 4)
    ax.set_ylim(0.10, 1.05)
    ax.set_yticks([])
    ax.set_xlabel(f"position in {args.locus} parent {short(P)} (bp)", fontsize=8.5)
    ax.tick_params(labelsize=7.5)
    for sp in ("left", "right", "top"):
        ax.spines[sp].set_visible(False)

    # ── the sequence of each window ───────────────────────────────────────────
    # the gene's real reading frame, not an assumed frame 0: two of the parents
    # in this set are in frame 1 and a fixed frame silently mistranslates them
    frame, aa_full = __import__("gc_conversion_peaks").frame_and_protein(seq)
    for k, (s, e) in enumerate(wins):
        ax = fig.add_subplot(gs[k + 1])
        w = e - s + 1
        names = ["parent " + short(P)] + [short(g) for g in
                                          sorted(donors[(s, e)], key=short)]
        seqs = [seq[s:e + 1]] + ["".join(proj[g][i] or "-" for i in range(s, e + 1))
                                 for g in sorted(donors[(s, e)], key=short)]
        for r, (nm, sq) in enumerate(zip(names, seqs)):
            y = len(names) - r
            ax.text(-0.012, y, nm, ha="right", va="center", fontsize=7.6,
                    fontweight="bold" if r == 0 else "normal",
                    color=INK if r == 0 else GREY_DARK, family="monospace",
                    transform=ax.get_yaxis_transform())
            for i, b in enumerate(sq):
                diff = r > 0 and b != seqs[0][i]
                if diff:
                    ax.add_patch(Rectangle((i - 0.47, y - 0.42), 0.94, 0.84,
                                           facecolor=YES, alpha=0.32,
                                           edgecolor="none"))
                ax.text(i, y, b, ha="center", va="center", fontsize=9.0,
                        family="monospace",
                        color=INK if (r == 0 or diff) else GREY,
                        fontweight="bold" if diff else "normal")
        # translation of the parent row, under its own codons
        for i in range(s, e + 1):
            if (i - frame) % 3 == 0 and 0 <= (i - frame) // 3 < len(aa_full):
                ax.text(i - s + 1, 0.30, aa_full[(i - frame) // 3], ha="center",
                        va="center", fontsize=8.0, color=INK, family="monospace")
        ax.text(-0.012, 0.30, "aa", ha="right", va="center", fontsize=7.4,
                color=INK, family="monospace", transform=ax.get_yaxis_transform())
        ax.set_xlim(-0.6, w - 0.4)
        ax.set_ylim(-0.1, len(names) + 0.75)
        ax.set_xticks([])
        ax.set_yticks([])
        for sp in ("left", "right", "top", "bottom"):
            ax.spines[sp].set_visible(False)
        ax.set_title(f"{s}–{e}  ({w} bp, {region_at(s)})   "
                     f"{calls[(s,e)]} calls from {len(donors[(s,e)])} donor"
                     f"{'s' if len(donors[(s,e)]) != 1 else ''}",
                     fontsize=8.6, loc="left", x=0.0, color=LOCUS.get(args.locus, INK))

    handles = [Patch(facecolor=SEGMENT["J"], alpha=0.85, label="CDR"),
               Patch(facecolor=GREY, alpha=0.45, label="framework"),
               Patch(facecolor=YES, label="conversion window / base differing from parent")]
    fig.legend(handles=handles, loc="lower center", ncol=3, fontsize=8,
               frameon=True, framealpha=0.95, bbox_to_anchor=(0.5, 0.004))
    save_figure(fig, args.out_figure)
    print(f"{args.locus} {short(P)}: {len(wins)} windows drawn", file=sys.stderr)


if __name__ == "__main__":
    main()
