"""
Worked examples of overlapping V gene annotations.

What the figure shows
---------------------
Some V "genes" in the annotation are not two genes at all.  Two entries can sit
at OVERLAPPING genomic coordinates on OPPOSITE strands, and over the overlap the
forward-strand sequence of one is byte-identical to the forward-strand sequence
of the other -- because it IS the same DNA, written into the file twice.  A V
gene's reverse complement still scores against V profiles, so an annotator that
scans both strands emits both.

This is a different thing from the tandem duplication that fills a V array.
Tandem duplicates are two genes at DIFFERENT coordinates that happen to be 90-99%
similar.  These are one locus at ONE set of coordinates, 100% identical to
itself, counted twice.

Why it matters, and how much
----------------------------
It matters for DENOMINATORS: "25 of 162 IGH V genes carry an RSS" counts 59 loci
twice, and merging the intervals gives 25 of 103.  It does NOT matter for any
tract call -- no donor->parent pair in either locus is an overlap pair -- nor for
expression, because the overlap is partial and the antisense copy covers less of
the transcript and so loses on coverage.
"""
import argparse
import csv
import sys
from collections import Counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow, Patch, Rectangle

from gc_lib import read_fasta
from gc_palette import save_figure, SEGMENT, GREY, GREY_DARK, INK, NO, YES

RC = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def revcomp(s):
    return s[::-1].translate(RC)


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vgene-fasta", required=True)
    ap.add_argument("--rss-annotation", required=True)
    ap.add_argument("--usage-assignments", required=True)
    ap.add_argument("--locus", required=True)
    ap.add_argument("--n-examples", type=int, default=6)
    ap.add_argument("--out-figure", required=True)
    ap.add_argument("--out-table", required=True)
    args = ap.parse_args()

    rss = {r["gene"]: r for r in read_tsv(args.rss_annotation)
           if r["locus"] == args.locus}
    seqs = read_fasta(args.vgene_fasta)
    usage = Counter(r["best_gene"] for r in read_tsv(args.usage_assignments)
                    if r["locus"] == args.locus)

    def genomic(n):
        """Forward-strand sequence, undoing however the entry was stored.

        A minus-strand entry is written as the reverse complement of the genome,
        so its index 0 is the RIGHT-hand coordinate. Comparing two entries
        without undoing that compares a sequence to its own reverse and reports
        them as unrelated.
        """
        return seqs[n] if rss[n]["strand"] == "+" else revcomp(seqs[n])

    iv = sorted((int(rss[n]["pos"]), int(rss[n]["pos"]) + len(seqs[n]), n)
                for n in seqs if n in rss)
    pairs = []
    for i in range(len(iv) - 1):
        for j in range(i + 1, len(iv)):
            if iv[j][0] >= iv[i][1]:
                break
            a, b = iv[i], iv[j]
            o, c = max(a[0], b[0]), min(a[1], b[1])
            same = genomic(a[2])[o - a[0]:c - a[0]] == genomic(b[2])[o - b[0]:c - b[0]]
            pairs.append((a, b, o, c, same))

    # distinct loci = overlapping intervals merged; needs no choice of which
    # member is "the real one", which is the point -- in 41 of 59 IGH pairs
    # neither member has an RSS, so there is no principled winner.
    merged = []
    for s, e, n in iv:
        if merged and s < merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
            merged[-1][2].append(n)
        else:
            merged.append([s, e, [n]])

    with open(args.out_table, "w") as fh:
        fh.write("gene_a\tstrand_a\tstart_a\tend_a\trss_a\ttranscripts_a\t"
                 "gene_b\tstrand_b\tstart_b\tend_b\trss_b\ttranscripts_b\t"
                 "overlap_bp\tidentical_over_overlap\n")
        for a, b, o, c, same in pairs:
            fh.write("\t".join(str(x) for x in [
                a[2], rss[a[2]]["strand"], a[0], a[1], rss[a[2]]["rss_state"],
                usage.get(a[2], 0),
                b[2], rss[b[2]]["strand"], b[0], b[1], rss[b[2]]["rss_state"],
                usage.get(b[2], 0), c - o, same]) + "\n")

    def has(n):
        return rss[n]["rss_state"] == "rss_present"

    # Show the informative cases first: pairs where the two members disagree on
    # RSS or on expression are the ones where keeping the wrong member would
    # actually change something.
    def interest(p):
        a, b, o, c, _ = p
        return (has(a[2]) != has(b[2])) * 2 + \
               (bool(usage.get(a[2], 0)) != bool(usage.get(b[2], 0)))
    show = sorted(pairs, key=lambda p: (-interest(p), -(p[3] - p[2])))[:args.n_examples]

    n = len(show)
    fig, axes = plt.subplots(n, 1, figsize=(11.0, 1.30 * n + 1.9), squeeze=False)
    axes = axes[:, 0]

    for ax, (a, b, o, c, same) in zip(axes, show):
        lo = min(a[0], b[0])
        hi = max(a[1], b[1])
        pad = (hi - lo) * 0.06
        # the shared block, drawn first so both arrows sit on top of it
        ax.add_patch(Rectangle((o, -0.95), c - o, 1.9, facecolor=SEGMENT["J"],
                               alpha=0.22, edgecolor="none", zorder=1))
        for k, (s, e, name) in enumerate((a, b)):
            y = 0.42 if k == 0 else -0.42
            fwd = rss[name]["strand"] == "+"
            col = SEGMENT["V"] if has(name) else GREY_DARK
            head = (e - s) * 0.16
            ax.add_patch(FancyArrow(
                s if fwd else e, y, (e - s - head) * (1 if fwd else -1), 0,
                width=0.26, head_width=0.40, head_length=head,
                length_includes_head=True, facecolor=col,
                edgecolor="black", linewidth=0.6, zorder=3))
            tx = usage.get(name, 0)
            ax.text(hi + pad * 0.4, y,
                    f"  {rss[name]['strand']}  {e - s} bp   "
                    f"{'RSS' if has(name) else 'no RSS'}"
                    f"{'   ' + str(tx) + ' transcripts' if tx else ''}",
                    va="center", ha="left", fontsize=8,
                    color=INK if has(name) else GREY_DARK)
        ax.annotate(f"{c - o} bp shared\n{'identical' if same else 'NOT identical'}",
                    ((o + c) / 2, 0), ha="center", va="center", fontsize=7.6,
                    color=SEGMENT["J"], fontweight="bold", zorder=5,
                    bbox=dict(boxstyle="round,pad=0.18", fc="white",
                              ec="none", alpha=0.88))
        ax.set_xlim(lo - pad, hi + pad * 5.2)
        ax.set_ylim(-1.05, 1.05)
        ax.set_yticks([])
        ax.set_xticks([lo, o, c, hi])
        ax.set_xticklabels([f"{v:,}" for v in (lo, o, c, hi)], fontsize=7.5)
        ax.tick_params(axis="x", length=2, pad=1.5)
        for sp in ("left", "right", "top"):
            ax.spines[sp].set_visible(False)

    n_same = sum(1 for p in pairs if p[4])
    axes[0].set_title(
        f"$\\bf{{{args.locus}}}$ — overlapping V gene annotations: the same DNA "
        f"entered twice, once per strand\n"
        f"{len(pairs)} overlapping pairs among {len(iv)} annotated entries; "
        f"{n_same}/{len(pairs)} identical over the shared block   ·   "
        f"merging intervals gives {len(merged)} distinct loci",
        fontsize=11, loc="left", x=0.0, pad=12)
    axes[-1].set_xlabel("contig position (bp)", fontsize=9)

    handles = [
        Patch(facecolor=SEGMENT["V"], edgecolor="black", label="entry has an RSS"),
        Patch(facecolor=GREY_DARK, edgecolor="black", label="entry has no RSS"),
        Patch(facecolor=SEGMENT["J"], alpha=0.22, label="shared genomic block"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, fontsize=8.5,
               frameon=True, framealpha=0.95, bbox_to_anchor=(0.5, 0.005))
    fig.tight_layout(rect=(0, 0.045, 1, 1))
    save_figure(fig, args.out_figure)

    print(f"{args.locus}: {len(pairs)} overlapping pairs, {n_same} identical over "
          f"the overlap; {len(iv)} entries -> {len(merged)} distinct loci",
          file=sys.stderr)
    both_rss = sum(1 for a, b, *_ in pairs if has(a[2]) and has(b[2]))
    one_rss = sum(1 for a, b, *_ in pairs if has(a[2]) != has(b[2]))
    print(f"  pairs where both have an RSS {both_rss}, exactly one {one_rss}, "
          f"neither {len(pairs) - both_rss - one_rss}", file=sys.stderr)


if __name__ == "__main__":
    main()
