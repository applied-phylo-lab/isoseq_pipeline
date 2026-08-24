"""
All-vs-all identity between the V genes of a locus.

Why this figure exists
----------------------
Every downstream statement of the form "this transcript came from gene X" rests
on gene X being distinguishable from gene Y in the first place.  That is an
assumption about the REFERENCE, not about the data, and it is testable directly:
align every V gene to every other one and look at how close the closest pair is.

The two loci behave completely differently, which is the point of the figure.
IGL's 23 genes top out around 93% identity, so a ~300 bp transcript has ~20
diagnostic positions to work with and the assignment is decisive.  IGH's 162
genes contain near-twins -- pairs separated by two or three bases over the whole
gene -- and for those the assignment is arbitrary no matter how good the data
are.  Ordering the matrix by contig position also shows WHERE the similar genes
are: near-identical pairs are adjacent, i.e. recent tandem duplicates.

Identity is computed over aligned, ungapped columns only, so it is "of the
positions these two genes share, how many agree" -- the same quantity the
assignment step maximises, and not diluted by indels the aligner had to open.
"""
import argparse
import itertools
import sys
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from Bio import Align
from matplotlib.lines import Line2D

from gc_lib import read_fasta
from gc_palette import save_figure, SEGMENT, GREY_DARK, INK, LOCUS


def read_tsv(path):
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(hdr, ln.rstrip("\n").split("\t")))
                for ln in fh if ln.strip()]


def short(name):
    """Gene label: the contig position, which is what orders the array."""
    parts = name.split(".")
    return parts[1] if len(parts) > 1 else name


RC = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")

# Shortest overlap that counts as real homology. Below this an alignment is a
# chance core rather than a shared gene, and its identity is not interpretable.
MIN_OVERLAP_BP = 100


def revcomp(s):
    return s[::-1].translate(RC)


def _local(al, a, b):
    """(identity over the aligned region, length of that region).

    LOCAL, not global, and normalised by the aligned columns rather than by gene
    length -- the standard BLAST-style percent identity. Both choices are forced
    by how this V gene FASTA is written out:

      * entries are not all in coding orientation. Some genes are stored as the
        reverse complement of others, so every pair has to be tried both ways.
      * entry boundaries are not consistent. Neighbouring genes can be offset by
        >100 bp, so two genuinely homologous V genes may share only two thirds of
        their length. A global aligner charges for the overhang and reports ~0.67
        for a pair that is ~0.90 identical where they actually overlap.

    Getting either wrong splits the IGH array into two apparent families sitting
    at chance similarity to each other, which is an artefact of the file rather
    than anything about the genes.
    """
    aln = al.align(a, b)[0]
    sa, sb = str(aln[0]), str(aln[1])
    same = tot = 0
    for x, y in zip(sa, sb):
        if x == "-" or y == "-":
            continue
        tot += 1
        same += x == y
    return (same / tot if tot else np.nan), tot


def pairwise_identity(seqs, names):
    """All-vs-all identity, orientation-agnostic and offset-tolerant.

    Returns (identity matrix, overlap matrix). Pairs overlapping by less than
    MIN_OVERLAP_BP are left NaN: no homology worth quoting an identity for.
    """
    al = Align.PairwiseAligner(mode="local", match_score=1, mismatch_score=-1,
                              open_gap_score=-4, extend_gap_score=-1)
    n = len(names)
    mat = np.full((n, n), np.nan)
    cov = np.zeros((n, n))
    np.fill_diagonal(mat, 1.0)
    rcs = {g: revcomp(seqs[g]) for g in names}
    for i, j in itertools.combinations(range(n), 2):
        a, b = names[i], names[j]
        best = max(_local(al, seqs[a], seqs[b]),
                   _local(al, seqs[a], rcs[b]),
                   key=lambda t: (t[1] >= MIN_OVERLAP_BP, t[0] * t[1]))
        v, ov = best
        if ov < MIN_OVERLAP_BP:
            v = np.nan
        mat[i, j] = mat[j, i] = v
        cov[i, j] = cov[j, i] = ov
    return mat, cov


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vgene-fasta", required=True)
    ap.add_argument("--locus", required=True)
    ap.add_argument("--rss-annotation", required=True)
    ap.add_argument("--usage-assignments",
                    help="unconstrained_assignments.tsv; marks which genes are "
                         "the best match for at least one transcript")
    ap.add_argument("--out-figure", required=True)
    ap.add_argument("--out-table", required=True)
    args = ap.parse_args()

    seqs = read_fasta(args.vgene_fasta)
    rss = {r["gene"]: r for r in read_tsv(args.rss_annotation)
           if r["locus"] == args.locus}
    # Order by contig position, so adjacency in the matrix is adjacency in the
    # genome and a block of high identity reads as a tandem duplication.
    names = sorted([g for g in seqs if g in rss],
                   key=lambda g: int(rss[g]["pos"]))
    if not names:
        names = sorted(seqs)
    mat, cov = pairwise_identity(seqs, names)

    usage = defaultdict(int)
    if args.usage_assignments:
        for r in read_tsv(args.usage_assignments):
            if r["locus"] == args.locus:
                usage[r["best_gene"]] += 1

    off = mat.copy()
    np.fill_diagonal(off, np.nan)

    # Some V gene annotations OVERLAP each other in genomic coordinates -- in IGH,
    # 59 pairs do, every one of them on opposite strands. Those two entries are
    # not two genes; they are the same DNA annotated twice, once per strand, so
    # an "identity" between them is 1.0 by construction and says nothing about
    # how distinguishable the array's genes are. They are kept in the matrix (the
    # duplication is worth seeing) but excluded from every summary statistic.
    spans = {}
    for g in names:
        r = rss.get(g)
        if r and r.get("pos", "NA") != "NA":
            spans[g] = (int(r["pos"]), int(r["pos"]) + len(seqs[g]))
    dup = np.zeros_like(off, dtype=bool)
    for i, a in enumerate(names):
        for j, b in enumerate(names):
            if i >= j or a not in spans or b not in spans:
                continue
            (s1, e1), (s2, e2) = spans[a], spans[b]
            if min(e1, e2) - max(s1, s2) > 0:
                dup[i, j] = dup[j, i] = True

    clean = np.where(dup, np.nan, off)

    def safe_max(row):
        return np.nan if np.all(np.isnan(row)) else np.nanmax(row)

    def safe_arg(row):
        return None if np.all(np.isnan(row)) else int(np.nanargmax(row))

    nn = np.array([safe_max(clean[i]) for i in range(len(names))])
    nn_all = np.array([safe_max(off[i]) for i in range(len(names))])

    with open(args.out_table, "w") as fh:
        fh.write("gene\tpos\trss_state\tn_transcripts\tnearest_neighbour\t"
                 "nearest_identity\tn_overlapping_annotations\n")
        for i, g in enumerate(names):
            j = safe_arg(clean[i])
            fh.write(f"{g}\t{rss.get(g, {}).get('pos', 'NA')}\t"
                     f"{rss.get(g, {}).get('rss_state', 'NA')}\t{usage.get(g, 0)}\t"
                     f"{names[j] if j is not None else 'NA'}\t"
                     f"{nn[i]:.4f}\t{int(dup[i].sum())}\n")

    n = len(names)
    big = n > 40
    fig, ax = plt.subplots(figsize=(11.5, 9.6) if big else (9.6, 8.2))

    # These are all real paralogues -- after deduplication there are no
    # chance-level pairs left, and every off-diagonal value sits between ~0.63
    # and ~0.99. Anchoring the ramp at chance (0.52) would therefore spend half
    # its range on values that never occur and compress every real difference
    # into the top of the scale. Instead the ramp is stretched over the range
    # that is actually present, so a 0.86 pair and a 0.95 pair are visibly
    # different colours. vmin is the ~1st percentile floored to a round value, so
    # a single outlier pair does not pull the whole scale down; a shared vmin
    # across both loci keeps the IGH and IGL panels comparable.
    # 0.70 sits just below the 1st percentile of both loci (IGH 0.75, IGL 0.73),
    # so the ramp covers 0.70-1.00 and the handful of more distant pairs land on
    # the extend triangle. Fixed rather than per-figure so the IGH and IGL panels
    # use the same scale and can be read against each other.
    im = ax.imshow(mat, cmap="magma_r", vmin=0.70, vmax=1.0,
                   interpolation="nearest")

    cb = fig.colorbar(im, ax=ax, fraction=0.040, pad=0.02, extend="min")
    cb.set_label("identity over the aligned region (local, either orientation)",
                 fontsize=9)
    cb.ax.tick_params(labelsize=8)

    step = 1 if n <= 30 else (5 if n <= 80 else 10)
    ticks = list(range(0, n, step))
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    lab = [short(names[i]) for i in ticks]
    ax.set_xticklabels(lab, rotation=90, fontsize=7 if not big else 6)
    ax.set_yticklabels(lab, fontsize=7 if not big else 6)

    # RSS state and expression as marker strips outside the matrix. Both are
    # properties of a gene, not of a pair, so they do not belong in the cells --
    # but they are what turns "these two are similar" into "and it matters",
    # because a near-twin only causes trouble if transcripts are assigned to it.
    for i, g in enumerate(names):
        has = rss.get(g, {}).get("rss_state") == "rss_present"
        if has:
            ax.plot([-1.6], [i], marker="o", ms=3.6, color=SEGMENT["V"],
                    markeredgecolor="none", clip_on=False, zorder=5)
        if usage.get(g):
            ax.plot([-3.0], [i], marker="s", ms=3.2,
                    color=LOCUS.get(args.locus, INK),
                    markeredgecolor="none", clip_on=False, zorder=5)

    worst_i = int(np.nanargmax(nn))
    worst_j = safe_arg(clean[worst_i])
    ax.set_title(
        f"{args.locus} — all-vs-all V gene identity ({n} genes, "
        f"{n * (n - 1) // 2:,} pairs)\n"
        f"closest distinct pair {nn[worst_i]:.3f} "
        f"({short(names[worst_i])} / {short(names[worst_j])})   ·   "
        f"{int(np.nansum(nn >= 0.95))} genes have a neighbour ≥ 0.95"
        + (f"   ·   {int(dup.sum() // 2)} pairs are the SAME DNA annotated on "
           f"both strands" if dup.any() else ""),
        fontsize=11, fontweight="bold", loc="left", x=0.0, pad=10)
    ax.set_xlabel("V gene (contig position, ordered)", fontsize=9)

    handles = [
        Line2D([], [], marker="o", ls="none", ms=6, color=SEGMENT["V"],
               label="has an RSS"),
        Line2D([], [], marker="s", ls="none", ms=6,
               color=LOCUS.get(args.locus, INK),
               label="best match for ≥1 transcript"),
    ]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.10),
              ncol=2, fontsize=8.5, frameon=True, framealpha=0.95)

    fig.tight_layout()
    save_figure(fig, args.out_figure)

    print(f"{args.locus}: {n} annotated entries, "
          f"{int(dup.sum() // 2)} overlapping (same-DNA) pairs involving "
          f"{int((dup.any(axis=1)).sum())} entries", file=sys.stderr)
    print(f"  closest DISTINCT pair {nn[worst_i]:.4f}", file=sys.stderr)
    for thr in (0.90, 0.95, 0.97, 0.99):
        print(f"  neighbour ≥ {thr}: excluding same-DNA pairs "
              f"{int(np.nansum(nn >= thr))}/{n}   "
              f"(including them {int(np.nansum(nn_all >= thr))}/{n})",
              file=sys.stderr)


if __name__ == "__main__":
    main()
