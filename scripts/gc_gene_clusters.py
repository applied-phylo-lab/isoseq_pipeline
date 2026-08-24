"""
Cluster the V genes of a locus by sequence similarity and show the families.

Difference from the identity heatmap
------------------------------------
`gc_gene_similarity.py` orders genes by contig position, which shows WHERE the
similar genes are. This orders them by similarity instead, which shows WHICH
genes form families: hierarchical clustering on the same pairwise-identity
matrix, with the dendrogram drawn alongside the reordered heatmap and the
clusters boxed. A gene family is a block on the diagonal; a lineage-specific
expansion is a large one.

Handling pairs that do not align
--------------------------------
Some pairs share less than MIN_OVERLAP_BP and have no defined identity (white in
the heatmap). Clustering needs a complete distance matrix, so those are given the
maximum distance: two V genes that share almost no sequence ARE maximally
dissimilar, and a gene that aligns to nothing (IGL 6329906) correctly falls out
as its own singleton rather than being forced into a family.

Cutting the tree
----------------
Clusters are defined by cutting the dendrogram at a fixed identity
(`--cut-identity`, default 0.90): genes joined below 1 - 0.90 = 0.10 distance are
one cluster. The cut is drawn on the dendrogram so the grouping is not a black
box, and the cluster of each gene is written to the table with its RSS state and
expression, since a family that contains the one functional gene is the donor
array for that gene.
"""
import argparse
import csv
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle
from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform

from gc_gene_similarity import pairwise_identity, short, read_tsv
from gc_lib import read_fasta
from gc_palette import save_figure, SEGMENT, GREY, GREY_DARK, INK, LOCUS, cycle


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vgene-fasta", required=True)
    ap.add_argument("--locus", required=True)
    ap.add_argument("--rss-annotation", required=True)
    ap.add_argument("--usage-assignments")
    ap.add_argument("--cut-identity", type=float, default=0.90,
                    help="genes joined above this identity are one cluster")
    ap.add_argument("--out-figure", required=True)
    ap.add_argument("--out-table", required=True)
    args = ap.parse_args()

    seqs = read_fasta(args.vgene_fasta)
    rss = {r["gene"]: r for r in read_tsv(args.rss_annotation)
           if r["locus"] == args.locus}
    names = sorted([g for g in seqs if g in rss],
                   key=lambda g: int(rss[g]["pos"])) or sorted(seqs)
    mat, _cov = pairwise_identity(seqs, names)

    usage = {}
    if args.usage_assignments:
        for r in read_tsv(args.usage_assignments):
            if r["locus"] == args.locus:
                usage[r["best_gene"]] = usage.get(r["best_gene"], 0) + 1

    n = len(names)
    # distance = 1 - identity; unalignable pairs get the maximum distance so a
    # gene that shares nothing with the rest clusters alone rather than being
    # forced into whichever family it happens to touch.
    dist = 1.0 - mat
    np.fill_diagonal(dist, 0.0)
    dist[np.isnan(dist)] = 1.0
    dist = (dist + dist.T) / 2                 # enforce exact symmetry
    Z = hierarchy.linkage(squareform(dist, checks=False), method="average")

    cut = 1.0 - args.cut_identity
    labels = hierarchy.fcluster(Z, t=cut, criterion="distance")
    # order clusters by size, largest first, for stable colouring
    order_by_size = {c: i for i, (c, _) in enumerate(
        sorted(((c, np.sum(labels == c)) for c in set(labels)),
               key=lambda kv: -kv[1]))}
    n_clusters = len(set(labels))
    n_singletons = sum(1 for c in set(labels) if np.sum(labels == c) == 1)

    big = n > 40
    fig = plt.figure(figsize=(12.5, 10.5) if big else (10.5, 8.8))
    gs = fig.add_gridspec(1, 3, width_ratios=[0.20, 1.0, 0.045],
                          wspace=0.02, left=0.04, right=0.93,
                          top=0.90, bottom=0.14)

    # ── dendrogram ────────────────────────────────────────────────────────────
    axd = fig.add_subplot(gs[0, 0])
    dn = hierarchy.dendrogram(Z, orientation="left", no_labels=True,
                              color_threshold=cut, above_threshold_color=GREY_DARK,
                              ax=axd)
    axd.axvline(cut, color=INK, ls="--", lw=1.0)
    axd.set_xlim(dist.max() * 1.02, 0)
    axd.set_xticks([0, cut, round(dist.max(), 1)])
    axd.set_xticklabels([f"{1 - x:.2f}" for x in [0, cut, round(dist.max(), 1)]],
                        fontsize=7)
    axd.set_xlabel("identity", fontsize=8)
    axd.tick_params(labelsize=7)
    for sp in ("top", "right", "left"):
        axd.spines[sp].set_visible(False)
    axd.set_yticks([])

    idx = dn["leaves"]                          # dendrogram leaf order (top→bottom)
    ordered = mat[np.ix_(idx, idx)]

    # ── reordered heatmap ─────────────────────────────────────────────────────
    axm = fig.add_subplot(gs[0, 1])
    im = axm.imshow(ordered[::-1], cmap="magma_r", vmin=0.70, vmax=1.0,
                    interpolation="nearest", aspect="auto")
    # cluster boxes: contiguous runs of the same cluster along the leaf order
    lead = [labels[i] for i in idx]
    start = 0
    for k in range(1, n + 1):
        if k == n or lead[k] != lead[start]:
            if k - start > 1:
                y0 = n - k
                axm.add_patch(Rectangle((start - 0.5, y0 - 0.5), k - start, k - start,
                                        fill=False, edgecolor=SEGMENT["J"], lw=1.8))
            start = k
    axm.set_xticks([])
    axm.set_yticks([])
    axm.set_title(
        f"{args.locus} — V genes clustered by identity "
        f"({n} genes, {n_clusters} clusters at ≥ {args.cut_identity:.0%}: "
        f"{n_clusters - n_singletons} families + {n_singletons} singletons)",
        fontsize=12, fontweight="bold", loc="left", x=0.0, pad=10)

    cb = fig.colorbar(im, cax=fig.add_subplot(gs[0, 2]), extend="min")
    cb.set_label("pairwise identity", fontsize=9)
    cb.ax.tick_params(labelsize=8)

    # ── per-gene annotation strip along the bottom ────────────────────────────
    # cluster colour, RSS, and expression under each column, so a family that
    # holds the functional gene is identifiable at a glance
    palette = cycle()
    ccol = {c: palette[order_by_size[c] % len(palette)] for c in set(labels)}
    for col, leaf in enumerate(idx):
        g = names[leaf]
        # cluster colour bar
        if np.sum(labels == labels[leaf]) > 1:
            axm.add_patch(Rectangle((col - 0.5, n - 0.5), 1, n * 0.018,
                                    color=ccol[labels[leaf]], clip_on=False))
        if rss.get(g, {}).get("rss_state") == "rss_present":
            axm.plot([col], [n + n * 0.045], marker="o", ms=5, color=SEGMENT["V"],
                     markeredgecolor="black", markeredgewidth=0.5, clip_on=False)
        if usage.get(g):
            axm.plot([col], [n + n * 0.075], marker="s", ms=4.2,
                     color=LOCUS.get(args.locus, INK), clip_on=False)

    from matplotlib.lines import Line2D
    handles = [
        Rectangle((0, 0), 1, 1, color=GREY_DARK, label="cluster (≥2 genes)"),
        Line2D([], [], marker="o", ls="none", ms=6, color=SEGMENT["V"],
               markeredgecolor="black", label="has an RSS"),
        Line2D([], [], marker="s", ls="none", ms=6, color=LOCUS.get(args.locus, INK),
               label="best match for ≥1 transcript"),
    ]
    axm.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.03),
               ncol=3, fontsize=8.5, frameon=True, framealpha=0.95)

    with open(args.out_table, "w") as fh:
        fh.write("gene\tpos\tcluster\tcluster_size\trss_state\tn_transcripts\n")
        for g in names:
            i = names.index(g)
            c = labels[i]
            fh.write(f"{g}\t{rss.get(g, {}).get('pos', 'NA')}\t"
                     f"{order_by_size[c] + 1}\t{int(np.sum(labels == c))}\t"
                     f"{rss.get(g, {}).get('rss_state', 'NA')}\t{usage.get(g, 0)}\n")

    save_figure(fig, args.out_figure)

    print(f"{args.locus}: {n} genes -> {n_clusters} clusters at >= "
          f"{args.cut_identity:.0%} ({n_singletons} singletons)", file=sys.stderr)
    for c in sorted(set(labels), key=lambda c: -np.sum(labels == c)):
        members = [names[i] for i in range(n) if labels[i] == c]
        if len(members) == 1:
            continue
        has_rss = sum(1 for g in members if rss.get(g, {}).get("rss_state") == "rss_present")
        expr = sum(1 for g in members if usage.get(g))
        print(f"  cluster {order_by_size[c] + 1}: {len(members)} genes, "
              f"{has_rss} with RSS, {expr} expressed  "
              f"[{', '.join(short(g) for g in members[:6])}"
              f"{'...' if len(members) > 6 else ''}]", file=sys.stderr)


if __name__ == "__main__":
    main()
