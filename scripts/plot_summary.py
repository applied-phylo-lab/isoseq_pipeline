"""
Generate summary plots:
  1. Per-gene exact match counts (bar chart, faceted by locus)
  2. Identity distribution for rank-1 vs rank-2 hits (gene conversion signal)
  3. Delta-identity heatmap: rank1 - rank2 per transcript (lower delta = more
     gene-conversion-like, i.e. multiple similar donors)
"""
import argparse
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pandas as pd
import numpy as np
from gc_palette import save_figure



def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exact-matches", required=True)
    parser.add_argument("--top-alignments", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    exact = pd.read_csv(args.exact_matches, sep="\t")
    top = pd.read_csv(args.top_alignments, sep="\t")

    fig = plt.figure(figsize=(16, 14))
    gs = gridspec.GridSpec(3, 2, figure=fig, hspace=0.45, wspace=0.35)

    # ── Plot 1: Exact match counts per gene, coloured by locus ────────────────
    ax1 = fig.add_subplot(gs[0, :])
    exact_hits = exact[exact["transcripts_with_exact_match"] > 0].copy()
    if not exact_hits.empty:
        loci = exact_hits["locus"].unique()
        palette = {l: c for l, c in zip(loci, plt.cm.Set2.colors)}
        colors = exact_hits["locus"].map(palette)
        ax1.bar(range(len(exact_hits)), exact_hits["transcripts_with_exact_match"],
                color=colors, edgecolor="white", linewidth=0.4)
        ax1.set_xticks(range(len(exact_hits)))
        ax1.set_xticklabels(exact_hits["gene"], rotation=90, fontsize=6)
        ax1.set_ylabel("# IG transcripts with exact match")
        ax1.set_title("Germline V genes found verbatim in transcripts (functional candidates)")
        handles = [plt.Rectangle((0, 0), 1, 1, color=palette[l]) for l in loci]
        ax1.legend(handles, loci, title="Locus", loc="upper right")
    else:
        ax1.text(0.5, 0.5, "No exact matches found", ha="center", va="center",
                 transform=ax1.transAxes, fontsize=12)
        ax1.set_title("Germline V genes found verbatim in transcripts")

    # ── Plot 2: Identity distribution rank1 vs rank2 ──────────────────────────
    ax2 = fig.add_subplot(gs[1, 0])
    rank1 = top[top["rank"] == 1]["identity"]
    rank2 = top[top["rank"] == 2]["identity"]
    bins = np.linspace(0, 1, 51)
    ax2.hist(rank1, bins=bins, alpha=0.7, label="Rank 1 (best)", color="#2196F3")
    ax2.hist(rank2, bins=bins, alpha=0.7, label="Rank 2", color="#FF9800")
    ax2.set_xlabel("Alignment identity")
    ax2.set_ylabel("# alignments")
    ax2.set_title("Identity distribution: rank 1 vs rank 2 hits")
    ax2.legend()

    # ── Plot 3: Delta-identity (rank1 - rank2) histogram ──────────────────────
    ax3 = fig.add_subplot(gs[1, 1])
    r1 = top[top["rank"] == 1][["transcript", "identity"]].rename(columns={"identity": "id1"})
    r2 = top[top["rank"] == 2][["transcript", "identity"]].rename(columns={"identity": "id2"})
    merged = r1.merge(r2, on="transcript")
    merged["delta"] = merged["id1"] - merged["id2"]
    ax3.hist(merged["delta"], bins=50, color="#9C27B0", edgecolor="white", linewidth=0.3)
    ax3.axvline(0.05, color="red", linestyle="--", linewidth=1,
                label="δ=0.05 (gene conversion zone)")
    ax3.set_xlabel("Identity rank1 − rank2")
    ax3.set_ylabel("# transcripts")
    ax3.set_title("Identity gap between best and second-best V gene hit\n"
                  "(small gap → possible gene conversion)")
    ax3.legend(fontsize=8)

    # ── Plot 4: V gene coverage distribution (rank 1) ─────────────────────────
    ax4 = fig.add_subplot(gs[2, 0])
    ax4.hist(top[top["rank"] == 1]["v_coverage"], bins=50, color="#4CAF50",
             edgecolor="white", linewidth=0.3)
    ax4.axvline(1.0, color="red", linestyle="--", linewidth=1, label="100% coverage")
    ax4.set_xlabel("V gene coverage (fraction of gene aligned)")
    ax4.set_ylabel("# transcripts")
    ax4.set_title("Germline V gene coverage in best hit")
    ax4.legend()

    # ── Plot 5: Top-N identity per transcript (line per transcript sample) ─────
    ax5 = fig.add_subplot(gs[2, 1])
    sample_ids = top["transcript"].unique()
    rng = np.random.default_rng(42)
    sample_ids = rng.choice(sample_ids, size=min(200, len(sample_ids)), replace=False)
    for tid in sample_ids:
        sub = top[top["transcript"] == tid].sort_values("rank")
        ax5.plot(sub["rank"], sub["identity"], color="steelblue", alpha=0.07, linewidth=0.8)
    # Overlay median
    med = top.groupby("rank")["identity"].median()
    ax5.plot(med.index, med.values, color="red", linewidth=2, label="Median")
    ax5.set_xlabel("Alignment rank")
    ax5.set_ylabel("Identity")
    ax5.set_title(f"Identity decay across top-{top['rank'].max()} alignments\n"
                  "(sample of transcripts; red = median)")
    ax5.set_xticks(sorted(top["rank"].unique()))
    ax5.legend()

    fig.suptitle("IsoSeq IG Transcript V Gene Alignment Summary", fontsize=14, fontweight="bold")
    save_figure(fig, args.output, dpi=150)
    print(f"Saved {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
