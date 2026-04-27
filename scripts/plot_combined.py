"""
Combined cross-sample summary plot for the IsoSeq IG pipeline.

Panels
------
1. IGH V gene usage heatmap  — genes × samples, colour = transcript count
2. IGL V gene usage heatmap  — genes × samples, colour = transcript count
3. Locus assignment overview — stacked bar per sample (IGH-only / IGL-only /
                               both loci / exact-match only)
4. IGH × IGL co-occurrence   — per sample count of transcripts with hits to
                               BOTH loci (= strong true-antibody candidates)
5. IGH × IGL combination heatmap — for co-occurring transcripts, which IGH
                               gene pairs with which IGL gene (pooled across
                               all samples, top genes only)

Input
-----
--top-alns   : one top-N alignments TSV per sample (space-separated)
--exact-dirs : one exact_match_per_transcript TSV per sample (space-separated)
--samples    : matching sample labels (same order as --top-alns)
--output     : output PDF path
"""
import argparse
import sys
from pathlib import Path
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd


# ─── helpers ──────────────────────────────────────────────────────────────────

def load_top_alns(paths, sample_labels):
    frames = []
    for path, label in zip(paths, sample_labels):
        try:
            df = pd.read_csv(path, sep="\t")
            df["sample"] = label
            frames.append(df)
        except Exception as e:
            print(f"  Warning: could not load {path}: {e}", file=sys.stderr)
    if not frames:
        raise RuntimeError("No top-alignment files could be loaded.")
    return pd.concat(frames, ignore_index=True)


def load_exact(paths, sample_labels):
    frames = []
    for path, label in zip(paths, sample_labels):
        try:
            df = pd.read_csv(path, sep="\t")
            df["sample"] = label
            frames.append(df)
        except Exception as e:
            print(f"  Warning: could not load {path}: {e}", file=sys.stderr)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def load_exact_summaries(paths, sample_labels):
    frames = []
    for path, label in zip(paths, sample_labels):
        try:
            df = pd.read_csv(path, sep="\t")
            df["sample"] = label
            frames.append(df)
        except Exception as e:
            print(f"  Warning: could not load {path}: {e}", file=sys.stderr)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def plot_pooled_summary(exact_summary_all, top_all, output_path):
    """
    Reproduce the five per-sample plots from plot_summary.py on data pooled
    across all samples.
    """
    fig = plt.figure(figsize=(16, 14))
    gs = gridspec.GridSpec(3, 2, figure=fig, hspace=0.45, wspace=0.35)

    # ── Plot 1: Exact match counts per gene (summed across samples) ────────────
    ax1 = fig.add_subplot(gs[0, :])
    if not exact_summary_all.empty:
        pooled_exact = (
            exact_summary_all.groupby(["gene", "locus"], as_index=False)
            ["transcripts_with_exact_match"].sum()
        )
        pooled_exact = pooled_exact[pooled_exact["transcripts_with_exact_match"] > 0]
        pooled_exact = pooled_exact.sort_values("transcripts_with_exact_match", ascending=False)
    else:
        pooled_exact = pd.DataFrame()

    if not pooled_exact.empty:
        loci = pooled_exact["locus"].unique()
        palette = {l: c for l, c in zip(loci, plt.cm.Set2.colors)}
        colors = pooled_exact["locus"].map(palette)
        ax1.bar(range(len(pooled_exact)), pooled_exact["transcripts_with_exact_match"],
                color=colors, edgecolor="white", linewidth=0.4)
        ax1.set_xticks(range(len(pooled_exact)))
        ax1.set_xticklabels(pooled_exact["gene"], rotation=90, fontsize=6)
        ax1.set_ylabel("# IG transcripts with exact match (all samples)")
        ax1.set_title("Germline V genes found verbatim in transcripts — pooled across all samples")
        handles = [plt.Rectangle((0, 0), 1, 1, color=palette[l]) for l in loci]
        ax1.legend(handles, loci, title="Locus", loc="upper right")
    else:
        ax1.text(0.5, 0.5, "No exact matches found", ha="center", va="center",
                 transform=ax1.transAxes, fontsize=12)
        ax1.set_title("Germline V genes found verbatim in transcripts")

    # ── Plot 2: Identity distribution rank-1 vs rank-2 ────────────────────────
    ax2 = fig.add_subplot(gs[1, 0])
    rank1 = top_all[top_all["rank"] == 1]["identity"]
    rank2 = top_all[top_all["rank"] == 2]["identity"]
    bins = np.linspace(0, 1, 51)
    ax2.hist(rank1, bins=bins, alpha=0.7, label="Rank 1 (best)", color="#2196F3")
    ax2.hist(rank2, bins=bins, alpha=0.7, label="Rank 2", color="#FF9800")
    ax2.set_xlabel("Alignment identity")
    ax2.set_ylabel("# alignments")
    ax2.set_title("Identity distribution: rank 1 vs rank 2 hits (all samples)")
    ax2.legend()

    # ── Plot 3: Delta-identity histogram ──────────────────────────────────────
    ax3 = fig.add_subplot(gs[1, 1])
    r1 = (top_all[top_all["rank"] == 1][["sample", "transcript", "identity"]]
          .rename(columns={"identity": "id1"}))
    r2 = (top_all[top_all["rank"] == 2][["sample", "transcript", "identity"]]
          .rename(columns={"identity": "id2"}))
    merged = r1.merge(r2, on=["sample", "transcript"])
    merged["delta"] = merged["id1"] - merged["id2"]
    ax3.hist(merged["delta"], bins=50, color="#9C27B0", edgecolor="white", linewidth=0.3)
    ax3.axvline(0.05, color="red", linestyle="--", linewidth=1,
                label="δ=0.05 (gene conversion zone)")
    ax3.set_xlabel("Identity rank1 − rank2")
    ax3.set_ylabel("# transcripts")
    ax3.set_title("Identity gap between best and second-best V gene hit — all samples\n"
                  "(small gap → possible gene conversion)")
    ax3.legend(fontsize=8)

    # ── Plot 4: V gene coverage distribution (rank 1) ─────────────────────────
    ax4 = fig.add_subplot(gs[2, 0])
    ax4.hist(top_all[top_all["rank"] == 1]["v_coverage"], bins=50, color="#4CAF50",
             edgecolor="white", linewidth=0.3)
    ax4.axvline(1.0, color="red", linestyle="--", linewidth=1, label="100% coverage")
    ax4.set_xlabel("V gene coverage (fraction of gene aligned)")
    ax4.set_ylabel("# transcripts")
    ax4.set_title("Germline V gene coverage in best hit (all samples)")
    ax4.legend()

    # ── Plot 5: Top-N identity decay per transcript ───────────────────────────
    ax5 = fig.add_subplot(gs[2, 1])
    all_tids = top_all[["sample", "transcript"]].drop_duplicates()
    rng = np.random.default_rng(42)
    n_sample = min(200, len(all_tids))
    idx = rng.choice(len(all_tids), size=n_sample, replace=False)
    sampled = all_tids.iloc[idx]
    sampled_aln = sampled.merge(top_all, on=["sample", "transcript"])
    for (s, tid), grp in sampled_aln.groupby(["sample", "transcript"]):
        grp = grp.sort_values("rank")
        ax5.plot(grp["rank"], grp["identity"], color="steelblue", alpha=0.05, linewidth=0.8)
    med = top_all.groupby("rank")["identity"].median()
    ax5.plot(med.index, med.values, color="red", linewidth=2, label="Median")
    ax5.set_xlabel("Alignment rank")
    ax5.set_ylabel("Identity")
    ax5.set_title(f"Identity decay across top-{top_all['rank'].max()} alignments — all samples\n"
                  "(sample of transcripts; red = median)")
    ax5.set_xticks(sorted(top_all["rank"].unique()))
    ax5.legend()

    fig.suptitle("IsoSeq IG Transcript V Gene Alignment Summary — All Samples Pooled",
                 fontsize=13, fontweight="bold")
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"Saved {output_path}", file=sys.stderr)


def usage_heatmap(ax, counts_df, title, top_n=20):
    """
    counts_df : DataFrame with columns [sample, gene, count]
    Rows = top genes (sorted by total count), columns = samples.
    """
    pivot = counts_df.pivot_table(index="gene", columns="sample",
                                  values="count", aggfunc="sum", fill_value=0)
    # Keep top_n genes by total count across all samples
    pivot = pivot.loc[pivot.sum(axis=1).nlargest(top_n).index]
    pivot = pivot.sort_values(pivot.columns.tolist(), ascending=False)

    im = ax.imshow(pivot.values, aspect="auto", cmap="YlOrRd",
                   norm=mcolors.PowerNorm(gamma=0.5,
                                          vmin=0, vmax=pivot.values.max() or 1))
    ax.set_xticks(range(len(pivot.columns)))
    ax.set_xticklabels(pivot.columns, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(pivot.index)))
    ax.set_yticklabels(pivot.index, fontsize=7)
    ax.set_title(title, fontsize=10, fontweight="bold")
    plt.colorbar(im, ax=ax, label="# transcripts", pad=0.02)

    # Annotate cells with count (suppress zeros)
    for r in range(pivot.shape[0]):
        for c in pivot.shape[1] and range(pivot.shape[1]):
            val = pivot.values[r, c]
            if val > 0:
                ax.text(c, r, str(val), ha="center", va="center",
                        fontsize=6, color="black" if val < pivot.values.max() * 0.6 else "white")


def write_identity_tsv(aln, out_path):
    """
    Per-gene alignment identity statistics pooled across all samples and ranks.

    identity      = n_matches / aln_block   (fraction of aligned positions that match)
    v_coverage    = aln_block / target_len  (fraction of V gene that is aligned)
    full_exact    = identity == 1.0 AND v_coverage == 1.0  (verbatim full-gene match)

    identity == 1.0 with v_coverage < 1.0 means: a perfect-match *fragment* of
    the gene was found in the transcript — common with gene conversion.
    """
    rows = []
    for (gene, locus), grp in aln.groupby(["gene", "locus"]):
        n_total          = len(grp)
        n_id_1_0         = (grp["identity"] == 1.0).sum()
        n_id_gt_0_99     = (grp["identity"] >= 0.99).sum()
        n_id_gt_0_95     = (grp["identity"] >= 0.95).sum()
        n_id_gt_0_90     = (grp["identity"] >= 0.90).sum()
        n_full_exact     = ((grp["identity"] == 1.0) & (grp["v_coverage"] == 1.0)).sum()
        n_full_cov       = (grp["v_coverage"] >= 0.95).sum()
        n_id95_cov95     = ((grp["identity"] >= 0.95) & (grp["v_coverage"] >= 0.95)).sum()
        rows.append({
            "gene":                          gene,
            "locus":                         locus,
            "n_alignments_total":            n_total,
            "n_identity_1.0":                int(n_id_1_0),
            "n_identity_gte0.99":            int(n_id_gt_0_99),
            "n_identity_gte0.95":            int(n_id_gt_0_95),
            "n_identity_gte0.90":            int(n_id_gt_0_90),
            "n_full_exact_match":            int(n_full_exact),
            "n_full_gene_coverage":          int(n_full_cov),
            "n_identity_gte0.95_cov_gte0.95": int(n_id95_cov95),
            "max_identity":          round(float(grp["identity"].max()), 6),
            "median_identity":       round(float(grp["identity"].median()), 6),
            "mean_identity":         round(float(grp["identity"].mean()), 6),
            "max_v_coverage":        round(float(grp["v_coverage"].max()), 6),
            "median_v_coverage":     round(float(grp["v_coverage"].median()), 6),
        })
    df = (pd.DataFrame(rows)
            .sort_values(["n_identity_1.0", "n_identity_gte0.99"], ascending=False)
            .reset_index(drop=True))
    df.to_csv(out_path, sep="\t", index=False)
    print(f"Saved {out_path} ({len(df)} genes)", file=sys.stderr)
    return df


def plot_identity_analysis(aln, identity_df, out_path, top_n=25):
    """
    Four-panel figure clarifying why identity≈1.0 ≠ exact match.

    Panel 1 — 2D scatter: identity vs V gene coverage for all rank-1 alignments.
              The top-right corner (identity=1 & coverage=1) is the only zone that
              satisfies the strict exact-match criterion.  The top-left zone
              (identity=1, coverage<1) is a perfect-identity FRAGMENT — typical of
              gene conversion.

    Panel 2 — Zoomed identity histogram (0.85–1.0) with bars stacked by
              v_coverage band.  Separates "looks exact in the histogram" from
              "is actually a full-gene exact match".

    Panel 3 — Top genes by count of identity=1.0 alignments (any coverage).
              Shows which genes produce perfect-identity fragment hits.

    Panel 4 — Same but restricted to identity=1.0 AND v_coverage≥0.95.
              These are near-full-gene verbatim matches; if any gene reaches
              v_coverage=1.0 it will be the true exact match.
    """
    rank1 = aln[aln["rank"] == 1].copy()
    locus_palette = {"IGH": "#2196F3", "IGL": "#FF9800"}
    default_color = "#9E9E9E"

    # Dynamic height: fixed top section + equal-height rows for each bar panel
    scatter_h  = 5.0                          # inches for scatter + histogram row
    bar_h      = max(top_n * 0.38 + 1.5, 4)  # inches per bar panel (min 4)
    total_h    = scatter_h + 3 * bar_h + 1.5  # 1.5 for suptitle breathing room
    fig = plt.figure(figsize=(16, total_h))
    gs = gridspec.GridSpec(
        4, 2, figure=fig,
        height_ratios=[scatter_h, bar_h, bar_h, bar_h],
        hspace=0.55, wspace=0.35,
    )
    fig.suptitle("Alignment identity vs V gene coverage — all samples pooled",
                 fontsize=13, fontweight="bold")

    def _bar_panel(ax, data, value_col, xlabel, title, annotation_fn):
        """Shared logic for horizontal bar panels."""
        sub = data[data[value_col] > 0].nlargest(top_n, value_col)
        if sub.empty:
            ax.text(0.5, 0.5, f"No data for: {value_col}",
                    ha="center", va="center", transform=ax.transAxes)
            ax.set_title(title)
            return
        colors = [locus_palette.get(l, default_color) for l in sub["locus"]]
        ax.barh(range(len(sub)), sub[value_col], color=colors, edgecolor="white", height=0.7)
        ax.set_yticks(range(len(sub)))
        ax.set_yticklabels(sub["gene"], fontsize=7)
        ax.invert_yaxis()
        ax.set_xlabel(xlabel, fontsize=8)
        ax.set_title(title, fontsize=9)
        handles = [plt.Rectangle((0, 0), 1, 1, color=locus_palette.get(l, default_color))
                   for l in data["locus"].unique()]
        ax.legend(handles, data["locus"].unique(), fontsize=8, loc="lower right")
        xmax = sub[value_col].max()
        for i, (_, row) in enumerate(sub.iterrows()):
            label = annotation_fn(row)
            if label:
                ax.text(xmax * 0.01, i, label,
                        va="center", ha="left", fontsize=6, color="gray")

    # ── Panel 1: 2D scatter identity vs v_coverage ────────────────────────────
    ax = fig.add_subplot(gs[0, 0])
    for locus, grp in rank1.groupby("locus"):
        ax.scatter(grp["identity"], grp["v_coverage"],
                   alpha=0.15, s=4,
                   color=locus_palette.get(locus, default_color),
                   label=locus, rasterized=True)
    ax.axvline(1.0, color="red",   linestyle="--", linewidth=0.8, alpha=0.7)
    ax.axhline(1.0, color="red",   linestyle="--", linewidth=0.8, alpha=0.7)
    ax.axhline(0.95, color="gray", linestyle=":",  linewidth=0.7, alpha=0.6,
               label="95% coverage")
    ax.set_xlabel("Alignment identity  (n_matches / aln_block)")
    ax.set_ylabel("V gene coverage  (target_span / target_len)")
    ax.set_title("Identity vs V gene coverage\n(rank-1 alignments)")
    ax.legend(markerscale=3, fontsize=8)
    ax.text(1.001, 1.001, "exact\nzone", color="red", fontsize=6,
            ha="left", va="bottom", transform=ax.transAxes)

    # ── Panel 2: Zoomed identity histogram stacked by coverage band ───────────
    ax = fig.add_subplot(gs[0, 1])
    bins = np.linspace(0.85, 1.0, 61)
    cov_bands = [
        (0.95, 1.01, "#1B5E20", "coverage ≥ 95%  (near-full gene)"),
        (0.50, 0.95, "#81C784", "50% ≤ coverage < 95%"),
        (0.00, 0.50, "#E0E0E0", "coverage < 50%  (fragment)"),
    ]
    bottoms = np.zeros(len(bins) - 1)
    for lo, hi, color, label in cov_bands:
        mask = (rank1["identity"] >= 0.85) & \
               (rank1["v_coverage"] >= lo) & (rank1["v_coverage"] < hi)
        counts, _ = np.histogram(rank1.loc[mask, "identity"], bins=bins)
        ax.bar(bins[:-1], counts, width=np.diff(bins),
               bottom=bottoms, color=color, label=label,
               align="edge", edgecolor="none")
        bottoms += counts
    ax.axvline(1.0, color="red", linestyle="--", linewidth=1, label="identity = 1.0")
    ax.set_xlabel("Alignment identity")
    ax.set_ylabel("# rank-1 alignments")
    ax.set_title("Identity distribution (≥0.85), stacked by V gene coverage")
    ax.legend(fontsize=7, loc="upper left")
    ax.set_xlim(0.85, 1.005)

    # ── Panel 3: Top genes by identity = 1.0 (any coverage) ──────────────────
    _bar_panel(
        fig.add_subplot(gs[1, :]),
        identity_df, "n_identity_1.0",
        xlabel="# alignments with identity = 1.0  (any V gene coverage)",
        title=f"Top {top_n} genes by perfect-identity alignments  (includes partial-gene matches)",
        annotation_fn=lambda r: f"max_cov={r['max_v_coverage']:.2f}",
    )

    # ── Panel 4: Top genes by coverage ≥ 95% (any identity) ──────────────────
    _bar_panel(
        fig.add_subplot(gs[2, :]),
        identity_df, "n_full_gene_coverage",
        xlabel="# alignments with V gene coverage ≥ 95%  (any identity)",
        title=f"Top {top_n} genes by near-full-gene-coverage alignments",
        annotation_fn=lambda r: (f"exact={int(r['n_full_exact_match'])}"
                                  if r["n_full_exact_match"] > 0 else ""),
    )

    # ── Panel 5: identity ≥ 0.95 AND coverage ≥ 0.95 ─────────────────────────
    _bar_panel(
        fig.add_subplot(gs[3, :]),
        identity_df, "n_identity_gte0.95_cov_gte0.95",
        xlabel="# alignments with identity ≥ 0.95 AND V gene coverage ≥ 0.95",
        title=(f"Top {top_n} genes: identity ≥ 95% AND coverage ≥ 95%  "
               "(best combined signal for near-exact full-length matches)"),
        annotation_fn=lambda r: (f"exact={int(r['n_full_exact_match'])}"
                                  if r["n_full_exact_match"] > 0 else ""),
    )

    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Saved {out_path}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--top-alns", nargs="+", required=True,
                        help="top-N alignment TSV files, one per sample")
    parser.add_argument("--exact-per-transcript", nargs="+", required=True,
                        help="exact_match_per_transcript TSV files, one per sample")
    parser.add_argument("--samples", nargs="+", required=True,
                        help="Sample labels matching --top-alns order")
    parser.add_argument("--exact-summaries", nargs="+", required=True,
                        help="exact_match_summary TSV files, one per sample")
    parser.add_argument("--output", required=True)
    parser.add_argument("--out-pooled-summary", required=True,
                        help="PDF: same 5 plots as per-sample plot but pooled across all samples")
    parser.add_argument("--out-igh-usage", required=True,
                        help="TSV: total IGH gene usage across all samples")
    parser.add_argument("--out-igl-usage", required=True,
                        help="TSV: total IGL gene usage across all samples")
    parser.add_argument("--out-identity-plot", required=True,
                        help="PDF: alignment identity analysis figure")
    parser.add_argument("--out-identity-tsv", required=True,
                        help="TSV: per-gene alignment identity statistics")
    args = parser.parse_args()

    if len(args.top_alns) != len(args.samples):
        raise ValueError("--top-alns and --samples must have the same number of entries")

    print("Loading data...", file=sys.stderr)
    aln = load_top_alns(args.top_alns, args.samples)
    exact = load_exact(args.exact_per_transcript, args.samples)
    exact_summary = load_exact_summaries(args.exact_summaries, args.samples)

    # ── Derive per-transcript primary gene per locus ───────────────────────────
    # For each (sample, transcript, locus), find the highest-identity hit.
    # This is the "assigned" gene for that locus in that transcript.
    best_per_locus = (
        aln.sort_values("identity", ascending=False)
           .groupby(["sample", "transcript", "locus"], as_index=False)
           .first()[["sample", "transcript", "locus", "gene", "identity", "is_exact"]]
    )

    # ── Gene usage tables ──────────────────────────────────────────────────────
    def gene_counts(locus_name):
        sub = best_per_locus[best_per_locus["locus"] == locus_name]
        return (sub.groupby(["sample", "gene"])
                   .size()
                   .reset_index(name="count"))

    igh_counts = gene_counts("IGH")
    igl_counts = gene_counts("IGL")

    # ── Locus assignment per transcript ───────────────────────────────────────
    # Classify each transcript: IGH-only, IGL-only, both, or neither
    loci_per_transcript = (
        best_per_locus.groupby(["sample", "transcript"])["locus"]
        .apply(set)
        .reset_index()
        .rename(columns={"locus": "loci_set"})
    )
    loci_per_transcript["category"] = loci_per_transcript["loci_set"].map(
        lambda s: "Both (IGH + IGL)" if {"IGH", "IGL"} <= s
        else ("IGH only" if "IGH" in s
              else ("IGL only" if "IGL" in s else "Other"))
    )
    locus_summary = (
        loci_per_transcript.groupby(["sample", "category"])
        .size()
        .reset_index(name="count")
    )

    # ── Co-occurrence: transcripts with BOTH loci ──────────────────────────────
    both = loci_per_transcript[loci_per_transcript["category"] == "Both (IGH + IGL)"]
    cooccur_counts = both.groupby("sample").size().reset_index(name="count")

    # For co-occurring transcripts, find their IGH and IGL gene assignments
    cooccur_ids = both[["sample", "transcript"]]
    cooccur_genes = cooccur_ids.merge(best_per_locus, on=["sample", "transcript"])
    igh_in_cooccur = (cooccur_genes[cooccur_genes["locus"] == "IGH"]
                      [["sample", "transcript", "gene"]]
                      .rename(columns={"gene": "IGH_gene"}))
    igl_in_cooccur = (cooccur_genes[cooccur_genes["locus"] == "IGL"]
                      [["sample", "transcript", "gene"]]
                      .rename(columns={"gene": "IGL_gene"}))
    combinations = igh_in_cooccur.merge(igl_in_cooccur, on=["sample", "transcript"])
    combo_counts = (combinations.groupby(["IGH_gene", "IGL_gene"])
                    .size()
                    .reset_index(name="count"))

    # ── Figure layout ──────────────────────────────────────────────────────────
    fig = plt.figure(figsize=(20, 22))
    gs = gridspec.GridSpec(3, 2, figure=fig, hspace=0.55, wspace=0.45,
                           height_ratios=[1.4, 0.8, 1.4])

    SAMPLES_ORDERED = args.samples

    # Panel 1 — IGH gene usage heatmap
    ax1 = fig.add_subplot(gs[0, 0])
    if not igh_counts.empty:
        usage_heatmap(ax1, igh_counts, "IGH V gene usage across samples\n(top 20 genes)")
    else:
        ax1.text(0.5, 0.5, "No IGH alignments found", ha="center", va="center",
                 transform=ax1.transAxes)
        ax1.set_title("IGH V gene usage")

    # Panel 2 — IGL gene usage heatmap
    ax2 = fig.add_subplot(gs[0, 1])
    if not igl_counts.empty:
        usage_heatmap(ax2, igl_counts, "IGL V gene usage across samples\n(top 20 genes)")
    else:
        ax2.text(0.5, 0.5, "No IGL alignments found", ha="center", va="center",
                 transform=ax2.transAxes)
        ax2.set_title("IGL V gene usage")

    # Panel 3 — Locus assignment stacked bar
    ax3 = fig.add_subplot(gs[1, 0])
    cat_order = ["IGH only", "IGL only", "Both (IGH + IGL)", "Other"]
    cat_colors = {"IGH only": "#2196F3", "IGL only": "#FF9800",
                  "Both (IGH + IGL)": "#4CAF50", "Other": "#9E9E9E"}
    bottom = np.zeros(len(SAMPLES_ORDERED))
    sample_idx = {s: i for i, s in enumerate(SAMPLES_ORDERED)}
    for cat in cat_order:
        vals = []
        for s in SAMPLES_ORDERED:
            row = locus_summary[(locus_summary["sample"] == s) &
                                (locus_summary["category"] == cat)]
            vals.append(row["count"].values[0] if not row.empty else 0)
        ax3.bar(range(len(SAMPLES_ORDERED)), vals, bottom=bottom,
                label=cat, color=cat_colors[cat])
        bottom += np.array(vals)
    ax3.set_xticks(range(len(SAMPLES_ORDERED)))
    ax3.set_xticklabels(SAMPLES_ORDERED, rotation=45, ha="right", fontsize=8)
    ax3.set_ylabel("# IG transcripts")
    ax3.set_title("Locus assignment per sample", fontweight="bold")
    ax3.legend(loc="upper right", fontsize=8)

    # Panel 4 — Co-occurrence counts per sample
    ax4 = fig.add_subplot(gs[1, 1])
    if not cooccur_counts.empty:
        bar_vals = [cooccur_counts.set_index("sample")["count"].get(s, 0)
                    for s in SAMPLES_ORDERED]
        bars = ax4.bar(range(len(SAMPLES_ORDERED)), bar_vals, color="#4CAF50", edgecolor="white")
        ax4.bar_label(bars, padding=2, fontsize=8)
    ax4.set_xticks(range(len(SAMPLES_ORDERED)))
    ax4.set_xticklabels(SAMPLES_ORDERED, rotation=45, ha="right", fontsize=8)
    ax4.set_ylabel("# transcripts")
    ax4.set_title("Transcripts with BOTH IGH and IGL hits\n(true antibody transcript candidates)",
                  fontweight="bold")

    # Panel 5 — IGH × IGL combination heatmap (bottom row, full width)
    ax5 = fig.add_subplot(gs[2, :])
    if not combo_counts.empty:
        # Limit to top genes to keep the heatmap readable
        top_igh = (combo_counts.groupby("IGH_gene")["count"].sum()
                   .nlargest(15).index.tolist())
        top_igl = (combo_counts.groupby("IGL_gene")["count"].sum()
                   .nlargest(15).index.tolist())
        sub = combo_counts[combo_counts["IGH_gene"].isin(top_igh) &
                           combo_counts["IGL_gene"].isin(top_igl)]
        pivot = sub.pivot_table(index="IGH_gene", columns="IGL_gene",
                                values="count", aggfunc="sum", fill_value=0)
        # Order genes by marginal total
        pivot = pivot.loc[pivot.sum(axis=1).sort_values(ascending=False).index,
                          pivot.sum(axis=0).sort_values(ascending=False).index]

        im5 = ax5.imshow(pivot.values, aspect="auto", cmap="Blues",
                         norm=mcolors.PowerNorm(gamma=0.5,
                                                vmin=0, vmax=pivot.values.max() or 1))
        ax5.set_xticks(range(len(pivot.columns)))
        ax5.set_xticklabels(pivot.columns, rotation=60, ha="right", fontsize=8)
        ax5.set_yticks(range(len(pivot.index)))
        ax5.set_yticklabels(pivot.index, fontsize=8)
        ax5.set_xlabel("IGL gene", fontsize=9)
        ax5.set_ylabel("IGH gene", fontsize=9)
        plt.colorbar(im5, ax=ax5, label="# co-occurring transcripts", pad=0.01)
        for r in range(pivot.shape[0]):
            for c in range(pivot.shape[1]):
                val = pivot.values[r, c]
                if val > 0:
                    ax5.text(c, r, str(val), ha="center", va="center",
                             fontsize=7,
                             color="white" if val > pivot.values.max() * 0.6 else "black")
        ax5.set_title(
            "IGH × IGL gene combinations in co-occurring transcripts\n"
            "(pooled across all samples; top 15 genes per locus shown)",
            fontweight="bold"
        )
    else:
        ax5.text(0.5, 0.5,
                 "No transcripts found with both IGH and IGL hits",
                 ha="center", va="center", transform=ax5.transAxes, fontsize=12)
        ax5.set_title("IGH × IGL gene combinations", fontweight="bold")
        ax5.axis("off")

    fig.suptitle("Cross-sample IG transcript summary", fontsize=14, fontweight="bold", y=1.01)
    fig.savefig(args.output, dpi=150, bbox_inches="tight")
    print(f"Saved {args.output}", file=sys.stderr)

    # ── Pooled per-sample-style summary ───────────────────────────────────────
    plot_pooled_summary(exact_summary, aln, args.out_pooled_summary)

    # ── Gene usage TSVs ───────────────────────────────────────────────────────
    for locus_name, counts_df, out_path in (
        ("IGH", igh_counts, args.out_igh_usage),
        ("IGL", igl_counts, args.out_igl_usage),
    ):
        totals = (
            counts_df.groupby("gene")["count"]
            .sum()
            .reset_index(name="total_transcripts")
            .sort_values("total_transcripts", ascending=False)
            .reset_index(drop=True)
        )
        # Also include per-sample breakdown as extra columns
        per_sample = (
            counts_df.pivot_table(index="gene", columns="sample",
                                  values="count", aggfunc="sum", fill_value=0)
            .reset_index()
        )
        merged = totals.merge(per_sample, on="gene")
        merged.insert(0, "locus", locus_name)
        merged.to_csv(out_path, sep="\t", index=False)
        print(f"Saved {out_path} ({len(merged)} genes)", file=sys.stderr)

    # ── Identity analysis ─────────────────────────────────────────────────────
    identity_df = write_identity_tsv(aln, args.out_identity_tsv)
    plot_identity_analysis(aln, identity_df, args.out_identity_plot)

    # ── Print co-occurrence summary to stderr ──────────────────────────────────
    total_cooccur = len(both)
    total_transcripts = len(loci_per_transcript)
    print(
        f"\nCo-occurrence summary (pooled across all samples):\n"
        f"  Total IG transcripts : {total_transcripts}\n"
        f"  Both IGH + IGL hits  : {total_cooccur} ({total_cooccur/total_transcripts*100:.1f}%)",
        file=sys.stderr,
    )
    if not combinations.empty:
        print("\n  Top IGH×IGL combinations:", file=sys.stderr)
        for _, row in combo_counts.nlargest(10, "count").iterrows():
            print(f"    {row['IGH_gene']} × {row['IGL_gene']} : {row['count']}", file=sys.stderr)


if __name__ == "__main__":
    main()
