"""
Summarise, across several reference haplotypes, what changes when the germline
comes from the wrong bird.

The design that makes this interpretable
----------------------------------------
The same transcripts are scored against every reference, so all differences come
from the reference alone.  The references split into two kinds:

  * a DIFFERENT individual (bAgePho0, bAgePho1) -- what you get when you borrow
    a published assembly.
  * the SAME individual's other haplotype (bAgePho2_alt) -- same bird, so any
    discordance here is assembly and haplotype noise, not individual mismatch.

That second one is the control the analysis needs.  Without it, a discordance
rate has no baseline: you cannot tell how much is "wrong bird" and how much is
simply "two assemblies of anything never agree perfectly".  With it, the excess
over the same-bird row is the part attributable to using a different animal.
"""
import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from gc_palette import save_figure, LOCUS, NO, YES, GREY, GREY_DARK, INK


def read_kv(path):
    d = {}
    with open(path) as fh:
        fh.readline()
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2:
                d[p[0]] = p[1]
    return d


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stats", nargs="+", required=True,
                    help="cmp_<reference>_<locus>.tsv files")
    ap.add_argument("--same-bird", default="bAgePho2_alt",
                    help="reference that is the SAME individual as the RNA "
                         "(its other haplotype); used as the baseline")
    ap.add_argument("--matched-label", default="bAgePho2_pri")
    ap.add_argument("--out-figure", required=True)
    ap.add_argument("--out-table", required=True)
    args = ap.parse_args()

    rows = []
    for path in args.stats:
        base = os.path.basename(path).replace("cmp_", "").replace(".tsv", "")
        ref, locus = base.rsplit("_", 1)
        d = read_kv(path)
        try:
            rows.append({
                "ref": ref, "locus": locus,
                "n": int(d["transcripts_compared"]),
                "discord": float(d["parent_discordance_rate"]) * 100,
                "gain": float(d["median_identity_gain"]) * 100,
                "used_a": int(d["genes_called_used_a"]),
                "used_b": int(d["genes_called_used_b"]),
                "same_bird": ref == args.same_bird,
            })
        except KeyError:
            print(f"skipping {path}: missing fields", file=sys.stderr)

    if not rows:
        raise SystemExit("no usable stats files")
    rows.sort(key=lambda r: (r["locus"], r["same_bird"], -r["discord"]))

    with open(args.out_table, "w") as fh:
        fh.write("reference\tlocus\tsame_individual\ttranscripts\t"
                 "parent_discordance_pct\tmedian_identity_gain_pct\t"
                 "genes_look_expressed_ref\tgenes_look_expressed_matched\n")
        for r in rows:
            fh.write(f"{r['ref']}\t{r['locus']}\t{r['same_bird']}\t{r['n']}\t"
                     f"{r['discord']:.1f}\t{r['gain']:.2f}\t"
                     f"{r['used_a']}\t{r['used_b']}\n")

    loci = sorted({r["locus"] for r in rows})
    fig, axes = plt.subplots(1, len(loci) + 1,
                             figsize=(5.2 * (len(loci) + 1), 4.6))
    if len(loci) + 1 == 1:
        axes = [axes]

    for ax, locus in zip(axes, loci):
        sub = [r for r in rows if r["locus"] == locus]
        labels = [r["ref"] for r in sub]
        vals = [r["discord"] for r in sub]
        cols = [YES if r["same_bird"] else NO for r in sub]
        ax.barh(range(len(sub)), vals, color=cols, edgecolor="black", lw=.6)
        ax.set_yticks(range(len(sub)))
        ax.set_yticklabels([f"{l}{'  (same bird)' if r['same_bird'] else ''}"
                            for l, r in zip(labels, sub)], fontsize=8.5)
        for i, r in enumerate(sub):
            ax.text(r["discord"], i, f"  {r['discord']:.0f}%  (n={r['n']})",
                    va="center", fontsize=8)
        base = next((r["discord"] for r in sub if r["same_bird"]), None)
        if base is not None:
            ax.axvline(base, color=INK, ls="--", lw=1.4, zorder=5)
            ax.text(base, len(sub) - 0.35,
                    " same-bird baseline", fontsize=7.5, color=INK, va="top")
        ax.set_xlabel("transcripts given a DIFFERENT parent (%)", fontsize=9)
        ax.set_xlim(0, max(vals) * 1.45)
        ax.set_title(f"{locus} — parent assignment changes",
                     fontsize=11, fontweight="bold", loc="left",
                     color=LOCUS.get(locus, INK))
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)

    ax = axes[-1]
    for locus, mark in zip(loci, ("o", "s")):
        sub = [r for r in rows if r["locus"] == locus]
        ax.scatter([r["gain"] for r in sub], [r["discord"] for r in sub],
                   s=90, marker=mark, edgecolor="black", lw=.6,
                   color=[YES if r["same_bird"] else LOCUS.get(locus, NO) for r in sub],
                   label=locus, zorder=3)
        for r in sub:
            ax.annotate(r["ref"].replace("bAgePho", "b"),
                        (r["gain"], r["discord"]), fontsize=6.5,
                        xytext=(4, 4), textcoords="offset points")
    ax.axhline(0, color=GREY, lw=1)
    ax.axvline(0, color=GREY, lw=1)
    ax.set_xlabel("median identity gained with the matched germline (%)", fontsize=9)
    ax.set_ylabel("transcripts given a different parent (%)", fontsize=9)
    ax.legend(fontsize=8, title="locus", title_fontsize=8)
    ax.set_title("identity barely moves,\nparent assignment does",
                 fontsize=11, fontweight="bold", loc="left")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    fig.suptitle("Effect of the germline reference on what you infer "
                 f"(all references scored against {args.matched_label}, "
                 "identical transcripts throughout)",
                 fontsize=12.5, fontweight="bold", y=1.03)
    fig.tight_layout()
    save_figure(fig, args.out_figure)

    print(f"{len(rows)} comparisons written to {args.out_table}", file=sys.stderr)
    for r in rows:
        tag = " (SAME BIRD)" if r["same_bird"] else ""
        print(f"  {r['ref']:14s} {r['locus']}  discordance {r['discord']:5.1f}%  "
              f"identity gain {r['gain']:+.2f}%{tag}", file=sys.stderr)


if __name__ == "__main__":
    main()
