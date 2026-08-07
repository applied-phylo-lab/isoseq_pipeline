"""
Summarise, across several reference haplotypes, what changes when the germline
comes from the wrong bird.

The design that makes this interpretable
----------------------------------------
The same transcripts are scored against every reference, so all differences come
from the reference alone.  The references split into two kinds:

  * a DIFFERENT individual (bAgePho0, bAgePho1) -- what you get when you borrow
    a published assembly.
  * the SAME individual's other haplotype (bAgePho2_alt).

What the same-bird row is, and is not
-------------------------------------
It is NOT a technical noise floor, and must not be described as one.  The bird is
diploid and both haplotypes are transcribed, so a V gene on the alt haplotype is a
genuinely expressible gene with its own alleles.  When a transcript is assigned to
an alt gene rather than the pri orthologue, that assignment can be *correct* --
the transcript may really have come from the alt allele.  So the discordance in
this row mixes three things that cannot be separated here:

  * real allelic origin: the transcript came from the other haplotype's copy;
  * genuine copy-number and content differences between the two haplotypes;
  * assembly differences (fragmentation, collapsed or duplicated genes).

What it therefore measures is the floor imposed by scoring a DIPLOID animal
against a HAPLOID reference -- the irreducible cost of picking one haplotype,
which applies to the matched reference too.  The matched run is likewise only
half the bird's germline.

Why it is still the right control
---------------------------------
Every different-bird reference carries that same haploid-reference cost, and adds
between-individual divergence on top.  So the excess over the same-bird row still
isolates the part attributable to using a different animal, which is the claim
being made.  If anything the baseline is generous: some of the 21.5% is real
biology rather than error, so the true technical floor is lower and the different-
bird penalty correspondingly larger than the subtraction suggests.
"""
import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

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
    fig, axes = plt.subplots(1, len(loci), figsize=(5.6 * len(loci), 4.6),
                             squeeze=False)
    axes = axes[0]

    titles = []
    for panel, (ax, locus) in enumerate(zip(axes, loci)):
        sub = [r for r in rows if r["locus"] == locus]
        labels = [r["ref"] for r in sub]
        vals = [r["discord"] for r in sub]
        cols = [YES if r["same_bird"] else NO for r in sub]
        ax.barh(range(len(sub)), vals, color=cols, edgecolor="black", lw=.6)
        ax.set_yticks(range(len(sub)))
        ax.set_yticklabels(
            [f"{l}{'  (same bird, other haplotype)' if r['same_bird'] else ''}"
             for l, r in zip(labels, sub)], fontsize=8.5)
        for i, r in enumerate(sub):
            ax.text(r["discord"], i, f"  {r['discord']:.0f}%  (n={r['n']})",
                    va="center", fontsize=8)
        base = next((r["discord"] for r in sub if r["same_bird"]), None)
        if base is not None:
            ax.axvline(base, color=INK, ls="--", lw=1.4, zorder=5)
            # Deliberately not called a "noise floor": both haplotypes are
            # transcribed, so part of this bar is a transcript genuinely coming
            # from the other allele, not an error.
            ax.annotate(" haploid-reference floor (same bird)",
                        xy=(base, 1.0), xycoords=("data", "axes fraction"),
                        xytext=(3, -4), textcoords="offset points",
                        fontsize=7.5, color=INK, va="top", ha="left")
        ax.set_xlabel("transcripts given a DIFFERENT parent (%)", fontsize=9)
        ax.set_xlim(0, max(vals) * 1.45)
        # Letter bold, the rest plain, and black rather than the locus colour:
        # the locus is already named in the title, so colouring it too makes the
        # heading carry data it does not have.
        letter = chr(ord("A") + panel)
        titles.append(r"$\bf{" + letter + r"}$" +
                      f"  {locus} — parent assignment changes")
        ax.set_title(titles[-1], fontsize=11, loc="left", color="black")
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)

    fig.tight_layout()

    # Push each title out to the visual left edge of its panel. loc="left" aligns
    # to the AXES BOX, but the reference names on the y-axis are long (one runs to
    # "bAgePho2_alt  (same bird, other haplotype)") and sit outside that box, so
    # the heading ends up looking centred over the panel. Measuring the rendered
    # tick labels and aligning to the leftmost of them handles both panels without
    # hard-coding an offset per panel.
    # set_title(x=...) sticks; Text.set_x() alone does not, because matplotlib
    # re-applies the loc-based position on every draw.
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    for ax, title in zip(axes, titles):
        labs = [t for t in ax.get_yticklabels() if t.get_text()]
        if not labs:
            continue
        box = ax.get_window_extent(rend)
        left = min(t.get_window_extent(rend).x0 for t in labs)
        ax.set_title(title, fontsize=11, loc="left", color="black",
                     x=(left - box.x0) / box.width)

    save_figure(fig, args.out_figure)

    print(f"{len(rows)} comparisons written to {args.out_table}", file=sys.stderr)
    for r in rows:
        tag = " (SAME BIRD)" if r["same_bird"] else ""
        print(f"  {r['ref']:14s} {r['locus']}  discordance {r['discord']:5.1f}%  "
              f"identity gain {r['gain']:+.2f}%{tag}", file=sys.stderr)


if __name__ == "__main__":
    main()
