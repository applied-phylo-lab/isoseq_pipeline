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
        ax.set_title(f"{locus} — parent assignment changes",
                     fontsize=11, fontweight="bold", loc="left",
                     color=LOCUS.get(locus, INK))
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)

    ax = axes[-1]
    # Colour is the ONLY channel carrying locus.  Encoding locus with colour AND
    # shape, then overriding the colour on the control point, put a blue square
    # on the plot that matched no legend entry.  Shape now carries the one other
    # distinction that matters -- whether the reference is the same bird.
    for r in rows:
        same = r["same_bird"]
        ax.scatter(r["gain"], r["discord"],
                   s=150 if same else 90,
                   marker="D" if same else "o",
                   color=LOCUS.get(r["locus"], INK),
                   edgecolor="black", lw=1.3 if same else .6, zorder=3)
    # Several references land on top of each other at ~0% identity gain, which is
    # the substantive point of the panel -- so the labels have to be fanned out
    # rather than left overlapping.
    # Displaced labels go DOWN AND LEFT.  Pushing them down-right walks them onto
    # the next marker in the cluster, which is what hid "b0_pri" behind "b1_pri".
    used = []
    for r in sorted(rows, key=lambda r: (-r["gain"], -r["discord"])):
        dx, dy, ha, step = 6, 5, "left", 0
        while any(abs(r["gain"] - gx) < 0.13
                  and abs(r["discord"] + dy / 3.0 - gy) < 2.6
                  for gx, gy in used):
            step += 1
            dy -= 11
            dx, ha = (-7, "right") if step % 2 else (6, "left")
        used.append((r["gain"], r["discord"] + dy / 3.0))
        ax.annotate(r["ref"].replace("bAgePho", "b"),
                    (r["gain"], r["discord"]), fontsize=6.5, ha=ha,
                    xytext=(dx, dy), textcoords="offset points")
    ax.axhline(0, color=GREY, lw=1)
    ax.axvline(0, color=GREY, lw=1)
    # Most references cluster at ~0% identity gain, so without a left margin a
    # displaced label lands on top of the y-axis tick labels.
    gains = [r["gain"] for r in rows]
    span = max(max(gains) - min(gains), 0.1)
    ax.set_xlim(min(gains) - 0.16 * span, max(gains) + 0.12 * span)
    ax.set_xlabel("median identity gained with matched germline (%)", fontsize=9)
    ax.set_ylabel("transcripts given a different parent (%)", fontsize=9)

    def handle(marker, colour, label, size=7, lw=.6):
        return Line2D([], [], marker=marker, linestyle="none", markersize=size,
                      color=colour, markeredgecolor="black",
                      markeredgewidth=lw, label=label)

    locus_leg = ax.legend(
        handles=[handle("o", LOCUS.get(l, INK), l) for l in loci],
        fontsize=8, title="locus", title_fontsize=8, loc="upper left")
    ax.add_artist(locus_leg)
    ax.legend(handles=[handle("o", GREY_DARK, "different bird"),
                       handle("D", GREY_DARK, "same bird, other haplotype",
                              size=8, lw=1.3)],
              fontsize=8, title="reference", title_fontsize=8, loc="lower right")
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
