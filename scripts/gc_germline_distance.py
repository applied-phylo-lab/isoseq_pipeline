"""
How far every transcript sits from the closest germline V gene.

Why this figure exists, and why it is not the figure that was asked for
----------------------------------------------------------------------
The question was "show me the transcripts that match their parent perfectly, and
which genes they belong to".  The answer is that there are none: not one of the
595 transcripts matches ANY germline V gene in the assembly exactly, and the
closest sits 1 mismatch away.  A figure of the perfect matches would be an empty
panel, so this shows the distribution instead, with the closest cases named.

Distance is to the closest germline gene of ANY kind, NOT to the RSS-restricted
parent.  The figure is used to ask whether the RSS annotation is missing genes,
and restricting the comparison to RSS-bearing genes would assume the answer.

That null result is worth more than the figure would have been.  Two things
follow from it.

  * There is no unmutated compartment in this library.  Every transcript has
    diverged from every germline V gene in the assembly, which is what a
    bursa-diversified repertoire should look like -- peripheral B cells in birds
    have all been through gene conversion before they are ever sampled.
  * It bounds what any assignment can be.  If the nearest germline is ~11-14
    mismatches away, the difference between the best and second-best gene is
    always going to be a handful of bases, so the assignment step is choosing
    between two distant options rather than recognising one exact source.  That
    is the same fact the margin analysis reports from the other direction.

The distance shown is over the assessable region only (junction-trimmed), so it
is mutations in V, not junctional diversity at the 3' end.
"""
import argparse
import csv
import statistics as st
import sys
from collections import Counter, defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

from gc_palette import save_figure, LOCUS, SEGMENT, GREY, GREY_DARK, INK


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def short(name):
    parts = name.split(".")
    return parts[1] if len(parts) > 1 else name


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--assignments", required=True,
                    help="transcript_assignments.tsv (parent, RSS-restricted)")
    ap.add_argument("--unconstrained",
                    help="unconstrained_assignments.tsv (closest gene of any kind)")
    ap.add_argument("--rss-annotation", required=True)
    ap.add_argument("--out-figure", required=True)
    ap.add_argument("--out-table", required=True)
    args = ap.parse_args()

    par = read_tsv(args.assignments)
    unc = read_tsv(args.unconstrained) if args.unconstrained else []
    rss = {r["gene"]: r for r in read_tsv(args.rss_annotation)}
    loci = sorted({r["locus"] for r in par})

    with open(args.out_table, "w") as fh:
        fh.write("locus\tset\tn_transcripts\tn_exact\tmin_mismatches\t"
                 "median_mismatches\tmedian_identity\tclosest_gene\n")
        for tag, rows in (("parent", par), ("closest_any_gene", unc)):
            for L in loci:
                sub = [r for r in rows if r["locus"] == L]
                if not sub:
                    continue
                mm = [int(r["mismatches"]) for r in sub]
                best = min(sub, key=lambda r: int(r["mismatches"]))
                fh.write(f"{L}\t{tag}\t{len(sub)}\t{sum(1 for m in mm if m == 0)}\t"
                         f"{min(mm)}\t{st.median(mm):.0f}\t"
                         f"{st.median(float(r['identity']) for r in sub):.4f}\t"
                         f"{best['best_gene']}\n")

    fig, axes = plt.subplots(1, 2, figsize=(12.4, 4.6),
                             gridspec_kw={"width_ratios": [1.15, 1.0]})

    # ── A: how far from germline ──────────────────────────────────────────────
    ax = axes[0]
    # Distance to the CLOSEST germline gene of any kind -- deliberately not the
    # RSS-restricted parent. The question this figure exists to answer is whether
    # the RSS annotation is missing genes, and restricting to RSS-bearing genes
    # would build the answer into the measurement.
    #
    # Bin width is exactly 1, centred on integers. A wider bin starting at 0 puts
    # the 1-mismatch transcripts in a bar spanning [0, 2), which reads as though
    # something sits at zero. Nothing does, and the figure has to make that
    # unambiguous rather than nearly-unambiguous.
    src = unc if unc else par
    hi = max(int(r["mismatches"]) for r in src)
    bins = np.arange(-0.5, hi + 1.5, 1.0)
    for L in loci:
        cg = [int(r["mismatches"]) for r in src if r["locus"] == L]
        ax.hist(cg, bins=bins, histtype="step", lw=1.8,
                color=LOCUS.get(L, INK), label=f"{L}  (n={len(cg)})")
    ax.axvline(0, color=GREY_DARK, lw=1.0, ls=":", zorder=1)
    ax.set_xlabel("mismatches to the closest germline V gene of any kind "
                  "(assessable region)", fontsize=9)
    ax.set_ylabel("transcripts", fontsize=9)
    ax.legend(fontsize=8.5, frameon=False, loc="upper right")
    ax.set_title(r"$\bf{A}$  every transcript has diverged from every germline V",
                 fontsize=10.5, loc="left", x=0.0)
    ax.set_xlim(-1.2, hi + 1.2)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    # ── B: the closest transcript per gene ────────────────────────────────────
    # Per gene rather than per transcript: the question was which GENES the
    # near-germline transcripts belong to, and a gene with one near-germline
    # transcript is the interesting object even if its other 40 are far away.
    #
    # Same unconstrained measure as panel A, and the same one the main figure's
    # architecture panels use to count genes, so the three agree on what a gene
    # being "used" means. RSS state is marked rather than filtered on: a gene
    # whose closest transcript sits 1 mismatch away but which carries no RSS is
    # exactly the object this figure is meant to surface.
    ax = axes[1]
    closest = defaultdict(lambda: 10 ** 9)
    count = Counter()
    locus_of = {}
    for r in src:
        g = r["best_gene"]
        closest[g] = min(closest[g], int(r["mismatches"]))
        count[g] += 1
        locus_of[g] = r["locus"]
    # Selected by closeness -- these are the genes the figure is about -- but then
    # ORDERED by locus and contig position, which is how every other figure in the
    # set lays genes out. Leaving them in rank order made the bars monotonic,
    # which reads as a trend and is really just the sort.
    def pos_of(g):
        try:
            return int(rss.get(g, {}).get("pos", 0))
        except ValueError:
            return 0
    order = sorted(closest, key=lambda g: (closest[g], -count[g]))[:20]
    order.sort(key=lambda g: (locus_of[g], pos_of(g)))
    y = np.arange(len(order))
    ax.barh(y, [closest[g] for g in order],
            color=[LOCUS.get(locus_of[g], INK) for g in order],
            edgecolor="black", lw=0.5)
    for i, g in enumerate(order):
        ax.text(closest[g] + 0.4, i, f" n={count[g]}", va="center", fontsize=7.4,
                color=GREY_DARK)
        if rss.get(g, {}).get("rss_state") == "rss_present":
            ax.plot([-0.5], [i], marker="o", ms=4.6, color=SEGMENT["V"],
                    markeredgecolor="black", markeredgewidth=0.5,
                    clip_on=False, zorder=6)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{locus_of[g]} {short(g)}" for g in order], fontsize=7.6)
    ax.invert_yaxis()
    ax.set_xlabel("mismatches of that gene's CLOSEST transcript "
                  "(same measure as A)", fontsize=9)
    ax.set_xlim(-1.0, max(closest[g] for g in order) * 1.22)
    n_g = len(closest)
    ax.set_title(r"$\bf{B}$  the 20 genes closest to germline",
                 fontsize=10.5, loc="left", x=0.0)
    ax.legend(handles=[Line2D([], [], marker="o", ls="none", ms=5.5,
                              color=SEGMENT["V"], markeredgecolor="black",
                              markeredgewidth=0.5, label="has an RSS")],
              fontsize=8, frameon=False, loc="lower right")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    n_tot = len(src)
    n_exact = sum(1 for r in src if int(r["mismatches"]) == 0)
    mn = min(int(r["mismatches"]) for r in src)
    fig.tight_layout()
    save_figure(fig, args.out_figure)

    print(f"{n_exact}/{n_tot} transcripts match any germline V exactly; "
          f"closest {mn} mismatches", file=sys.stderr)
    for tag, rows in (("parent", par), ("closest gene", unc or par)):
        for L in loci:
            mm = [int(r["mismatches"]) for r in rows if r["locus"] == L]
            print(f"  {tag:13s} {L}: n={len(mm)} min={min(mm)} "
                  f"median={st.median(mm):.0f} max={max(mm)}", file=sys.stderr)


if __name__ == "__main__":
    main()
