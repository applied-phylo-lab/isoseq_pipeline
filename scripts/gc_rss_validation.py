"""
Is the RSS annotation any good?  Test it against expression, which knows nothing
about it.

The logic
---------
An RSS is what lets a V gene be rearranged.  A gene without one cannot be the
parent of a transcript -- it can only ever donate sequence into somebody else's
rearranged gene.  So if the RSS calls are right, expression should be lopsided:

    RSS present + expressed      expected, the working genes
    RSS present + silent         fine: functional but not used in this bird
    RSS absent  + silent         expected, the donor array
    RSS absent  + EXPRESSED      the problem cell -- either the RSS was missed,
                                 or those transcripts are misassigned

The last cell is the one that measures annotation quality, and it is why this is
worth plotting rather than tabulating: a gene there is either a false negative in
the RSS screen or a donor so similar to the real parent that transcripts land on
it by mistake.  Those two explanations look different in the data, and panel C
separates them -- a missed-RSS gene should carry transcripts that match it
BETTER than they match any RSS-bearing gene, while a similarity artefact should
not.

Expression is deliberately taken from the UNCONSTRAINED assignment (every V gene
allowed to receive transcripts).  Using the constrained one would be circular:
it only lets RSS genes be parents, so the interesting cell would be empty by
construction.
"""
import argparse
import statistics
import sys
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

from gc_palette import (save_figure, LOCUS, SEGMENT, GREY, GREY_DARK, INK,
                        YES, NO, BOTH)


def read_tsv(path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(header, line.rstrip("\n").split("\t")))
                for line in fh if line.strip()]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rss-annotation", required=True)
    ap.add_argument("--assignments", required=True,
                    help="unconstrained_assignments.tsv -- every gene eligible")
    ap.add_argument("--constrained-assignments",
                    help="transcript_assignments.tsv -- best RSS-bearing parent "
                         "per transcript; enables panel C")
    ap.add_argument("--min-identity", type=float, default=0.98,
                    help="identity for a transcript to count as expressing a gene")
    ap.add_argument("--out-figure", required=True)
    ap.add_argument("--out-table", required=True)
    args = ap.parse_args()

    rss = {r["gene"]: r for r in read_tsv(args.rss_annotation)}

    # Expression = the gene is some transcript's best hit.  Deliberately NOT
    # gated on the 0.98 identity bar: in a gene-conversion locus every transcript
    # is a mosaic sitting ~5% from any single germline gene, so that bar would
    # discard almost everything and leave the panels empty.  The strict count is
    # kept alongside as the "near-perfect match" evidence.
    n_tx = defaultdict(int)
    n_tx_strict = defaultdict(int)
    best_id = defaultdict(float)
    for r in read_tsv(args.assignments):
        g = r.get("best_gene") or r.get("gene")
        if not g:
            continue
        try:
            idt = float(r.get("identity", r.get("best_identity", 0)) or 0)
        except ValueError:
            idt = 0.0
        best_id[g] = max(best_id[g], idt)
        n_tx[g] += 1
        if idt >= args.min_identity:
            n_tx_strict[g] += 1

    loci = sorted({r["locus"] for r in rss.values()})
    rows = []
    for g, r in rss.items():
        has = r["rss_state"] == "rss_present"
        expressed = n_tx[g] > 0
        rows.append({"gene": g, "locus": r["locus"], "pos": int(r["pos"]),
                     "rss": has, "expressed": expressed,
                     "n_tx": n_tx[g], "n_tx_strict": n_tx_strict[g],
                     "best_id": best_id[g],
                     "cell": ("RSS+/expr+" if has and expressed else
                              "RSS+/expr-" if has else
                              "RSS-/expr+" if expressed else "RSS-/expr-")})

    with open(args.out_table, "w") as fh:
        fh.write("gene\tlocus\tpos\trss_present\texpressed\tn_transcripts\t"
                 "n_transcripts_near_perfect\tbest_identity\tcell\n")
        for r in sorted(rows, key=lambda r: (r["locus"], -r["n_tx"])):
            fh.write(f"{r['gene']}\t{r['locus']}\t{r['pos']}\t{r['rss']}\t"
                     f"{r['expressed']}\t{r['n_tx']}\t{r['n_tx_strict']}\t"
                     f"{r['best_id']:.4f}\t"
                     f"{r['cell']}\n")

    # ── figure ────────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.8),
                             gridspec_kw={"width_ratios": [1.15, 1.35, 1.3]})

    CELLS = ["RSS+/expr+", "RSS+/expr-", "RSS-/expr+", "RSS-/expr-"]
    CELL_C = {"RSS+/expr+": YES, "RSS+/expr-": BOTH,
              "RSS-/expr+": NO, "RSS-/expr-": GREY}

    # A: the 2x2 itself, per locus
    ax = axes[0]
    width = 0.38
    for i, locus in enumerate(loci):
        sub = [r for r in rows if r["locus"] == locus]
        counts = [sum(1 for r in sub if r["cell"] == c) for c in CELLS]
        xs = np.arange(len(CELLS)) + (i - (len(loci) - 1) / 2) * width
        ax.bar(xs, counts, width=width * 0.92,
               color=LOCUS.get(locus, INK), edgecolor="black", lw=.6,
               label=f"{locus}  (n={len(sub)})")
        for x, c in zip(xs, counts):
            if c:
                ax.text(x, c, f"{c}", ha="center", va="bottom", fontsize=8)
    ax.set_xticks(range(len(CELLS)))
    ax.set_xticklabels(["RSS\n+expressed", "RSS\nsilent",
                        "no RSS\n+EXPRESSED", "no RSS\nsilent"], fontsize=8)
    ax.set_ylabel("V genes")
    ax.legend(fontsize=8)
    ax.set_title("A  RSS state vs expression\nthird bar is the annotation's error budget",
                 fontsize=10, fontweight="bold", loc="left")

    # B: how much expression does each cell actually carry?
    ax = axes[1]
    for i, locus in enumerate(loci):
        sub = [r for r in rows if r["locus"] == locus]
        tot = sum(r["n_tx"] for r in sub) or 1
        share = [sum(r["n_tx"] for r in sub if r["cell"] == c) / tot * 100
                 for c in CELLS]
        xs = np.arange(len(CELLS)) + (i - (len(loci) - 1) / 2) * width
        ax.bar(xs, share, width=width * 0.92, color=LOCUS.get(locus, INK),
               edgecolor="black", lw=.6, label=locus)
        for x, s in zip(xs, share):
            if s > 0.5:
                ax.text(x, s, f"{s:.0f}%", ha="center", va="bottom", fontsize=8)
    ax.set_xticks(range(len(CELLS)))
    ax.set_xticklabels(["RSS\n+expressed", "RSS\nsilent",
                        "no RSS\n+EXPRESSED", "no RSS\nsilent"], fontsize=8)
    ax.set_ylabel("% of transcripts")
    ax.legend(fontsize=8)
    ax.set_title("B  where the transcripts actually go\n"
                 "most land on a gene with no RSS — panel C asks why",
                 fontsize=10, fontweight="bold", loc="left")

    # C: for every transcript that prefers a no-RSS gene, how much better is that
    # gene than the best gene that CAN actually rearrange?
    #
    #   large margin  -> the no-RSS gene really is the better explanation, which
    #                    points at a missed RSS (or a missing allele);
    #   margin ~ 0    -> the two are interchangeable, so preferring the donor is
    #                    a coin flip, exactly what heavy conversion produces when
    #                    a transcript drifts toward the sequence it copied.
    #
    # This is the discriminating measurement; counting genes cannot separate them.
    ax = axes[2]
    margins = defaultdict(list)
    if args.constrained_assignments:
        constrained = {r["transcript"]: r
                       for r in read_tsv(args.constrained_assignments)}
        for r in read_tsv(args.assignments):
            g = r["best_gene"]
            if g not in rss or rss[g]["rss_state"] == "rss_present":
                continue
            c = constrained.get(r["transcript"])
            if not c or not c.get("best_gene"):
                continue
            try:
                margins[r["locus"]].append(
                    (float(r["identity"]) - float(c["identity"])) * 100)
            except ValueError:
                continue

    if any(margins.values()):
        bins = np.arange(0, 8.25, 0.25)
        for locus in loci:
            if not margins.get(locus):
                continue
            ax.hist(margins[locus], bins=bins, alpha=.72,
                    color=LOCUS.get(locus, INK), edgecolor="black", lw=.4,
                    label=f"{locus}  (n={len(margins[locus])}, "
                          f"median {statistics.median(margins[locus]):.2f}%)")
        ax.axvline(1.0, color=INK, ls="--", lw=1.5)
        ax.annotate("← interchangeable:\n    conversion artefact",
                    xy=(1.0, 0.02), xycoords=("data", "axes fraction"),
                    xytext=(-6, 0), textcoords="offset points",
                    ha="right", va="bottom", fontsize=7.5, color=INK)
        ax.annotate("clearly better:\nlikely MISSED RSS →",
                    xy=(1.0, 0.02), xycoords=("data", "axes fraction"),
                    xytext=(6, 0), textcoords="offset points",
                    ha="left", va="bottom", fontsize=7.5, color=NO)
        ax.set_xlabel("how much better the no-RSS gene fits than the best\n"
                      "RSS-bearing gene (percentage points of identity)",
                      fontsize=8.5)
        ax.set_ylabel("transcripts")
        ax.legend(fontsize=7.5)
        ax.set_title("C  missed RSS, or conversion drift?",
                     fontsize=10, fontweight="bold", loc="left")
    else:
        ax.axis("off")
        ax.set_title("C  needs --constrained-assignments", fontsize=10,
                     fontweight="bold", loc="left")

    for a in axes:
        for s in ("top", "right"):
            a.spines[s].set_visible(False)

    fig.suptitle("Do the RSS calls agree with expression?  "
                 "(expression scored with every V gene eligible, so the test is "
                 "independent of the RSS annotation)",
                 fontsize=12, fontweight="bold", y=1.02)
    fig.tight_layout()
    save_figure(fig, args.out_figure)

    for locus in loci:
        sub = [r for r in rows if r["locus"] == locus]
        tot_tx = sum(r["n_tx"] for r in sub) or 1
        print(f"{locus}:", file=sys.stderr)
        for c in CELLS:
            n = sum(1 for r in sub if r["cell"] == c)
            tx = sum(r["n_tx"] for r in sub if r["cell"] == c)
            print(f"   {c:<12} {n:>3} genes  {tx:>5} transcripts "
                  f"({tx/tot_tx*100:.1f}%)", file=sys.stderr)


if __name__ == "__main__":
    main()
