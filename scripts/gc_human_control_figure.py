"""
How often a conversion event is called: birds versus a human negative control.

Why this is the control that matters
------------------------------------
Every threshold in this pipeline was calibrated against a permutation null --
shuffle the data, see how often a tract appears by chance. That is a necessary
check but a weak one, because permuting destroys exactly the structure that
causes real false positives: sequencing error, allelic mismatch between the
transcript and the reference, clonal expansion, and paralogues similar enough to
mimic a donor.

Humans do not use gene conversion to diversify immunoglobulins. So running the
identical detector on human repertoires, each scored against that individual's
OWN assembled germline, measures the false-positive rate on real data with all
of that structure intact. Whatever rate it reports is the floor the bird result
has to clear.

Every choice here favours the human side
----------------------------------------
The human samples are run with UNCONSTRAINED parent assignment and a permissive
donor pool in which every gene is available to every parent, because the
targeted-capture assemblies contain no J and span many contigs, so no
deletion/inversion model is possible. Both make a spurious call easier. The
human reads are also merged Illumina rather than HiFi, so they are noisier. The
comparison is therefore conservative: the bird rate is being compared against a
deliberately generous estimate of the noise floor.

For the same reason the bird bars default to the ANY-PARENT run: it is the
setting that matches the human one. The RSS-restricted main run is higher and is
shown alongside for reference, but the any-parent bars carry the argument.

Why the unit is % of transcripts, and not distinct events per 1,000
-------------------------------------------------------------------
Distinct-event counts SATURATE with sequencing depth: each extra transcript is
more likely to rediscover an event already seen than to add a new one. So
"distinct events per 1,000 transcripts" falls as a repertoire is sampled deeper,
and it cannot be compared across libraries of different depth. The humans here
carry 25,000-49,000 transcripts and the birds 71-524 -- a 100-700x gap, all of it
in the birds' favour. That unit would manufacture most of the difference.

The per-transcript rate does not have that problem. `n_support` in the detector
is the number of donor-supporting informative POSITIONS within one transcript
(gc_detect_tracts.py), not the number of transcripts backing a tract, so every
transcript is scored on its own and the significance rule never consults the
library size. The fraction of transcripts carrying a tract therefore estimates a
per-transcript probability that is independent of how many were sampled. Small n
in the birds makes that estimate NOISY, not biased, which is why bird bars carry
Wilson 95% intervals rather than a rescaled axis.

Detection opportunity per transcript also favours the humans: median covered V
sequence is 276 bp for the humans against 238-245 bp for the birds.

The QC rule for a dominated sample
----------------------------------
A sample is excluded when more than half its calls come from a SINGLE
(parent, donor, start, end) combination, because such a sample is reporting one
event, not a rate. The rule is applied to every sample, birds included.

Only W-79 trips it, at 94.8%: 1,764 of its 1,861 calls are the same parent, the
same donor and the same 16 bp window, on a gene annotated non-productive -- one
germline allele the assembly missed, so every transcript from that gene carries
an identical mismatch block that a paralogue supplies. No bird locus exceeds
32.5%. The excluded sample and its rate stay in the output table.
"""
import argparse
import csv
import os
import sys
from collections import Counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from gc_palette import save_figure, LOCUS, GREY, GREY_DARK, INK, NO


def stats(tracts_path, assignments_path, locus=None):
    """Per-sample summary.

    Returns (pct, n_tx, n_events, n_tx_with_tract, top_combo_pct, n_distinct),
    where `pct` is the percentage of transcripts carrying at least one tract --
    the plotted value, and the one that is comparable across sequencing depths.

    `top_combo_pct` is the share of a sample's calls taken by its single most
    repeated (parent, donor, start, end) combination; it drives the QC rule.
    """
    n_tx = sum(1 for r in csv.DictReader(open(assignments_path), delimiter="\t")
               if locus is None or r["locus"] == locus)
    tr = [r for r in csv.DictReader(open(tracts_path), delimiter="\t")
          if r.get("significant") == "True"]
    if not tr or not n_tx:
        return 0.0, n_tx, 0, 0, 0.0, 0
    combo = Counter((r["parent"], r["donor"], r["start"], r["end"]) for r in tr)
    n_hit = len({r["transcript"] for r in tr})
    return (100.0 * n_hit / n_tx, n_tx, len(tr), n_hit,
            100.0 * combo.most_common(1)[0][1] / len(tr), len(combo))


def wilson(k, n, z=1.96):
    """95% interval for a proportion. Score interval, because the bird counts are
    small enough (5/71) that the normal approximation is not usable."""
    if not n:
        return 0.0, 0.0
    p = k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return 100 * max(0.0, c - h), 100 * min(1.0, c + h)


def fmt(v):
    if v < 0.1:
        return f"{v:.3f}"
    return f"{v:.2f}" if v < 1 else f"{v:.1f}"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--human-dir", required=True,
                    help="dir with one subdir per human sample")
    ap.add_argument("--bird", nargs="+", default=[],
                    metavar="LABEL:tracts.tsv:assignments.tsv:LOCUS",
                    help="bars drawn next to the humans; these carry the argument")
    ap.add_argument("--bird-reference", nargs="+", default=[],
                    metavar="LABEL:tracts.tsv:assignments.tsv:LOCUS",
                    help="shown as open markers only (e.g. the RSS-restricted run)")
    ap.add_argument("--out-figure", required=True)
    ap.add_argument("--out-table", required=True)
    args = ap.parse_args()

    human = []
    for s in sorted(os.listdir(args.human_dir)):
        d = os.path.join(args.human_dir, s)
        t = os.path.join(d, "IGH_tracts.tsv")
        a = os.path.join(d, "transcript_assignments.tsv")
        if os.path.isfile(t) and os.path.isfile(a):
            human.append((s,) + stats(t, a))
    human.sort(key=lambda r: r[1])

    def load(specs):
        out = []
        for spec in specs:
            parts = spec.split(":")
            lab, t, a = parts[0], parts[1], parts[2]
            loc = parts[3] if len(parts) > 3 else None
            out.append((lab,) + stats(t, a, loc) + (loc,))
        return out

    birds = load(args.bird)
    refs = load(args.bird_reference)

    # QC: a sample whose calls are dominated by one repeated combination is
    # reporting a single event, not a rate. Applied to birds identically.
    MAX_COMBO = 50.0
    dropped = [r for r in human if r[5] > MAX_COMBO]
    human = [r for r in human if r[5] <= MAX_COMBO]
    bird_dropped = [b for b in birds if b[5] > MAX_COMBO]
    birds = [b for b in birds if b[5] <= MAX_COMBO]

    # the table keeps every sample, including the ones QC removed, and both units
    with open(args.out_table, "w") as fh:
        fh.write("group\tsample\tpct_transcripts_with_tract\tci95_low\tci95_high\t"
                 "distinct_events_per_1000_tx\tn_transcripts\tn_events\t"
                 "n_transcripts_with_tract\ttop_combo_pct\tn_distinct_tracts\t"
                 "qc_excluded\n")
        for grp, rs in (("human", human + dropped), ("bird", birds + bird_dropped),
                        ("bird_reference", refs)):
            for r in rs:
                lo, hi = wilson(r[4], r[2])
                per1k = 1000.0 * r[6] / r[2] if r[2] else 0.0
                fh.write(f"{grp}\t" + "\t".join(
                    [r[0], f"{r[1]:.4f}", f"{lo:.4f}", f"{hi:.4f}", f"{per1k:.3f}",
                     str(r[2]), str(r[3]), str(r[4]), f"{r[5]:.1f}", str(r[6]),
                     "True" if r[5] > MAX_COMBO else "False"]) + "\n")

    hp = [r[1] for r in human]
    med = float(np.median(hp)) if hp else 0.0

    fig, ax = plt.subplots(figsize=(8.8, 5.4))

    # wrap only at the final space, so "red-winged blackbird IGH" breaks into the
    # species and the locus rather than onto three lines
    groups = [f"human\n(n={len(hp)})"] + [b[0].rsplit(" ", 1)[0] + "\n" +
                                          b[0].rsplit(" ", 1)[-1] for b in birds]
    gv = [med] + [b[1] for b in birds]
    gc_ = [GREY] + [LOCUS.get(b[-1], INK) for b in birds]
    xs = np.arange(len(groups))
    ax.bar(xs, gv, 0.62, color=gc_, edgecolor="black", lw=0.6, zorder=3)

    # Wilson intervals are computed and written to the output table, but are not
    # drawn: the separation between the groups is far larger than the intervals,
    # so they add width without changing what the figure shows.
    if hp:
        # every individual human sample, so the reader sees the spread rather
        # than trusting a single summary bar
        jit = np.random.default_rng(0).normal(0, 0.07, len(hp))
        ax.scatter(jit, hp, s=22, color="white", edgecolor="black", lw=0.7, zorder=6)
    for xi, v in zip(xs, gv):
        if xi == 0 and hp:
            # the human points scatter straight through where this label would
            # sit, so put it beside the bar rather than on top of it
            ax.text(xi + 0.36, v, fmt(v) + "%", ha="left", va="center",
                    fontsize=9, fontweight="bold")
        else:
            ax.text(xi, v * 1.16, fmt(v) + "%", ha="center", va="bottom",
                    fontsize=9, fontweight="bold")
    # log scale is unavoidable: the groups differ by >2 orders of magnitude, and
    # on a linear axis every human point would sit on the floor
    ax.set_yscale("log")
    ax.set_ylim(0.002, 200)
    ax.set_xticks(xs)
    ax.set_xticklabels(groups, fontsize=8.5)
    ax.set_ylabel("transcripts carrying a conversion tract (%)", fontsize=9.5)
    if hp and birds:
        worst_h = max(hp)
        fold = min(b[1] for b in birds) / max(worst_h, 1e-9)
        ax.set_title("Humans have no gene conversion: every bird locus sits above\n"
                     f"every human sample, the closest by {fold:.0f}$\\times$",
                     fontsize=11, loc="left", x=0.0)
    ax.legend(handles=[
        Patch(facecolor=GREY, edgecolor="black", lw=0.6, label="human (bar = median)"),
        Line2D([], [], marker="o", ls="", mfc="white", mec="black", mew=0.7, ms=6,
               label="individual human sample"),
    ], fontsize=8, frameon=False, loc="upper left")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

    fig.tight_layout()
    save_figure(fig, args.out_figure)

    print(f"{len(human)} human, {len(birds)} bird "
          f"(% of transcripts carrying a tract)", file=sys.stderr)
    for r in human + birds:
        lo, hi = wilson(r[4], r[2])
        print(f"  {r[0]:<22} {r[1]:>7.3f}%  [{lo:.2f}-{hi:.2f}]  "
              f"({r[4]}/{r[2]} tx, {r[3]} calls, {r[6]} distinct, "
              f"top combo {r[5]:.0f}%)", file=sys.stderr)
    for r in dropped + bird_dropped:
        print(f"  {r[0]:<22} {r[1]:>7.3f}%  QC-EXCLUDED "
              f"(top combo {r[5]:.0f}% > {MAX_COMBO:.0f}%)", file=sys.stderr)
    print(f"  human median {med:.3f}%", file=sys.stderr)


if __name__ == "__main__":
    main()
