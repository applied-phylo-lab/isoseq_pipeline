"""
Conversion tracts shown across the WHOLE V gene, not cropped to the tract.

Why the cropped version is not enough
-------------------------------------
The main figure's sequence panel shows the tract plus 12 bp either side.  That is
enough to see that a block of donor-matching bases exists, but not enough to rule
out the alternative that worries anyone reading it: that the "parent" was simply
misassigned and the transcript is really a descendant of the donor throughout.
The crop cannot distinguish those, because it never shows what the transcript
does away from the tract.

This figure shows every informative position across the full assessable region.
The pattern that matters is parent - parent - parent - DONOR BLOCK - parent -
parent: the transcript follows its parent on BOTH flanks and switches only over
the tract.  A misassignment produces a completely different picture -- donor
agreement spread across the whole gene, with no flanks.

Only informative positions are drawn
------------------------------------
A position where parent and donor are identical says nothing about which of them
the transcript came from, and roughly 90% of positions are like that.  Drawing
them would push the informative ones apart until the switch is invisible.  So the
x axis is the informative positions in order, with their real V coordinates on
the axis and the tract shaded, and a schematic strip above showing where those
positions actually sit along the gene.

Events are chosen for FLANK EVIDENCE, not for support
-----------------------------------------------------
Ranking by the number of donor-diagnostic positions would pick the events with
the longest tracts, which are the least informative for this question -- a tract
covering most of the gene leaves no flank to check.  Instead events are ranked by
how much parent-following evidence exists on BOTH sides of the tract, which is
exactly the evidence the figure exists to show.
"""
import argparse
import csv
import sys
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from Bio import Align
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle

from gc_lib import read_fasta, parse_paf, projected_query
from gc_palette import save_figure, LOCUS, SEGMENT, GREY, GREY_DARK, INK, NO, YES

DONOR_C = YES
PARENT_C = INK
OTHER_C = NO
# parent and donor identical here -- the transcript matches both, so the position
# carries no evidence either way. Deliberately the palest thing in the figure.
UNINF_C = "#E4E7EB"


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def short(name):
    p = name.split(".")
    return p[1] if len(p) > 1 else name


def build(paf, vgenes, tracts, locus, min_flank, rss=None):
    """One entry per distinct (parent, donor, start, end), with full-gene columns."""
    sig = [t for t in tracts if t.get("significant") == "True"
           and t.get("primary_donor", "True") == "True"]
    groups = defaultdict(list)
    for t in sig:
        groups[(t["parent"], t["donor"], t["start"], t["end"])].append(t)
    reps = []
    for members in groups.values():
        r = dict(max(members, key=lambda t: int(t["n_support"])))
        r["n_transcripts"] = len(members)
        reps.append(r)

    want = {t["transcript"] for t in reps}
    parent_of = {t["transcript"]: t["parent"] for t in reps}
    proj_of = {}
    for rec in parse_paf(paf):
        if rec.query not in want:
            continue
        if rec.target != parent_of.get(rec.query):
            continue
        pr = projected_query(rec, vgenes[rec.target])
        if pr is None:
            continue
        prev = proj_of.get(rec.query)
        if prev is None or sum(x is not None for x in pr) > sum(x is not None for x in prev):
            proj_of[rec.query] = pr

    aligner = Align.PairwiseAligner(mode="global", open_gap_score=-10,
                                    extend_gap_score=-0.5, match_score=2,
                                    mismatch_score=-1)
    cache, events = {}, []
    for t in reps:
        tid, par, don = t["transcript"], t["parent"], t["donor"]
        proj = proj_of.get(tid)
        if proj is None:
            continue
        if (par, don) not in cache:
            aln = aligner.align(vgenes[par], vgenes[don])[0]
            dp = [None] * len(vgenes[par])
            for (ps, pe), (ds, de) in zip(aln.aligned[0], aln.aligned[1]):
                for off in range(pe - ps):
                    dp[ps + off] = vgenes[don][ds + off]
            cache[(par, don)] = dp
        dproj = cache[(par, don)]
        s, e = int(t["start"]), int(t["end"])

        # allcols keeps EVERY position, cols keeps only the informative ones.
        # Both are needed. The rows below are drawn from cols, because a position
        # where parent and donor agree cannot say which of them the transcript
        # came from and including them would spread the diagnostic positions so
        # far apart that the switch stops being visible. But dropping them
        # silently hides how much of the gene is uninformative in the first
        # place, which is the reader's main handle on whether the parent and
        # donor are even distinguishable -- so allcols is drawn at real scale in
        # the strip above, with its own colour for "parent = donor".
        allcols, cols = [], []
        for i in range(len(vgenes[par])):
            pb, db, tb = vgenes[par][i], dproj[i], proj[i]
            if tb is None or tb == "-":
                cls = "gap"
            elif db is None or db == pb:
                cls = "uninf"              # parent and donor agree: no evidence
            elif tb == db:
                cls = "donor"
            elif tb == pb:
                cls = "parent"
            else:
                cls = "other"
            allcols.append((i, cls))
            if cls not in ("uninf", "gap") or (db is not None and db != pb):
                if cls != "uninf":
                    cols.append((i, pb, db, tb, cls))
        if not cols:
            continue
        left = sum(1 for i, *_ , c in cols if i < s and c == "parent")
        right = sum(1 for i, *_ , c in cols if i > e and c == "parent")
        inside = sum(1 for i, *_ , c in cols if s <= i <= e and c == "donor")
        # Donor agreement OUTSIDE the called tract. This is the number that
        # decides whether the switch is clean: the detector calls a tract only
        # where donor-diagnostic positions are contiguous within max_gap and
        # numerous enough to pass the support cutoff, so a genuine conversion
        # with ragged edges leaves donor-following positions just outside it.
        # A LOT of them, spread the length of the gene, means something else --
        # either a second tract the detector split off, or a parent/donor pair
        # too close to tell apart. Reporting it is the only way the reader can
        # judge which, so it goes in the panel title rather than being left out.
        d_out = sum(1 for i, *_ , c in cols if not (s <= i <= e) and c == "donor")
        other = sum(1 for i, *_ , c in cols if c == "other")
        n_uninf = sum(1 for _i, c in allcols if c == "uninf")
        # 5' ANCHOR -- the test for which gene is the parent.
        #
        # NOT the 3' end. In chicken, conversion tracts can run past the V 3' end
        # and into the D region, so a transcript following its donor at the
        # junction is an expected outcome of conversion, not evidence that the
        # parent was misassigned. Anchoring there would reject real events.
        #
        # The 5' end (FR1) is the conserved end and is where the transcript
        # reliably still follows the gene it rearranged from, so that is where
        # parent identity is read. Overall identity is useless for this: a
        # heavily converted transcript ends up closer to its donor by
        # construction, which is what conversion does.
        head_cols = [c for _i, *_r, c in cols][:8]
        p5_par = head_cols.count("parent")
        p5_don = head_cols.count("donor")
        events.append({
            "parent": par, "donor": don, "start": s, "end": e, "locus": locus,
            "m": int(t["n_support"]), "n_tx": t["n_transcripts"],
            "cols": cols, "allcols": allcols, "left": left, "right": right,
            "inside": inside, "d_out": d_out, "other": other,
            "n_uninf": n_uninf, "p5_par": p5_par, "p5_don": p5_don,
            "p_rss": (rss or {}).get(par, "") == "rss_present",
            "d_rss": (rss or {}).get(don, "") == "rss_present",
            "glen": len(vgenes[par]),
        })
    # both flanks must carry evidence, and the weaker flank is what is ranked on
    events = [e for e in events if min(e["left"], e["right"]) >= min_flank]
    events.sort(key=lambda e: (-min(e["left"], e["right"]), -e["inside"]))
    return events


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--paf", required=True)
    ap.add_argument("--vgene-fasta", required=True)
    ap.add_argument("--tracts", required=True)
    ap.add_argument("--locus", required=True)
    ap.add_argument("--extra", action="append", default=[],
                    help="LOCUS=paf,fasta,tracts for a second locus")
    ap.add_argument("--rss-annotation", required=True,
                    help="marks whether parent and donor each carry an RSS")
    ap.add_argument("--min-flank", type=int, default=0,
                    help="parent-following positions required on EACH side; "
                         "0 draws every tract, which is the point of the figure")
    ap.add_argument("--out-figure", required=True)
    ap.add_argument("--out-table", required=True)
    args = ap.parse_args()

    sets = [(args.locus, args.paf, args.vgene_fasta, args.tracts)]
    for spec in args.extra:
        loc, rest = spec.split("=", 1)
        p, f, t = rest.split(",")
        sets.append((loc, p, f, t))

    rssmap = {r["gene"]: r["rss_state"] for r in read_tsv(args.rss_annotation)}
    allev = []
    for loc, paf, fa, tr in sets:
        ev = build(paf, read_fasta(fa), read_tsv(tr), loc, args.min_flank, rssmap)
        print(f"{loc}: {len(ev)} distinct tracts with >={args.min_flank} "
              f"parent-following positions on both flanks", file=sys.stderr)
        allev.append((loc, ev))

    with open(args.out_table, "w") as fh:
        fh.write("locus\tparent\tdonor\tstart\tend\tn_support\tn_transcripts\t"
                 "informative_positions\tparent_left\tdonor_inside\tparent_right\t"
                 "donor_outside_tract\tmatches_neither\t"
                 "first8_follow_parent\tfirst8_follow_donor\t"
                 "parent_has_rss\tdonor_has_rss\n")
        for loc, ev in allev:
            for e in ev:
                fh.write(f"{loc}\t{e['parent']}\t{e['donor']}\t{e['start']}\t"
                         f"{e['end']}\t{e['m']}\t{e['n_tx']}\t{len(e['cols'])}\t"
                         f"{e['left']}\t{e['inside']}\t{e['right']}\t"
                         f"{e['d_out']}\t{e['other']}\t"
                         f"{e['p5_par']}\t{e['p5_don']}\t"
                         f"{e['p_rss']}\t{e['d_rss']}\n")

    # ── proof score ───────────────────────────────────────────────────────────
    # "Proof" here means the parent - DONOR - parent signature, which needs two
    # things at once and is not captured by either alone:
    #
    #   anchor5 = fraction of the FIRST 8 informative positions following the
    #             parent. FR1 is the conserved end and is where a transcript
    #             still follows what it rearranged from, so this is the test for
    #             whether the parent is the parent. Not the 3' end: a chicken
    #             conversion tract can run into D, so donor agreement at the
    #             junction is expected rather than disqualifying.
    #   purity  = donor-following positions inside the tract / all donor-following
    #             positions. 1.0 means donor agreement stops at the tract; low
    #             means it leaks across the gene.
    #
    # Ranking on the product puts events that anchor to their parent at 5' AND
    # keep their donor agreement inside the tract at the top, and the ambiguous
    # ones at the bottom where they can be seen rather than quietly dropped.
    for _loc, ev in allev:
        for e in ev:
            e["flank"] = min(e["left"], e["right"])
            e["anchor5"] = e["p5_par"] / max(1, e["p5_par"] + e["p5_don"])
            e["purity"] = e["inside"] / max(1, e["inside"] + e["d_out"])
            e["proof"] = e["anchor5"] * e["purity"]
    show = [e for _loc, ev in allev for e in ev]
    show.sort(key=lambda e: (e["locus"], -e["proof"], -e["flank"]))
    if not show:
        raise SystemExit("no tracts to draw")

    colour = {"parent": PARENT_C, "donor": DONOR_C, "other": OTHER_C,
              "gap": GREY, "uninf": UNINF_C}
    n = len(show)
    maxlen = max(e["glen"] for e in show)
    fig, ax = plt.subplots(figsize=(12.2, 0.30 * n + 2.0))

    for k, e in enumerate(show):
        y = n - 1 - k                       # best proof at the top
        for i, cls in e["allcols"]:
            ax.add_patch(Rectangle((i, y - 0.30), 1.02, 0.60,
                                   facecolor=colour[cls], edgecolor="none",
                                   zorder=2))
        # the called tract, as a bar directly above its own row
        ax.add_patch(Rectangle((e["start"], y + 0.32),
                               e["end"] - e["start"] + 1, 0.13,
                               facecolor=SEGMENT["J"], edgecolor="none",
                               zorder=3))
        ax.text(maxlen + 6, y,
                f"{e['p5_par']:>2}/{e['p5_par'] + e['p5_don']:<2} "
                f"{e['inside']:>2} | {e['d_out']:>2}   "
                f"×{e['n_tx']}  m={e['m']}",
                va="center", ha="left", fontsize=6.8, color=GREY_DARK,
                family="monospace")

    # RSS state of BOTH genes, as a pair of dots left of each row. Which gene can
    # be rearranged at all is the other half of the parent call -- the sequence
    # evidence says which gene the transcript follows at 5', and this says which
    # gene was even eligible to be the parent. Where they disagree, the call is
    # worth a second look, so both have to be visible at once.
    for k, e in enumerate(show):
        y = n - 1 - k
        for dx, has in ((-30, e["p_rss"]), (-17, e["d_rss"])):
            ax.plot([dx], [y], marker="o", ms=4.0,
                    color=SEGMENT["V"] if has else "white",
                    markeredgecolor="black", markeredgewidth=0.55,
                    clip_on=False, zorder=6)
    ax.text(-30, n - 0.35, "P", fontsize=6.5, ha="center", va="bottom",
            color=INK, fontweight="bold")
    ax.text(-17, n - 0.35, "D", fontsize=6.5, ha="center", va="bottom",
            color=INK, fontweight="bold")

    ax.set_yticks(range(n))
    ax.set_yticklabels([f"{e['locus']}  {short(e['parent'])} ← {short(e['donor'])}"
                        for e in reversed(show)], fontsize=6.8)
    for lab, e in zip(ax.get_yticklabels(), reversed(show)):
        lab.set_color(LOCUS.get(e["locus"], INK))
    ax.set_ylim(-0.8, n - 0.2)
    ax.set_xlim(-38, maxlen + 96)
    ax.set_xlabel("position in the parent V gene (bp, true scale)", fontsize=9)
    ax.text(maxlen + 6, n - 0.35,
            "5'par  in | out   ×tx", fontsize=6.8, color=INK,
            family="monospace", va="bottom", ha="left")
    for sp in ("left", "right", "top"):
        ax.spines[sp].set_visible(False)
    ax.tick_params(axis="y", length=0)

    handles = [
        Patch(facecolor=PARENT_C, label="follows the PARENT"),
        Patch(facecolor=DONOR_C, label="follows the DONOR"),
        Patch(facecolor=OTHER_C, label="matches neither"),
        Patch(facecolor=UNINF_C, label="parent = donor (no evidence)"),
        Patch(facecolor=GREY, label="not covered"),
        Patch(facecolor=SEGMENT["J"], label="called tract"),
        Line2D([], [], marker="o", ls="none", ms=6, color=SEGMENT["V"],
               markeredgecolor="black", markeredgewidth=0.55,
               label="P/D dots: filled = that gene has an RSS"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=8.2,
               frameon=True, framealpha=0.95, bbox_to_anchor=(0.5, 0.004))
    fig.tight_layout(rect=(0, 0.028 + 0.32 / max(1, n), 1, 1))
    save_figure(fig, args.out_figure)

    print(f"{n} tracts drawn; 5' anchor >=0.75 and purity >=0.6: "
          f"{sum(1 for e in show if e['anchor5'] >= 0.75 and e['purity'] >= 0.6)}",
          file=sys.stderr)
    weak = [e for e in show if e["anchor5"] < 0.5]
    print(f"  tracts whose 5' end follows the DONOR more than the parent "
          f"(parent call suspect): {len(weak)}", file=sys.stderr)
    for e in weak:
        print(f"    {e['locus']} {short(e['parent'])} <- {short(e['donor'])}  "
              f"5' parent {e['p5_par']} donor {e['p5_don']}  "
              f"parentRSS={e['p_rss']} donorRSS={e['d_rss']}", file=sys.stderr)


if __name__ == "__main__":
    main()
