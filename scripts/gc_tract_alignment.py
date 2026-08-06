"""
Show the actual sequence evidence behind each called conversion tract.

Every other figure in this analysis is a summary -- counts, arcs, enrichments.
With only tens of confident events per locus, the individual alignments can
simply be shown, which is far more convincing than any statistic: the reader
sees a contiguous block where the transcript abandons its parent and follows a
donor, then returns to the parent on both sides.

Layout, per event
-----------------
    parent      the germline gene the transcript rearranged from
    donor       the gene supplying the tract
    transcript  coloured by what each base agrees with:

        navy   agrees with the parent (and the donor differs there)
        teal   agrees with the DONOR and not the parent  -- the evidence
        rose   agrees with neither -- an independent point mutation
        grey   parent and donor are identical here, so the base is uninformative

Only informative columns carry any weight, so the run of teal inside the tract
against navy outside it is the whole argument, visible directly.
"""
import argparse
import sys
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
from Bio import Align

from gc_lib import read_fasta, parse_paf, projected_query
from gc_palette import save_figure, SEGMENT, YES, NO, GREY, GREY_DARK, INK

PARENT_C = INK          # follows the parent
DONOR_C = YES           # follows the donor -- the conversion evidence
OTHER_C = NO            # follows neither -- point mutation
UNINF_C = GREY          # parent == donor, tells us nothing


def align_donor(parent_seq, donor_seq, aligner):
    aln = aligner.align(parent_seq, donor_seq)[0]
    out = [None] * len(parent_seq)
    for (ps, pe), (ds, de) in zip(aln.aligned[0], aln.aligned[1]):
        for off in range(pe - ps):
            out[ps + off] = donor_seq[ds + off]
    return out


def read_tsv(path):
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(hdr, line.rstrip("\n").split("\t")))
                for line in fh if line.strip()]


def short(name):
    try:
        return name.split(".")[1]
    except IndexError:
        return name


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--paf", required=True)
    ap.add_argument("--vgene-fasta", required=True)
    ap.add_argument("--assignments", required=True)
    ap.add_argument("--tracts", required=True)
    ap.add_argument("--locus", required=True)
    ap.add_argument("--flank", type=int, default=14,
                    help="bases of parent-following sequence to show either side, "
                         "so the tract has visible boundaries (default: 14)")
    ap.add_argument("--max-events", type=int, default=18)
    ap.add_argument("--primary-only", action="store_true",
                    help="skip tracts where a better-supported donor exists")
    ap.add_argument("--out-figure", required=True)
    args = ap.parse_args()

    vgenes = read_fasta(args.vgene_fasta)
    parent_of = {}
    for r in read_tsv(args.assignments):
        if r["locus"] == args.locus and r["best_gene"] in vgenes:
            parent_of[r["transcript"]] = r["best_gene"]

    tracts = [t for t in read_tsv(args.tracts) if t.get("significant") == "True"]
    if args.primary_only:
        tracts = [t for t in tracts if t.get("primary_donor", "True") == "True"]
    if not tracts:
        raise SystemExit("no significant tracts to draw")

    # Collapse identical events.  The same conversion, inherited through clonal
    # expansion, appears once per transcript in the tract table; drawing each
    # copy would fill the figure with sixteen rows of the same alignment.  What
    # is worth showing is the set of DISTINCT events, each labelled with how many
    # transcripts carry it -- that count is itself informative, since a tract
    # seen in many transcripts is either an early event or a recurrent one.
    groups = defaultdict(list)
    for t in tracts:
        groups[(t["parent"], t["donor"], t["start"], t["end"])].append(t)
    reps = []
    for members in groups.values():
        rep = dict(max(members, key=lambda t: int(t["n_support"])))
        rep["n_transcripts"] = len(members)
        reps.append(rep)
    n_distinct = len(reps)
    n_events = len(tracts)
    # strongest evidence first, then most widely shared
    reps.sort(key=lambda t: (-int(t["n_support"]), -t["n_transcripts"]))
    tracts = reps[:args.max_events]

    # projection of each needed transcript onto its parent
    want = {t["transcript"] for t in tracts}
    proj_of = {}
    for rec in parse_paf(args.paf):
        if rec.query not in want:
            continue
        p = parent_of.get(rec.query)
        if p is None or rec.target != p:
            continue
        pr = projected_query(rec, vgenes[p])
        if pr is None:
            continue
        prev = proj_of.get(rec.query)
        if prev is None or sum(x is not None for x in pr) > sum(x is not None for x in prev):
            proj_of[rec.query] = pr

    aligner = Align.PairwiseAligner(mode="global", open_gap_score=-10,
                                    extend_gap_score=-0.5, match_score=2,
                                    mismatch_score=-1)
    dcache = {}

    events = []
    for t in tracts:
        tid, par, don = t["transcript"], t["parent"], t["donor"]
        proj = proj_of.get(tid)
        if proj is None or par not in vgenes or don not in vgenes:
            continue
        key = (par, don)
        if key not in dcache:
            dcache[key] = align_donor(vgenes[par], vgenes[don], aligner)
        dproj = dcache[key]
        s, e = int(t["start"]), int(t["end"])
        lo, hi = max(0, s - args.flank), min(len(vgenes[par]) - 1, e + args.flank)
        cols = []
        for i in range(lo, hi + 1):
            pb = vgenes[par][i]
            db = dproj[i]
            tb = proj[i]
            if tb is None or tb == "-":
                cls = "gap"
            elif db is None or db == pb:
                cls = "uninf"
            elif tb == db:
                cls = "donor"
            elif tb == pb:
                cls = "parent"
            else:
                cls = "other"
            cols.append((i, pb, db if db else ".", tb if tb else "-", cls))
        events.append({"tid": tid, "parent": par, "donor": don,
                       "start": s, "end": e, "m": int(t["n_support"]),
                       "allowed": t.get("donor_allowed") == "True",
                       "n_alt": int(t.get("n_candidate_donors", 1) or 1),
                       "n_tx": t.get("n_transcripts", 1),
                       "cols": cols, "lo": lo, "hi": hi,
                       "glen": len(vgenes[par])})

    if not events:
        raise SystemExit("no events could be rendered")

    width = max(len(ev["cols"]) for ev in events)
    fig_h = 1.05 * len(events) + 1.9
    fig, ax = plt.subplots(figsize=(max(11.0, width * 0.135 + 4.2), fig_h))

    colour = {"parent": PARENT_C, "donor": DONOR_C, "other": OTHER_C,
              "uninf": UNINF_C, "gap": GREY}
    row_h = 1.0
    fs = 6.2 if width > 70 else 7.4

    for k, ev in enumerate(events):
        y = -k * row_h
        # tract span marker across the top of the block
        n = len(ev["cols"])
        s_off = ev["start"] - ev["lo"]
        e_off = ev["end"] - ev["lo"]
        ax.add_patch(Rectangle((s_off - .5, y - 0.40), e_off - s_off + 1, 0.80,
                               facecolor=DONOR_C, alpha=0.13, edgecolor=DONOR_C,
                               lw=0.9, zorder=1))
        for j, (i, pb, db, tb, cls) in enumerate(ev["cols"]):
            ax.text(j, y + 0.22, pb, ha="center", va="center", fontsize=fs,
                    family="monospace", color=GREY_DARK, zorder=3)
            ax.text(j, y + 0.02, db, ha="center", va="center", fontsize=fs,
                    family="monospace", color=GREY_DARK, zorder=3)
            ax.text(j, y - 0.20, tb, ha="center", va="center", fontsize=fs,
                    family="monospace", color=colour[cls], zorder=3,
                    fontweight="bold" if cls in ("donor", "other") else "normal")
        flag = "" if ev["allowed"] else "   ⚠ donor was DELETED"
        alt = f"   ·  {ev['n_alt']-1} competing donor(s)" if ev["n_alt"] > 1 else ""
        ax.text(-1.6, y + 0.22, f"parent {short(ev['parent'])}", ha="right",
                va="center", fontsize=7, color=GREY_DARK)
        ax.text(-1.6, y + 0.02, f"donor  {short(ev['donor'])}", ha="right",
                va="center", fontsize=7, color=GREY_DARK)
        seen = (f"seen in {ev['n_tx']} transcripts" if ev["n_tx"] > 1
                else "seen in 1 transcript")
        ax.text(-1.6, y - 0.20, "transcript", ha="right",
                va="center", fontsize=7, color=INK, fontweight="bold")
        # The count is the headline for each row, so it gets its own weight
        # rather than being buried in the trailing metadata.
        ax.text(n + 1.0, y + 0.14, f"×{ev['n_tx']}", ha="left", va="center",
                fontsize=9.5, fontweight="bold",
                color=DONOR_C if ev["allowed"] else NO)
        ax.text(n + 1.0, y - 0.10,
                f"{seen}\nm={ev['m']}  V:{ev['start']}–{ev['end']}{alt}{flag}",
                ha="left", va="center", fontsize=6.4,
                color=INK if ev["allowed"] else NO)

    ax.set_xlim(-13, width + 15)
    ax.set_ylim(-row_h * len(events) + 0.4, 0.85)
    ax.axis("off")

    handles = [
        Patch(facecolor=DONOR_C, label="transcript follows the DONOR (evidence)"),
        Patch(facecolor=PARENT_C, label="transcript follows the parent"),
        Patch(facecolor=OTHER_C, label="matches neither (point mutation)"),
        Patch(facecolor=UNINF_C, label="parent = donor here (uninformative)"),
    ]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.01),
              ncol=4, fontsize=8, framealpha=0.95)

    fig.suptitle(
        f"{args.locus} — sequence evidence for the {len(events)} best-supported "
        f"DISTINCT conversion tracts "
        f"({n_distinct} distinct events across {n_events} transcript-level calls)\n"
        f"identical events collapsed; ×N = how many transcripts carry that exact "
        f"tract · shaded = called tract · m = donor-diagnostic positions",
        fontsize=11, fontweight="bold", y=0.995)
    fig.tight_layout()
    save_figure(fig, args.out_figure)
    print(f"{args.locus}: drew {len(events)} events to {args.out_figure}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
