"""
Where along the V gene conversion actually lands, against how variable that
position is between genes.

The question this answers
------------------------
Tracts are not spread evenly along the V gene.  In IGL they fall in four windows
covering 57 of 283 positions; in IGH each parent contributes one window.  The
obvious question is whether those windows are anything more than "where the
detector can see", because a tract can only be called where the donor differs
from the parent -- so the divergence profile is plotted underneath the tract
windows rather than being left to the reader's imagination.

Landmarks
---------
The framework motifs are marked because they are what turn a coordinate into a
statement.  FR2 ends with ...L-E-W-(V/I)-R-Q..., and CDR2 begins immediately
after it; CDR1 sits immediately before the W-(V/I)-R-Q.  Nucleotide-level motif
search is used rather than translation, because most of these genes are
pseudogenes whose reading frame cannot be assumed.
"""
import argparse
import csv
import re
import sys
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from Bio import Align
from Bio.Seq import Seq
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from gc_lib import read_fasta
from gc_palette import save_figure, LOCUS, SEGMENT, GREY, GREY_DARK, INK, NO, YES

# Framework landmarks, found at PROTEIN level and mapped back to nucleotides.
#
# Nucleotide regexes were tried first and are not good enough: synonymous
# variation between paralogues breaks them, so three of eight parents went
# unannotated and the failure was easy to misread as pseudogene degeneracy.
# These are all real V genes with intact frameworks, and the motifs below are
# conserved at the amino acid level where the nucleotide is not.
FR2_START = re.compile(r"W[VILYFM][RKQ]Q")       # ...W-x-R-Q, first residues of FR2
# L-E-W-V, last residues of FR2 (heavy). Y is allowed at the W position:
# several of these genes carry L-E-Y-V and excluding them left real CDR2
# boundaries unmarked.
FR2_END_H = re.compile(r"[LIVM]E[WY][VILMA]")
FR2_END_L = re.compile(r"[LIVMT][LIVMT][IVL]Y")  # ...V-T-V-I-Y, FR2 -> CDR2 (light)
FR3_END = re.compile(r"Y[YFHC][CV]")             # ...Y-Y-C, the cysteine closing FR3


def frame_and_protein(seq):
    """Reading frame chosen by V-domain content, not by stop count.

    Picking the frame with fewest stops is wrong: a short V exon often has zero
    stops in more than one frame, and the tie then goes to whichever frame is
    listed first. That silently mis-framed two of the eight parents here and made
    their frameworks look degenerate. Scoring on how many V landmarks a frame
    actually contains picks the coding frame directly.
    """
    best = (-1, 0, 0, "")
    for f in range(3):
        aa = str(Seq(seq[f:len(seq) - ((len(seq) - f) % 3)]).translate())
        hits = sum(bool(r.search(aa)) for r in (FR2_START, FR2_END_H, FR2_END_L, FR3_END))
        cand = (hits, -aa.count("*"), -f, aa)
        if cand > best:
            best = cand
    # the FRAME itself is returned: it cannot be recovered afterwards from
    # len(seq) - 3*len(aa), which gives the TRAILING remainder, not the leading
    # offset, and silently shifted every landmark by up to 2 bp.
    return -best[2], best[3]


def landmarks(seq):
    """[(nt position, label)] for the framework boundaries this gene shows."""
    off, aa = frame_and_protein(seq)
    out = []
    m = FR2_START.search(aa)
    if m:
        out.append((off + m.start() * 3, "CDR1|FR2"))
    for rx in (FR2_END_H, FR2_END_L):
        m = rx.search(aa)
        if m:
            # heavy: CDR2 starts after the W. light: CDR2 starts AT the Y, so the
            # boundary is one residue earlier.
            end = m.end() if rx is FR2_END_H else m.end() - 1
            out.append((off + end * 3, "FR2|CDR2"))
            break
    m = FR3_END.search(aa)
    if m:
        out.append((off + m.end() * 3, "FR3|CDR3"))
    return out, aa


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def short(g):
    p = g.split(".")
    return p[1] if len(p) > 1 else g


def merge(spans):
    out = []
    for s, e in sorted(spans):
        if out and s <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def divergence(parent, seqs):
    """Per position: fraction of the locus's other genes that differ from parent."""
    al = Align.PairwiseAligner(mode="global", open_gap_score=-10,
                              extend_gap_score=-0.5, match_score=2,
                              mismatch_score=-1)
    ref = seqs[parent]
    diff = np.zeros(len(ref))
    n = np.zeros(len(ref))
    for g, s in seqs.items():
        if g == parent:
            continue
        aln = al.align(ref, s)[0]
        for (ps, pe), (ds, de) in zip(aln.aligned[0], aln.aligned[1]):
            for o in range(pe - ps):
                i = ps + o
                n[i] += 1
                diff[i] += ref[i] != s[ds + o]
    return np.divide(diff, np.maximum(n, 1))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--locus", action="append", required=True,
                    metavar="NAME=fasta,tracts")
    ap.add_argument("--rss-annotation", required=True,
                    help="so each panel can state whether that parent has an RSS")
    ap.add_argument("--out-figure", required=True)
    ap.add_argument("--out-table", required=True)
    args = ap.parse_args()

    rssmap = {r["gene"]: r["rss_state"] for r in read_tsv(args.rss_annotation)}
    panels = []
    for spec in args.locus:
        name, rest = spec.split("=", 1)
        fa, tr = rest.split(",")
        seqs = read_fasta(fa)
        tracts = [t for t in read_tsv(tr) if t.get("significant") == "True"]
        byp = defaultdict(list)
        calls = defaultdict(int)
        donors = defaultdict(set)
        for t in tracts:
            byp[t["parent"]].append((int(t["start"]), int(t["end"])))
        for p in byp:
            byp[p] = merge(byp[p])
        for t in tracts:
            for s, e in byp[t["parent"]]:
                if s <= int(t["start"]) <= e:
                    calls[(t["parent"], s, e)] += 1
                    donors[(t["parent"], s, e)].add(t["donor"])
        panels.append((name, seqs, byp, calls, donors))

    with open(args.out_table, "w") as fh:
        fh.write("locus\tparent\tgene_len\tstart\tend\tlength\tn_calls\tn_donors\t"
                 "parent_has_rss\tmean_divergence_in_window\tmean_divergence_gene\t"
                 "bp_after_FR2_LEW\tbp_before_FR2_WxRQ\n")
        for name, seqs, byp, calls, donors in panels:
            for p, wins in byp.items():
                dv = divergence(p, seqs)
                lm = landmarks(seqs[p])[0]
                lew = [x for x, lab in lm if lab == "FR2|CDR2"]
                wq = [x for x, lab in lm if lab == "CDR1|FR2"]
                for s, e in wins:
                    up = [h for h in lew if h <= e]
                    dn = [h for h in wq if h >= s]
                    fh.write(f"{name}\t{short(p)}\t{len(seqs[p])}\t{s}\t{e}\t{e-s+1}\t"
                             f"{calls[(p,s,e)]}\t{len(donors[(p,s,e)])}\t"
                             f"{rssmap.get(p)=='rss_present'}\t"
                             f"{dv[s:e+1].mean():.3f}\t{dv.mean():.3f}\t"
                             f"{s-max(up) if up else 'NA'}\t"
                             f"{min(dn)-e if dn else 'NA'}\n")

    nrow = sum(len(p[2]) for p in panels)
    # one shared x range so window positions are comparable between panels
    XMAX = max(len(x) for _n, sq, _b, _c, _d in panels for x in sq.values())
    fig, axes = plt.subplots(nrow, 1, figsize=(11.4, 1.30 * nrow + 1.3),
                             squeeze=False)
    axes = axes[:, 0]
    k = 0
    for name, seqs, byp, calls, donors in panels:
        col = LOCUS.get(name, INK)
        for p in sorted(byp, key=lambda g: -sum(calls[(g, s, e)] for s, e in byp[g])):
            ax = axes[k]; k += 1
            seq = seqs[p]
            dv = divergence(p, seqs)
            lim = len(seq) - 20
            # smooth only for the eye; the windows and stats use raw values
            w = 7
            sm = np.convolve(dv, np.ones(w) / w, mode="same")
            ax.fill_between(range(lim), sm[:lim], color=col, alpha=0.30, lw=0)
            ax.plot(range(lim), sm[:lim], color=col, lw=1.0)
            ax.axhline(dv[:lim].mean(), color=GREY_DARK, lw=0.8, ls=":", zorder=1)
            for s, e in byp[p]:
                ax.axvspan(s, e + 1, color=SEGMENT["J"], alpha=0.30, lw=0, zorder=0)
                ax.annotate(f"{calls[(p,s,e)]} calls\n{len(donors[(p,s,e)])} donors",
                            ((s + e) / 2, ax.get_ylim()[1] * 0.98),
                            ha="center", va="top", fontsize=6.4,
                            color=SEGMENT["J"], fontweight="bold")
            for x, lab in landmarks(seq)[0]:
                cdr2 = lab == "FR2|CDR2"
                ax.axvline(x, color=NO if cdr2 else INK, lw=1.2 if cdr2 else 1.0,
                           ls="-" if cdr2 else "--", zorder=4)
                ax.annotate(lab, (x, 0), textcoords="offset points",
                            xytext=(2 if cdr2 else -2, 1), fontsize=6,
                            color=NO if cdr2 else INK, va="bottom",
                            ha="left" if cdr2 else "right")
            ax.set_xlim(0, XMAX)
            ax.set_ylabel("divergence", fontsize=7)
            ax.tick_params(labelsize=7)
            ax.set_title(f"{name}  parent {short(p)}   "
                         f"({len(seq)} bp; mean divergence {dv[:lim].mean():.2f})",
                         fontsize=8.4, loc="left", x=0.0, color=col)
            for sp in ("top", "right"):
                ax.spines[sp].set_visible(False)

    axes[-1].set_xlabel("position in the parent V gene (bp)", fontsize=9)
    handles = [
        Patch(facecolor=SEGMENT["J"], alpha=0.35, label="conversion window"),
        Line2D([], [], color=NO, lw=1.5, label="end of FR2 (L-E-W) → CDR2 starts"),
        Line2D([], [], color=INK, lw=1.2, ls="--", label="start of FR2 (W-x-R-Q) → CDR1 ends"),
        Line2D([], [], color=GREY_DARK, lw=1, ls=":", label="gene mean divergence"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=7.8,
               frameon=True, framealpha=0.95, bbox_to_anchor=(0.5, 0.002))
    fig.tight_layout(rect=(0, 0.028 + 0.2 / max(1, nrow), 1, 1))
    save_figure(fig, args.out_figure)
    print(f"{nrow} parent genes plotted", file=sys.stderr)


if __name__ == "__main__":
    main()
