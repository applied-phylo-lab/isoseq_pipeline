"""
Quantify what changes when the germline reference comes from a DIFFERENT
individual than the sequenced bird.

The same transcripts are scored twice -- once against each germline set -- so
every difference reported is caused by the reference alone.  Four consequences
are measured, in increasing order of how badly they corrupt downstream biology:

  1. apparent divergence from germline.  A mismatched reference inflates this,
     and the inflation is read as somatic hypermutation or gene conversion that
     never happened.
  2. which genes are callable as "used".  Requiring near-perfect identity to
     call a gene expressed fails when the reference lacks that individual's
     alleles.
  3. gene usage profile.  If the ranking of V genes changes, any statement
     about repertoire composition changes with it.
  4. PARENT ASSIGNMENT.  This is the one that propagates: a transcript assigned
     to the wrong germline gene makes every downstream tract, donor and
     mechanism call wrong, because they are all defined relative to the parent.

Orthologues between the two germline sets are established by best reciprocal
alignment, so "assigned to a different gene" means genuinely different, not
merely differently named.
"""
import argparse
import statistics
import sys
from collections import Counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from Bio import Align

from gc_lib import read_fasta, parse_paf, projected_query

from gc_palette import save_figure, NO as C_A, YES as C_B, GREY, BOTH


def best_per_transcript(paf, vgenes, margin, min_cov):
    """Highest-identity germline gene per transcript, scored off the junction."""
    best = {}
    for rec in parse_paf(paf):
        vs = vgenes.get(rec.target)
        if vs is None:
            continue
        proj = projected_query(rec, vs)
        if proj is None:
            continue
        limit = max(0, rec.tlen - margin)
        cov = match = 0
        for i in range(limit):
            b = proj[i]
            if b is None or b == "-":
                continue
            cov += 1
            match += (b == vs[i])
        if cov < min_cov:
            continue
        idt = match / cov
        if idt > best.get(rec.query, (0.0, None))[0]:
            best[rec.query] = (idt, rec.target)
    return best


def orthologues(a_seqs, b_seqs):
    """Best reciprocal alignment pairing between two germline sets."""
    aligner = Align.PairwiseAligner(mode="global", open_gap_score=-10,
                                    extend_gap_score=-0.5, match_score=2,
                                    mismatch_score=-1)

    def ident(x, y):
        aln = aligner.align(x, y)[0]
        s1, s2 = str(aln[0]), str(aln[1])
        m = sum(1 for p, q in zip(s1, s2) if p == q and p != "-")
        L = sum(1 for p, q in zip(s1, s2) if p != "-" and q != "-")
        return m / L if L else 0.0

    a2b, b2a = {}, {}
    for a, sa in a_seqs.items():
        scores = [(ident(sa, sb), b) for b, sb in b_seqs.items()]
        a2b[a] = max(scores)[1] if scores else None
    for b, sb in b_seqs.items():
        scores = [(ident(sb, sa), a) for a, sa in a_seqs.items()]
        b2a[b] = max(scores)[1] if scores else None
    # keep only reciprocal pairs
    return {a: b for a, b in a2b.items() if b is not None and b2a.get(b) == a}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--paf-a", required=True, help="transcripts vs germline A")
    ap.add_argument("--paf-b", required=True, help="transcripts vs germline B")
    ap.add_argument("--germline-a", required=True)
    ap.add_argument("--germline-b", required=True)
    ap.add_argument("--label-a", default="different individual")
    ap.add_argument("--label-b", default="matched individual")
    ap.add_argument("--locus", required=True)
    ap.add_argument("--junction-margin", type=int, default=20)
    ap.add_argument("--min-covered-bp", type=int, default=200)
    ap.add_argument("--used-identity", type=float, default=0.98,
                    help="identity needed to call a gene expression-confirmed")
    ap.add_argument("--restrict-to",
                    help="File of transcript IDs, one per line. Scores only "
                         "these. Used to force every reference in a multi-way "
                         "comparison onto the SAME transcript set: the coverage "
                         "floor is applied per reference, so a borderline "
                         "transcript can clear it against one germline and not "
                         "another, leaving rows with slightly different n and "
                         "undermining the claim that only the reference changed.")
    ap.add_argument("--out-figure", required=True)
    ap.add_argument("--out-stats", required=True)
    args = ap.parse_args()

    a_seqs, b_seqs = read_fasta(args.germline_a), read_fasta(args.germline_b)
    best_a = best_per_transcript(args.paf_a, a_seqs, args.junction_margin,
                                 args.min_covered_bp)
    best_b = best_per_transcript(args.paf_b, b_seqs, args.junction_margin,
                                 args.min_covered_bp)
    shared = set(best_a) & set(best_b)
    if args.restrict_to:
        with open(args.restrict_to) as fh:
            keep = {ln.strip() for ln in fh if ln.strip()}
        shared &= keep
    shared = sorted(shared)
    if not shared:
        raise SystemExit("no transcripts scored against both germline sets")

    ida = [best_a[t][0] for t in shared]
    idb = [best_b[t][0] for t in shared]

    ortho = orthologues(a_seqs, b_seqs)
    concordant = discordant = unmappable = 0
    for t in shared:
        ga, gb = best_a[t][1], best_b[t][1]
        mapped = ortho.get(ga)
        if mapped is None:
            unmappable += 1
        elif mapped == gb:
            concordant += 1
        else:
            discordant += 1

    used_a = {g for g in a_seqs
              if max([best_a[t][0] for t in shared if best_a[t][1] == g], default=0)
              >= args.used_identity}
    used_b = {g for g in b_seqs
              if max([best_b[t][0] for t in shared if best_b[t][1] == g], default=0)
              >= args.used_identity}

    usage_a = Counter(best_a[t][1] for t in shared)
    usage_b = Counter(best_b[t][1] for t in shared)

    # ── figure ────────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 4, figsize=(17, 4.1))

    ax = axes[0]
    bins = np.arange(.85, 1.005, .01)
    ax.hist(ida, bins=bins, alpha=.62, color=C_A, edgecolor="black", lw=.4,
            label=f"{args.label_a}\nmedian {statistics.median(ida):.4f}")
    ax.hist(idb, bins=bins, alpha=.85, color=C_B, edgecolor="black", lw=.5,
            label=f"{args.label_b}\nmedian {statistics.median(idb):.4f}")
    ax.set_xlabel("identity to best germline V")
    ax.set_ylabel("transcripts")
    ax.legend(fontsize=7.5)
    ax.set_title("A  apparent divergence\nmismatched reference inflates it",
                 fontsize=9.5, fontweight="bold", loc="left")

    ax = axes[1]
    delta = [b - a for a, b in zip(ida, idb)]
    ax.hist(delta, bins=30, color=GREY, edgecolor="black", lw=.4)
    ax.axvline(0, color="black", lw=1)
    med = statistics.median(delta)
    ax.axvline(med, color=C_B, lw=2, label=f"median {med:+.4f}")
    ax.set_xlabel("identity gain with matched reference")
    ax.set_ylabel("transcripts")
    ax.legend(fontsize=8)
    ax.set_title("B  per-transcript gain\n"
                 f"{sum(1 for d in delta if d > 0)}/{len(delta)} improve",
                 fontsize=9.5, fontweight="bold", loc="left")

    ax = axes[2]
    tot = concordant + discordant + unmappable
    ax.bar(["same gene\n(agree)", "DIFFERENT gene\n(wrong parent)",
            "no equivalent gene\n(unanswerable)"],
           [concordant, discordant, unmappable],
           color=[BOTH, C_A, GREY], edgecolor="black", lw=.6)
    for i, v in enumerate([concordant, discordant, unmappable]):
        ax.text(i, v, f"\n{v}\n({v/tot*100:.0f}%)", ha="center", va="bottom",
                fontsize=8.5)
    ax.set_ylabel("transcripts")
    ax.set_ylim(0, max(concordant, discordant, unmappable) * 1.42)
    ax.set_title("C  did the transcript get the same parent?\n"
                 "genes paired between references by best reciprocal alignment",
                 fontsize=9.5, fontweight="bold", loc="left")
    ax.tick_params(axis="x", labelsize=7.5)

    ax = axes[3]
    ax.bar([f"look expressed\n(≥{args.used_identity:g} identity\nto a transcript)",
            "attract ≥1\ntranscript"],
           [len(used_a), len(usage_a)], width=.38, align="edge",
           color=C_A, edgecolor="black", lw=.6, label=args.label_a)
    ax.bar([f"look expressed\n(≥{args.used_identity:g} identity\nto a transcript)",
            "attract ≥1\ntranscript"],
           [len(used_b), len(usage_b)], width=-.38, align="edge",
           color=C_B, edgecolor="black", lw=.6, label=args.label_b)
    for i, (va, vb) in enumerate([(len(used_a), len(used_b)),
                                  (len(usage_a), len(usage_b))]):
        ax.text(i + .19, va, f"{va}", ha="center", va="bottom", fontsize=8.5)
        ax.text(i - .19, vb, f"{vb}", ha="center", va="bottom", fontsize=8.5)
    ax.set_ylabel("V genes")
    ax.legend(fontsize=7.5)
    ax.set_title("D  how many V genes look real\n"
                 "the reference changes the answer",
                 fontsize=9.5, fontweight="bold", loc="left")
    ax.tick_params(axis="x", labelsize=7.5)

    for a in axes:
        for s in ("top", "right"):
            a.spines[s].set_visible(False)
    fig.suptitle(f"{args.locus} — effect of using a germline from a different individual "
                 f"({len(shared)} transcripts, identical in both panels)",
                 fontsize=12, fontweight="bold", y=1.04)
    fig.tight_layout()
    save_figure(fig, args.out_figure)

    with open(args.out_stats, "w") as fh:
        fh.write("metric\tvalue\n")
        fh.write(f"locus\t{args.locus}\n")
        fh.write(f"transcripts_compared\t{len(shared)}\n")
        fh.write(f"median_identity_{args.label_a.replace(' ','_')}\t{statistics.median(ida):.4f}\n")
        fh.write(f"median_identity_{args.label_b.replace(' ','_')}\t{statistics.median(idb):.4f}\n")
        fh.write(f"median_identity_gain\t{med:.4f}\n")
        fh.write(f"transcripts_improved\t{sum(1 for d in delta if d > 0)}\n")
        fh.write(f"parent_same_gene\t{concordant}\n")
        fh.write(f"parent_different_gene\t{discordant}\n")
        fh.write(f"parent_no_orthologue\t{unmappable}\n")
        fh.write(f"parent_discordance_rate\t{discordant/tot:.4f}\n")
        fh.write(f"genes_called_used_a\t{len(used_a)}\n")
        fh.write(f"genes_called_used_b\t{len(used_b)}\n")

    print(f"{args.locus}: {len(shared)} transcripts compared\n"
          f"  median identity {statistics.median(ida):.4f} -> {statistics.median(idb):.4f} "
          f"({med:+.4f})\n"
          f"  parent assignment DIFFERS for {discordant}/{tot} "
          f"({discordant/tot*100:.1f}%)\n"
          f"  genes called used: {len(used_a)} -> {len(used_b)}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
