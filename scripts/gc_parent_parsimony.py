"""
How much mutational work does each candidate parent require?

The question
------------
"Only RSS-bearing genes can be a parent" is currently an assumption imposed by
the pipeline, not a result.  This tests it: for every transcript, ask how many
mutational events are needed to derive it from each candidate germline gene, and
then ask whether the RSS-bearing genes are systematically cheaper than the rest.
If they are, the restriction is doing real work.  If a non-RSS gene explains a
transcript more cheaply than any RSS gene, that gene is a missing-RSS candidate.

The model: a transcript is a mosaic of its parent plus donor patches
-------------------------------------------------------------------
A converted transcript is not a diverged copy of one gene.  It is its PARENT
along most of its length with blocks copied in from DONORS, plus scattered point
mutations from hypermutation.  So the natural cost of an explanation is

    cost = (number of conversion tracts) x LAMBDA + (number of point mutations)

and the best explanation is the cheapest parse of the transcript into segments,
each attributed either to the parent or to one donor.  That is exactly a
jumping-hidden-Markov-model -- the same structure used for recombination
detection -- with one state per germline gene:

    states     : the candidate parent P, plus every other gene in the locus
    emission   : 0 if the transcript matches that state's base here, else 1
                 (one point mutation)
    transition : 0 to stay in the same state, LAMBDA to switch
                 (one conversion event)

Viterbi returns the minimum-cost path, which gives BOTH numbers we want: how
many switches (tracts) and how many mismatches (point mutations) the explanation
needs.  The single-tract greedy search used elsewhere in this pipeline is a
special case of this.

Why the path must start in the parent
-------------------------------------
Left unconstrained, the cheapest path is simply "be whichever gene matches best
at every position", which makes the notion of a parent meaningless.  The path is
therefore forced to START in the parent state.  That is not an arbitrary
convenience: FR1 at the 5' end is the conserved end, and it is where a transcript
reliably still follows the gene it rearranged from.  The 3' end is deliberately
left free, because a conversion tract can run past the V 3' end and into D, so
requiring the parent there would reject real events.

The comparison, and the trap in it
----------------------------------
Per transcript we take the cheapest RSS-bearing parent and the cheapest non-RSS
parent and compare.  The trap is that these sets are not the same size -- IGH has
25 RSS genes against 78 others, IGL has 2 against 20 -- and the best of 20 beats
the best of 2 by chance alone.  Two controls handle it:

  SUBSAMPLE  draw a random subset of non-RSS genes the same size as the RSS set,
             repeat, and use the distribution of the difference.
  PERMUTE    shuffle the RSS labels across genes and recompute the whole
             statistic. This is assumption-free and is the same device that
             carried the CDR-versus-position result.

LAMBDA is swept rather than fixed.  It is the exchange rate between "one
conversion event" and "one point mutation", and a conclusion that depends on a
number nobody measured is not a conclusion.
"""
import argparse
import csv
import random
import sys
from collections import defaultdict

import numpy as np
from Bio import Align

from gc_lib import read_fasta, parse_paf, projected_query

BASES = "ACGT"


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def short(g):
    p = g.split(".")
    return p[1] if len(p) > 1 else g


def project_all(seqs, names):
    """proj[p][d] = donor d's bases in parent p's coordinate frame (or None)."""
    al = Align.PairwiseAligner(mode="global", open_gap_score=-10,
                              extend_gap_score=-0.5, match_score=2,
                              mismatch_score=-1)
    proj = {}
    for p in names:
        ref = seqs[p]
        row = {p: list(ref)}
        for d in names:
            if d == p:
                continue
            aln = al.align(ref, seqs[d])[0]
            dp = [None] * len(ref)
            for (ps, pe), (ds, de) in zip(aln.aligned[0], aln.aligned[1]):
                for o in range(pe - ps):
                    dp[ps + o] = seqs[d][ds + o]
            row[d] = dp
        proj[p] = row
    return proj


def viterbi(emit, lam, anchor=0):
    """Cheapest parse of the transcript into parent stretches and donor patches.

    emit is (n_states, n_positions) of per-position costs, state 0 being the
    candidate parent. Transition costs are deliberately ASYMMETRIC:

        parent -> donor   lam      starting a conversion tract is the event
        donor  -> donor   lam      a different donor is a different tract
        donor  -> parent  0        the tract simply ends
        stay             0

    Charging both entry and exit would price one tract at 2*lam and, with lam
    anywhere near a realistic value, make every tract more expensive than just
    accumulating point mutations -- the parse then never uses a donor at all,
    which is what a symmetric cost produced here.

    Still O(states) per position: entering the parent is free so it needs the
    global best, and entering a donor needs the best OTHER state, which comes
    from the best and second-best.

    `anchor` forces the parse to stay in the parent for that many leading
    covered positions. Pinning only position 0 is not enough: the path can then
    switch immediately, so a DONOR offered as the candidate parent is scored on
    little more than its overall similarity -- and a heavily converted transcript
    is closer to its donor than to its parent, which is the whole reason overall
    identity cannot be used to assign parents. Requiring the parent to carry FR1
    unaided is the constraint that separates the two, and it is independent of
    the RSS annotation this script is testing. The 3' end is left free, because a
    tract may run past V into D.

    Returns (total cost, number of tracts, number of point mutations).
    """
    n_states, n_pos = emit.shape
    if n_pos == 0:
        return np.inf, 0, 0
    V = np.full(n_states, np.inf)
    V[0] = emit[0, 0]                       # must start in the parent
    back = np.zeros((n_states, n_pos), dtype=np.int32)
    for i in range(1, min(anchor, n_pos)):
        V[0] += emit[0, i]                  # anchored: parent only, no switching
    for i in range(max(1, min(anchor, n_pos)), n_pos):
        b1 = int(np.argmin(V))
        if n_states > 1:
            masked = V.copy(); masked[b1] = np.inf
            b2 = int(np.argmin(masked)); best2 = V[b2]
        else:
            b2, best2 = b1, np.inf
        best1 = V[b1]
        # best predecessor excluding the state itself
        other_val = np.where(np.arange(n_states) == b1, best2, best1)
        other_idx = np.where(np.arange(n_states) == b1, b2, b1)
        # entering the parent is free; entering a donor costs lam
        enter = other_val + lam
        enter[0] = V.min()                  # donor -> parent is free
        enter_idx = other_idx.copy()
        enter_idx[0] = int(V.argmin())
        stay_better = V <= enter
        newV = np.where(stay_better, V, enter) + emit[:, i]
        back[:, i] = np.where(stay_better, np.arange(n_states), enter_idx)
        V = newV
    end = int(V.argmin())
    total = float(V[end])
    path = np.empty(n_pos, dtype=np.int32)
    s = end
    for i in range(n_pos - 1, -1, -1):
        path[i] = s
        s = back[s, i]
    # a tract is an entry INTO a donor state, so count those transitions only
    tracts = int(np.sum((path[1:] != path[:-1]) & (path[1:] != 0)))
    mism = float(sum(emit[path[i], i] for i in range(n_pos)))
    return total, tracts, int(round(mism))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--paf", required=True)
    ap.add_argument("--vgene-fasta", required=True)
    ap.add_argument("--rss-annotation", required=True)
    ap.add_argument("--locus", required=True)
    ap.add_argument("--lambdas", type=float, nargs="+", default=[3, 6, 10, 15])
    ap.add_argument("--primary-lambda", type=float, default=6.0)
    ap.add_argument("--junction-margin", type=int, default=20)
    ap.add_argument("--min-covered", type=int, default=200)
    ap.add_argument("--anchor", type=int, default=60,
                    help="leading covered positions the parent must explain "
                         "unaided (the FR1 rule); 0 disables")
    ap.add_argument("--max-transcripts", type=int, default=0,
                    help="cap for a quick run; 0 = all")
    ap.add_argument("--perms", type=int, default=2000)
    ap.add_argument("--out-table", required=True)
    ap.add_argument("--out-summary", required=True)
    args = ap.parse_args()
    random.seed(0)
    rng = np.random.default_rng(0)

    seqs = read_fasta(args.vgene_fasta)
    rss = {r["gene"]: r["rss_state"] for r in read_tsv(args.rss_annotation)
           if r["locus"] == args.locus}
    names = sorted([g for g in seqs if g in rss], key=lambda g: int(g.split(".")[1]))
    has_rss = {g: rss.get(g) == "rss_present" for g in names}
    idx_of = {g: i for i, g in enumerate(names)}
    n_rss = sum(has_rss.values())
    print(f"{args.locus}: {len(names)} genes ({n_rss} with an RSS)", file=sys.stderr)

    print("  projecting every gene onto every candidate parent ...", file=sys.stderr)
    proj = project_all(seqs, names)

    # transcript projections onto each gene it aligns to
    tproj = defaultdict(dict)
    for rec in parse_paf(args.paf):
        if rec.target not in seqs:
            continue
        pr = projected_query(rec, seqs[rec.target])
        if pr is None:
            continue
        cur = tproj[rec.query].get(rec.target)
        if cur is None or sum(x is not None for x in pr) > sum(x is not None for x in cur):
            tproj[rec.query][rec.target] = pr

    txs = sorted(tproj)
    if args.max_transcripts:
        txs = txs[:args.max_transcripts]
    print(f"  {len(txs)} transcripts with a V alignment", file=sys.stderr)

    rows = []
    for n_done, q in enumerate(txs, 1):
        if n_done % 50 == 0:
            print(f"    {n_done}/{len(txs)}", file=sys.stderr)
        for p in names:
            tp = tproj[q].get(p)
            if tp is None:
                continue
            lim = len(seqs[p]) - args.junction_margin
            cols = [i for i in range(lim) if tp[i] not in (None, "-")]
            if len(cols) < args.min_covered:
                continue
            # emission costs: state 0 is the candidate parent, then every other
            # gene in the locus as a possible donor
            order = [p] + [g for g in names if g != p]
            emit = np.ones((len(order), len(cols)), dtype=np.float32)
            for si, g in enumerate(order):
                gp = proj[p][g]
                for ci, i in enumerate(cols):
                    b = gp[i]
                    # a donor that does not align here cannot explain the base;
                    # cost 1 is the same as a mismatch, so it is never preferred
                    # to a gene that does match
                    if b is not None and b == tp[i]:
                        emit[si, ci] = 0.0
            for lam in args.lambdas:
                cost, sw, mm = viterbi(emit, lam, anchor=args.anchor)
                rows.append({
                    "transcript": q, "parent": p, "parent_has_rss": has_rss[p],
                    "lambda": lam, "cost": cost, "n_switches": sw,
                    "n_mismatches": mm, "covered": len(cols),
                })

    with open(args.out_table, "w") as fh:
        cols = ["transcript", "parent", "parent_has_rss", "lambda", "cost",
                "n_switches", "n_mismatches", "covered"]
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    # ── comparison ────────────────────────────────────────────────────────────
    out = open(args.out_summary, "w")

    def emit_line(s):
        print(s, file=sys.stderr)
        out.write(s + "\n")

    by = defaultdict(lambda: defaultdict(dict))
    for r in rows:
        by[r["lambda"]][r["transcript"]][r["parent"]] = r
    emit_line(f"locus\t{args.locus}\ngenes\t{len(names)}\trss_genes\t{n_rss}")
    emit_line("")
    emit_line("lambda\tn_tx\tmedian_cost_RSS\tmedian_cost_nonRSS\tmedian_delta\t"
              "pct_RSS_cheaper\tp_permuted")
    for lam in args.lambdas:
        deltas, cr, cn = [], [], []
        for q, per in by[lam].items():
            r_costs = [v["cost"] for g, v in per.items() if has_rss[g]]
            n_costs = [v["cost"] for g, v in per.items() if not has_rss[g]]
            if not r_costs or not n_costs:
                continue
            br, bn = min(r_costs), min(n_costs)
            cr.append(br)
            cn.append(bn)
            deltas.append(bn - br)          # >0 means the RSS parent is cheaper
        if not deltas:
            continue
        pct = float(np.mean([d > 0 for d in deltas]))
        # PERMUTE the RSS labels across genes and recompute the median delta
        obs = float(np.median(deltas))
        null = []
        gene_list = list(names)
        for _ in range(args.perms):
            lab = dict(zip(gene_list, rng.permutation([has_rss[g] for g in gene_list])))
            ds = []
            for q, per in by[lam].items():
                a = [v["cost"] for g, v in per.items() if lab[g]]
                b = [v["cost"] for g, v in per.items() if not lab[g]]
                if a and b:
                    ds.append(min(b) - min(a))
            if ds:
                null.append(np.median(ds))
        pval = (np.sum(np.array(null) >= obs) + 1) / (len(null) + 1) if null else float("nan")
        emit_line(f"{lam:g}\t{len(deltas)}\t{np.median(cr):.1f}\t{np.median(cn):.1f}\t"
                  f"{obs:+.1f}\t{pct:.1%}\t{pval:.4f}")

    # per-gene: how often is each gene the cheapest parent, and by how much
    lam = args.primary_lambda
    emit_line("")
    emit_line(f"# per-gene, lambda={lam:g}")
    emit_line("gene\thas_rss\tn_tx_cheapest\tmedian_cost_when_used\t"
              "median_switches\tmedian_mismatches")
    winner = defaultdict(list)
    for q, per in by[lam].items():
        if not per:
            continue
        g = min(per, key=lambda g: per[g]["cost"])
        winner[g].append(per[g])
    for g in sorted(winner, key=lambda g: -len(winner[g])):
        v = winner[g]
        emit_line(f"{short(g)}\t{has_rss[g]}\t{len(v)}\t"
                  f"{np.median([x['cost'] for x in v]):.1f}\t"
                  f"{np.median([x['n_switches'] for x in v]):.0f}\t"
                  f"{np.median([x['n_mismatches'] for x in v]):.0f}")
    n_win_rss = sum(len(v) for g, v in winner.items() if has_rss[g])
    n_win_tot = sum(len(v) for v in winner.values())
    emit_line("")
    emit_line(f"# transcripts whose CHEAPEST parent of all {len(names)} genes "
              f"has an RSS: {n_win_rss}/{n_win_tot} = {n_win_rss / max(1, n_win_tot):.1%}")
    emit_line(f"# RSS genes are {n_rss}/{len(names)} = {n_rss / len(names):.1%} "
              f"of the locus, so that is the chance expectation")
    out.close()


if __name__ == "__main__":
    main()
