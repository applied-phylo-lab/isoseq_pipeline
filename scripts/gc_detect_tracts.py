"""
Detect gene conversion tracts, scoring each one against an explicit
somatic-hypermutation null and recording whether its donor was still on the
chromosome after V(D)J recombination.

How a tract is called
---------------------
For a transcript T assigned to parent germline gene F, and each candidate
donor D, everything is projected onto F's coordinates.  Only *informative*
positions matter -- those where F and D actually differ, because a position
where F and D agree carries no evidence about which one was the template.

At each informative position the transcript either

    supports the donor   T == D
    contradicts it       T == F
    is unexplained       T is neither (independent mutation or error)

A tract is a maximal run of consecutive informative positions that all support
the same donor.  Runs are broken by a contradicting or unexplained position, so
a called tract is internally consistent by construction.

Scoring against the SHM null
----------------------------
The null is that the differences are independent point mutations.  Conditioned
on the transcript differing from the parent at a position, the chance it
happens to carry the donor's specific base is about 1/3, so a run of m
supporting positions has

    p_raw = 3^-m

Multiplying by the number of donors searched gives a Bonferroni-corrected
p-value.  This is the step BrepConvert omits: it accepts any stretch >=3bp of
clustered mismatches and then asks BLAT to name a donor, which with a hundred
donors essentially always succeeds.  Requiring m supporting positions with
p_corrected below a threshold is a much higher bar -- with ~130 donors a tract
needs roughly m>=8 to clear 0.05.

Donor pool
----------
Each tract is tagged with whether its donor survived recombination
(gc_donor_pool.py).  Under cis-acting conversion, a tract whose donor lay in
the deleted interval is impossible, so those calls measure the false positive
rate directly -- no simulation needed.  Detection itself is run against ALL
donors regardless, so the impossible ones stay observable rather than being
filtered away.
"""
import argparse
import sys
from collections import defaultdict

from Bio import Align

from gc_lib import read_fasta, parse_paf, projected_query, parse_gene_names


def align_donor_to_parent(parent_seq, donor_seq, aligner):
    """
    Return a list, one entry per parent position, holding the donor base
    aligned there (or None where the donor is gapped/absent).
    """
    aln = aligner.align(parent_seq, donor_seq)[0]
    out = [None] * len(parent_seq)
    # aligned blocks: pairs of (parent_range, donor_range)
    for (ps, pe), (ds, de) in zip(aln.aligned[0], aln.aligned[1]):
        for off in range(pe - ps):
            out[ps + off] = donor_seq[ds + off]
    return out


def find_tracts(parent_seq, donor_proj, tx_proj, min_informative, max_gap_bp):
    """
    Walk informative positions and emit maximal pure-support runs.

    A run is broken by a contradicting/unexplained informative position AND by a
    gap of more than max_gap_bp between consecutive supporting positions.  The
    gap rule matters: a conversion tract is a contiguous stretch of
    donor-derived sequence, so two supporting positions 50bp apart are not
    evidence of one tract, they are two unrelated differences that happen to
    suit the same donor.  Without it a "tract" degenerates into "this window
    resembles donor D overall", which is a statement about gene similarity
    rather than about conversion, and it is exactly what makes a detector pick
    donors that recombination had already deleted.

    Returns list of dicts with the parent-coordinate span and the counts.
    """
    informative = []
    for i, (pb, db, tb) in enumerate(zip(parent_seq, donor_proj, tx_proj)):
        if db is None or tb is None or tb == "-" or db == pb:
            continue
        if tb == db:
            state = "support"
        elif tb == pb:
            state = "contradict"
        else:
            state = "other"
        informative.append((i, state))

    def flush(run, out):
        if len(run) >= min_informative:
            out.append({"start": run[0], "end": run[-1], "n_support": len(run)})

    tracts, run = [], []
    for pos, state in informative + [(None, "break")]:
        if state == "support":
            if run and (pos - run[-1]) > max_gap_bp:
                flush(run, tracts)      # too far apart to be one tract
                run = []
            run.append(pos)
            continue
        flush(run, tracts)
        run = []
    return tracts, informative


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--paf", required=True, help="Detailed PAF for this locus")
    p.add_argument("--vgene-fasta", required=True)
    p.add_argument("--assignments", required=True,
                   help="transcript_assignments.tsv from gc_call_functional_genes.py")
    p.add_argument("--donor-pool", required=True)
    p.add_argument("--locus", required=True)
    p.add_argument("--min-informative", type=int, default=3,
                   help="Minimum donor-supporting informative positions in a tract "
                        "(default: 3; report everything and filter on p later)")
    p.add_argument("--max-gap-bp", type=int, default=20,
                   help="Maximum distance between consecutive donor-supporting "
                        "positions within one tract (default: 20). Larger gaps "
                        "split the run, since a conversion tract is a contiguous "
                        "stretch of donor sequence, not diffuse similarity.")
    p.add_argument("--p-threshold", type=float, default=0.05,
                   help="Threshold for calling a tract significant: a corrected "
                        "p-value under --correction bonferroni, or a q-value "
                        "under bh (default: 0.05)")
    p.add_argument("--min-support-significant", type=int, default=4,
                   help="Support cutoff used when --correction support "
                        "(default: 4, the value at which the permutation null "
                        "yields zero tracts in both loci)")
    p.add_argument("--correction",
                   choices=("bonferroni", "bh", "none", "support"),
                   default="bonferroni",
                   help="Multiple-testing correction. 'bonferroni' multiplies by "
                        "the number of donors searched and controls the chance "
                        "that ANY reported tract is spurious -- correct but very "
                        "conservative, especially in IGH where ~162 donors force "
                        "m>=8. 'bh' controls the false discovery RATE among "
                        "reported tracts, which is the appropriate target for a "
                        "catalogue of events rather than a single hypothesis. "
                        "'none' reports raw p (for sensitivity analysis only). "
                        "'support' calls a tract significant on the number of "
                        "supporting positions alone, with the cutoff set by "
                        "--permute. This is the honest option here: 3^-m models "
                        "the chance of matching a donor base but says nothing "
                        "about CONTIGUITY, which is the constraint actually "
                        "doing the work -- so the analytic p-value penalises the "
                        "wrong thing and is wildly over-conservative.")
    p.add_argument("--permute", type=int, default=0, metavar="N",
                   help="Run N permutation replicates to calibrate the threshold "
                        "empirically. Each replicate keeps every transcript's "
                        "differences from its parent -- same number, same "
                        "substituted bases -- but scatters them to random "
                        "positions along the gene, destroying contiguity while "
                        "preserving everything else. Tracts found in the "
                        "permuted data are false by construction, so the ratio "
                        "of permuted to real tracts at a given support cutoff IS "
                        "the false discovery rate. This is the only usable "
                        "control in IGL, where nothing is ever deleted and the "
                        "topology check therefore has no power.")
    p.add_argument("--permute-seed", type=int, default=0)
    p.add_argument("--out-null", help="Per-threshold FDR table from --permute")
    p.add_argument("--out-tracts", required=True)
    p.add_argument("--out-summary", required=True)
    args = p.parse_args()

    vgenes = read_fasta(args.vgene_fasta)
    info = parse_gene_names(vgenes)

    # parent assignment per transcript, restricted to this locus
    parent_of = {}
    with open(args.assignments) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            r = dict(zip(header, line.rstrip("\n").split("\t")))
            if r["locus"] == args.locus and r["best_gene"] in vgenes:
                parent_of[r["transcript"]] = r["best_gene"]

    # allowed donor sets
    allowed_of, mech_of = {}, {}
    with open(args.donor_pool) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            r = dict(zip(header, line.rstrip("\n").split("\t")))
            allowed_of[r["rearranged_gene"]] = (
                set() if r["allowed_donors"] == "NONE"
                else set(r["allowed_donors"].split(",")))
            mech_of[r["rearranged_gene"]] = r["mechanism"]

    # transcript projected onto its parent's coordinates
    tx_proj = {}
    for rec in parse_paf(args.paf):
        parent = parent_of.get(rec.query)
        if parent is None or rec.target != parent:
            continue
        proj = projected_query(rec, vgenes[parent])
        if proj is None:
            continue
        prev = tx_proj.get(rec.query)
        if prev is None or sum(b is not None for b in proj) > sum(b is not None for b in prev):
            tx_proj[rec.query] = proj

    aligner = Align.PairwiseAligner(mode="global", open_gap_score=-10,
                                    extend_gap_score=-0.5, match_score=2,
                                    mismatch_score=-1)

    # donor projections, cached per (parent, donor)
    donor_cache = {}

    def get_donor_proj(parent, donor):
        key = (parent, donor)
        if key not in donor_cache:
            donor_cache[key] = align_donor_to_parent(vgenes[parent], vgenes[donor], aligner)
        return donor_cache[key]

    n_donors = len(vgenes) - 1
    rows = []
    for tid, proj in sorted(tx_proj.items()):
        parent = parent_of[tid]
        allowed = allowed_of.get(parent, set())
        mech = mech_of.get(parent, "unknown")
        for donor in vgenes:
            if donor == parent:
                continue
            dproj = get_donor_proj(parent, donor)
            tracts, _ = find_tracts(vgenes[parent], dproj, proj,
                                    args.min_informative, args.max_gap_bp)
            for t in tracts:
                m = t["n_support"]
                p_raw = 3.0 ** (-m)
                rows.append({
                    "transcript": tid, "parent": parent, "donor": donor,
                    "mechanism": mech,
                    "donor_allowed": donor in allowed,
                    "start": t["start"], "end": t["end"],
                    "span_bp": t["end"] - t["start"] + 1,
                    "n_support": m, "p_raw": p_raw,
                    "p_corrected": min(1.0, p_raw * n_donors),
                })

    # ── permutation calibration ──────────────────────────────────────────────
    if args.permute:
        import random

        rng = random.Random(args.permute_seed)
        real_by_m = defaultdict(int)
        for r in rows:
            real_by_m[r["n_support"]] += 1
        null_by_m = defaultdict(int)

        for rep in range(args.permute):
            for tid, proj in tx_proj.items():
                parent = parent_of[tid]
                pseq = vgenes[parent]
                # collect this transcript's differences, then relocate them
                diffs, sites = [], []
                for i, (pb, tb) in enumerate(zip(pseq, proj)):
                    if tb is None or tb == "-":
                        continue
                    sites.append(i)
                    if tb != pb:
                        diffs.append(tb)
                if not diffs or len(sites) <= len(diffs):
                    continue
                shuffled = list(proj)
                for i in sites:
                    shuffled[i] = pseq[i]           # reset to parent
                for i, b in zip(rng.sample(sites, len(diffs)), diffs):
                    shuffled[i] = b                 # scatter the same substitutions
                for donor in vgenes:
                    if donor == parent:
                        continue
                    dproj = get_donor_proj(parent, donor)
                    tracts, _ = find_tracts(pseq, dproj, shuffled,
                                            args.min_informative, args.max_gap_bp)
                    for t in tracts:
                        null_by_m[t["n_support"]] += 1

        if args.out_null:
            with open(args.out_null, "w") as fh:
                fh.write("min_support\treal_tracts\tnull_tracts_per_replicate\t"
                         "estimated_fdr\n")
                for cut in range(args.min_informative, 25):
                    real = sum(v for k, v in real_by_m.items() if k >= cut)
                    null = sum(v for k, v in null_by_m.items() if k >= cut) / args.permute
                    if real == 0 and null == 0:
                        break
                    fdr = min(1.0, null / real) if real else float("nan")
                    fh.write(f"{cut}\t{real}\t{null:.2f}\t{fdr:.4f}\n")
        print(f"  permutation: {args.permute} replicates", file=sys.stderr)
        for cut in range(args.min_informative, 16):
            real = sum(v for k, v in real_by_m.items() if k >= cut)
            null = sum(v for k, v in null_by_m.items() if k >= cut) / args.permute
            if real == 0 and null == 0:
                break
            fdr = min(1.0, null / real) if real else float("nan")
            print(f"    m>={cut:<3} real={real:<6} null={null:<8.2f} FDR={fdr:.3f}",
                  file=sys.stderr)

    # ── multiple testing ─────────────────────────────────────────────────────
    # Bonferroni asks "could this single tract be spurious given I scanned every
    # donor?" -- the right question for one hypothesis, far too strict for a
    # catalogue.  BH instead fixes the expected FRACTION of reported tracts that
    # are wrong, which is what a table of events actually needs to be trustworthy.
    n_tests = len(rows)
    if args.correction == "bh":
        order = sorted(range(n_tests), key=lambda i: rows[i]["p_raw"])
        running = 1.0
        # step-up: walk from the largest p downwards, keeping the running minimum
        for rank, i in enumerate(reversed(order), start=1):
            k = n_tests - rank + 1
            running = min(running, rows[i]["p_raw"] * n_tests / k)
            rows[i]["q_value"] = min(1.0, running)
        for r in rows:
            r["significant"] = r["q_value"] <= args.p_threshold
    elif args.correction == "support":
        for r in rows:
            r["q_value"] = float("nan")
            r["significant"] = r["n_support"] >= args.min_support_significant
    else:
        key = "p_corrected" if args.correction == "bonferroni" else "p_raw"
        for r in rows:
            r["q_value"] = float("nan")
            r["significant"] = r[key] <= args.p_threshold

    # ── competing donors for the same stretch of transcript ──────────────────
    # Several donors can explain one tract when they share the diagnostic bases.
    # Those are alternative explanations of ONE event, not several events, so
    # counting them all inflates the donor network.  Rank them by support and
    # flag the best as primary; figures can then draw the primary solid and the
    # alternatives greyed out instead of silently over-counting.
    by_region = defaultdict(list)
    for r in rows:
        if r["significant"]:
            by_region[(r["transcript"], r["start"], r["end"])].append(r)
    for group in by_region.values():
        group.sort(key=lambda r: (-r["n_support"], r["donor"]))
        for rank, r in enumerate(group, start=1):
            r["donor_rank"] = rank
            r["n_candidate_donors"] = len(group)
            r["primary_donor"] = rank == 1
    for r in rows:
        r.setdefault("donor_rank", 0)
        r.setdefault("n_candidate_donors", 0)
        r.setdefault("primary_donor", False)

    with open(args.out_tracts, "w") as fh:
        fh.write("transcript\tparent\tdonor\tmechanism\tdonor_allowed\tstart\tend"
                 "\tspan_bp\tn_support\tp_raw\tp_corrected\tq_value\tsignificant"
                 "\tdonor_rank\tn_candidate_donors\tprimary_donor\n")
        for r in sorted(rows, key=lambda r: (r["transcript"], r["start"], -r["n_support"])):
            fh.write(f"{r['transcript']}\t{r['parent']}\t{r['donor']}\t{r['mechanism']}"
                     f"\t{r['donor_allowed']}\t{r['start']}\t{r['end']}\t{r['span_bp']}"
                     f"\t{r['n_support']}\t{r['p_raw']:.3g}\t{r['p_corrected']:.3g}"
                     f"\t{r['q_value']:.3g}\t{r['significant']}"
                     f"\t{r['donor_rank']}\t{r['n_candidate_donors']}"
                     f"\t{r['primary_donor']}\n")

    sig = [r for r in rows if r["significant"]]
    tx_with_sig = {r["transcript"] for r in sig}
    # deletional transcripts only: that is where "impossible" is meaningful
    del_sig = [r for r in sig if r["mechanism"] == "deletion"]
    del_bad = [r for r in del_sig if not r["donor_allowed"]]

    with open(args.out_summary, "w") as fh:
        fh.write("metric\tvalue\n")
        fh.write(f"locus\t{args.locus}\n")
        fh.write(f"transcripts_analysed\t{len(tx_proj)}\n")
        fh.write(f"candidate_tracts\t{len(rows)}\n")
        fh.write(f"significant_tracts\t{len(sig)}\n")
        fh.write(f"transcripts_with_significant_tract\t{len(tx_with_sig)}\n")
        fh.write(f"significant_tracts_deletional\t{len(del_sig)}\n")
        fh.write(f"significant_tracts_impossible_donor\t{len(del_bad)}\n")
        if del_sig:
            fh.write(f"impossible_donor_fraction\t{len(del_bad)/len(del_sig):.4f}\n")
        fh.write(f"p_threshold\t{args.p_threshold}\n")
        fh.write(f"min_informative\t{args.min_informative}\n")
        fh.write(f"max_gap_bp\t{args.max_gap_bp}\n")

    print(
        f"{args.locus}: {len(tx_proj)} transcripts, {len(rows)} candidate tracts\n"
        f"  significant (p_corr<={args.p_threshold}): {len(sig)} tracts "
        f"in {len(tx_with_sig)} transcripts\n"
        f"  of deletional significant tracts: {len(del_bad)}/{len(del_sig)} "
        f"use an IMPOSSIBLE (deleted) donor"
        + (f" = {len(del_bad)/len(del_sig)*100:.1f}% false-positive floor"
           if del_sig else ""),
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
