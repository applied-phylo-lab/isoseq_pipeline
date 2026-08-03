"""
Call which germline V genes are genuinely *used* (rearranged and expressed),
using the transcript evidence rather than the Productive flag in the annotation.

Rationale
---------
The Productive/pseudogene labels in a de novo germline annotation are
predictions.  A V gene that a transcript matches at near-perfect identity over
its whole length must have been rearranged and expressed, which is direct
evidence of functionality.  Nothing else in the analysis should depend on the
annotation's guess.

Why not require a literal 100% match
------------------------------------
Two things stop a transcript from matching its own germline V exactly:

  * the 3' end of V is chewed back and joined to D/J with N additions, so the
    last few codons never match germline.  We therefore score identity only
    over V positions before `--junction-margin` bases from the 3' end.
  * somatic hypermutation and gene conversion introduce real differences, and
    the reference germline may come from a different individual entirely.

So "used" is defined as: at least `--min-transcripts` transcripts pick this
gene as their single best match, and the best of them reaches
`--min-identity` over at least `--min-covered-bp` assessable positions.

Each transcript is assigned to exactly one gene (its argmax), so paralogues do
not all get called functional off the back of one transcript.
"""
import argparse
import statistics
import sys
from collections import defaultdict

from gc_lib import read_fasta, parse_paf, projected_query, parse_gene_names


def score_alignment(rec, vseq, junction_margin):
    """Identity of one transcript against one V gene over the assessable region."""
    proj = projected_query(rec, vseq)
    if proj is None:
        return None
    limit = max(0, rec.tlen - junction_margin)
    covered = matches = 0
    for i in range(limit):
        b = proj[i]
        if b is None or b == "-":
            continue
        covered += 1
        if b == vseq[i]:
            matches += 1
    if covered == 0:
        return None
    return {
        "covered": covered,
        "matches": matches,
        "mismatches": covered - matches,
        "identity": matches / covered,
    }


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pafs", nargs="+", required=True,
                   help="Detailed PAFs (minimap2 --cs -N ...) of transcripts vs V genes")
    p.add_argument("--vgene-fastas", nargs="+", required=True)
    p.add_argument("--junction-margin", type=int, default=20,
                   help="Ignore this many bases at the 3' end of each V gene, "
                        "where V(D)J junction formation destroys germline identity "
                        "(default: 20)")
    p.add_argument("--min-identity", type=float, default=0.98,
                   help="Best-transcript identity needed to call a gene used "
                        "(default: 0.98)")
    p.add_argument("--min-covered-bp", type=int, default=200,
                   help="Minimum assessable positions covered (default: 200)")
    p.add_argument("--min-transcripts", type=int, default=2,
                   help="Independent transcripts that must pick the gene as their "
                        "best match (default: 2)")
    p.add_argument("--candidate-parents",
                   help="TSV from gc_rss_annotation.py. When given, only genes "
                        "whose rss_state is in --parent-states may be assigned as "
                        "a transcript's parent; all others are donor-only. In a "
                        "gene-conversion locus a heavily converted transcript "
                        "best-matches a donor rather than its true parent, so "
                        "unconstrained argmax assignment is systematically wrong.")
    p.add_argument("--parent-states", nargs="+", default=["rss_present"],
                   help="rss_state values that qualify a gene as rearrangeable "
                        "(default: rss_present)")
    p.add_argument("--parent-expression-rescue", type=float, metavar="IDENTITY",
                   help="Also admit an RSS-negative gene as a candidate parent if "
                        "some transcript matches it at >= IDENTITY over the "
                        "assessable region. RSS screens are built for precision, "
                        "not recall, so a gene with a genuine but unannotated RSS "
                        "would otherwise be excluded forever. Use with care: a "
                        "transcript converted along its whole length also matches "
                        "a donor at high identity, so this can readmit donors as "
                        "parents. Compare germline-germline identities against "
                        "the transcript identity distribution before trusting it.")
    p.add_argument("--out-genes", required=True,
                   help="Per-gene functionality verdict TSV")
    p.add_argument("--out-assignments", required=True,
                   help="Per-transcript best-gene assignment TSV")
    args = p.parse_args()

    vgenes = {}
    for path in args.vgene_fastas:
        vgenes.update(read_fasta(path))
    info = parse_gene_names(vgenes)

    allowed_parents = None
    if args.candidate_parents:
        allowed_parents = set()
        with open(args.candidate_parents) as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            for line in fh:
                r = dict(zip(hdr, line.rstrip("\n").split("\t")))
                if r.get("rss_state") in args.parent_states:
                    allowed_parents.add(r["gene"])
        if not allowed_parents:
            raise SystemExit(
                f"No genes matched --parent-states {args.parent_states} in "
                f"{args.candidate_parents}; nothing could be a parent.")
        if args.parent_expression_rescue is not None:
            # Pre-scan for RSS-negative genes that some transcript matches at
            # very high identity; those may be real parents whose RSS the screen
            # missed. Scored on the same assessable region as everything else.
            rescued = set()
            for paf in args.pafs:
                for rec in parse_paf(paf):
                    if rec.target in allowed_parents or rec.target not in vgenes:
                        continue
                    sc = score_alignment(rec, vgenes[rec.target], args.junction_margin)
                    if (sc and sc["covered"] >= args.min_covered_bp
                            and sc["identity"] >= args.parent_expression_rescue):
                        rescued.add(rec.target)
            if rescued:
                print(f"Expression rescue (>={args.parent_expression_rescue}) added "
                      f"{len(rescued)} RSS-negative gene(s) as candidate parents: "
                      + ", ".join(sorted(rescued)), file=sys.stderr)
                allowed_parents |= rescued
            else:
                print(f"Expression rescue (>={args.parent_expression_rescue}) added "
                      "no genes", file=sys.stderr)
        print(f"Restricting parents to {len(allowed_parents)} rearrangeable gene(s): "
              + ", ".join(sorted(allowed_parents)), file=sys.stderr)

    # best (and runner-up) gene per transcript
    best = {}
    runner_up = {}
    for paf in args.pafs:
        for rec in parse_paf(paf):
            vseq = vgenes.get(rec.target)
            if vseq is None:
                continue
            if allowed_parents is not None and rec.target not in allowed_parents:
                continue
            sc = score_alignment(rec, vseq, args.junction_margin)
            if sc is None or sc["covered"] < args.min_covered_bp:
                continue
            cand = dict(sc, gene=rec.target, locus=info[rec.target].locus)
            cur = best.get(rec.query)
            if cur is None or (cand["identity"], cand["covered"]) > (cur["identity"], cur["covered"]):
                if cur is not None:
                    runner_up[rec.query] = cur
                best[rec.query] = cand
            else:
                ru = runner_up.get(rec.query)
                if ru is None or (cand["identity"], cand["covered"]) > (ru["identity"], ru["covered"]):
                    runner_up[rec.query] = cand

    # aggregate per gene
    per_gene = defaultdict(list)
    for tid, b in best.items():
        per_gene[b["gene"]].append(b)

    with open(args.out_assignments, "w") as fh:
        fh.write("transcript\tlocus\tbest_gene\tidentity\tmismatches\tcovered_bp"
                 "\trunner_up_gene\trunner_up_identity\tmargin\n")
        for tid in sorted(best):
            b = best[tid]
            ru = runner_up.get(tid)
            ru_gene = ru["gene"] if ru else "NA"
            ru_id = f"{ru['identity']:.4f}" if ru else "NA"
            margin = f"{b['identity'] - ru['identity']:.4f}" if ru else "NA"
            fh.write(f"{tid}\t{b['locus']}\t{b['gene']}\t{b['identity']:.4f}"
                     f"\t{b['mismatches']}\t{b['covered']}\t{ru_gene}\t{ru_id}\t{margin}\n")

    n_used = 0
    with open(args.out_genes, "w") as fh:
        fh.write("gene\tlocus\tpos\tstrand\tannotated_productive\tn_transcripts"
                 "\tbest_identity\tmedian_identity\tmin_mismatches\tused\n")
        for gene in sorted(vgenes, key=lambda g: (info[g].locus, info[g].pos)):
            hits = per_gene.get(gene, [])
            gi = info[gene]
            if hits:
                ids = [h["identity"] for h in hits]
                best_id = max(ids)
                med_id = statistics.median(ids)
                min_mm = min(h["mismatches"] for h in hits)
            else:
                best_id = med_id = 0.0
                min_mm = -1
            used = (len(hits) >= args.min_transcripts
                    and best_id >= args.min_identity)
            n_used += used
            fh.write(f"{gene}\t{gi.locus}\t{gi.pos}\t{gi.strand}\t{gi.productive}"
                     f"\t{len(hits)}\t{best_id:.4f}\t{med_id:.4f}\t{min_mm}\t{used}\n")

    # concordance with the annotation, printed so disagreement is visible
    agree = disagree_ann_only = disagree_expr_only = 0
    for gene in vgenes:
        hits = per_gene.get(gene, [])
        used = (len(hits) >= args.min_transcripts
                and (max((h["identity"] for h in hits), default=0.0) >= args.min_identity))
        ann = info[gene].productive
        if used and ann:
            agree += 1
        elif ann and not used:
            disagree_ann_only += 1
        elif used and not ann:
            disagree_expr_only += 1

    print(
        f"Transcripts assigned      : {len(best)}\n"
        f"Genes called used         : {n_used} / {len(vgenes)}\n"
        f"  annotated productive too: {agree}\n"
        f"  annotated pseudogene    : {disagree_expr_only}  <- expressed despite the annotation\n"
        f"  productive but unused   : {disagree_ann_only}\n"
        f"Thresholds: identity>={args.min_identity}, covered>={args.min_covered_bp}bp, "
        f"transcripts>={args.min_transcripts}, junction margin={args.junction_margin}bp",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
