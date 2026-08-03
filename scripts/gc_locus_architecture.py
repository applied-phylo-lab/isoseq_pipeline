"""
Test whether the arrangement of V genes along the locus is what cis-acting,
deletional gene conversion predicts.

The hypothesis
--------------
Under deletional recombination every gene between the rearranged V and J is
excised, so a V gene's position fixes two independent things:

  * how many donors it gets *as a recipient*  -> genes distal to it survive.
    A J-proximal gene keeps the whole array; the most distal gene keeps none.
  * how many recipients it can serve *as a donor* -> genes proximal to it.
    The most distal gene can donate to everything; the J-proximal gene to
    nothing.

So selection should push functional genes toward the J-proximal end (maximum
donor supply) and let pseudogenes accumulate at the distal end, where a gene
is still a useful donor but would be a poor recipient.  A pseudogene sitting
at the J-proximal end is evolutionarily inert -- it can neither be usefully
rearranged nor serve as a donor -- and should be depleted.

Tests performed
---------------
  1. Do functional genes have more allowed donors than pseudogenes?
     (Mann-Whitney U, one-sided.)  Run for the annotation's own
     Productive flag and, independently, for expression-confirmed usage.
  2. Does expression (transcripts per gene) increase with donor supply?
     (Spearman.)
  3. Does identity-to-germline *fall* as donor supply rises?  This is the
     conversion-load prediction: a gene with many available donors is
     converted more often, so its transcripts resemble its own germline
     less.  It is also the explanation for why the most-used gene need not
     be the one showing a 100% match.  (Spearman, per transcript.)

Note the direction of bias in test 3: calling a gene "used" requires high
identity, which works *against* finding a negative correlation.  A negative
result therefore understates the effect.
"""
import argparse
import sys

from scipy import stats


def read_tsv(path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(header, line.rstrip("\n").split("\t"))) for line in fh if line.strip()]


def mwu(a, b, alternative, label, fh):
    if len(a) < 3 or len(b) < 3:
        fh.write(f"{label}\tinsufficient_n\t{len(a)}\t{len(b)}\tNA\tNA\n")
        print(f"  {label}: too few genes (n={len(a)} vs {len(b)})", file=sys.stderr)
        return
    u, p = stats.mannwhitneyu(a, b, alternative=alternative)
    ma, mb = stats.tmean(a), stats.tmean(b)
    fh.write(f"{label}\tmannwhitney_{alternative}\t{len(a)}\t{len(b)}\t{u:.1f}\t{p:.4g}\n")
    print(f"  {label}: mean {ma:.1f} vs {mb:.1f}, U={u:.1f}, p={p:.4g}"
          f" (n={len(a)} vs {len(b)})", file=sys.stderr)


def spearman(x, y, label, fh):
    if len(x) < 4:
        fh.write(f"{label}\tspearman\t{len(x)}\tNA\tNA\tNA\n")
        return
    r, p = stats.spearmanr(x, y)
    fh.write(f"{label}\tspearman\t{len(x)}\tNA\t{r:.4f}\t{p:.4g}\n")
    print(f"  {label}: rho={r:.3f}, p={p:.4g} (n={len(x)})", file=sys.stderr)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--functional-genes", required=True)
    p.add_argument("--donor-pool", required=True)
    p.add_argument("--assignments", required=True)
    p.add_argument("--locus", required=True)
    p.add_argument("--out", required=True)
    args = p.parse_args()

    genes = {g["gene"]: g for g in read_tsv(args.functional_genes)
             if g["locus"] == args.locus}
    pool = {d["rearranged_gene"]: d for d in read_tsv(args.donor_pool)}
    assigns = [a for a in read_tsv(args.assignments) if a["locus"] == args.locus]

    # Restrict to deletional genes: an inversional gene retains the whole array
    # by a different mechanism and does not carry positional information.
    rows = []
    for name, g in genes.items():
        d = pool.get(name)
        if d is None or d["mechanism"] != "deletion":
            continue
        rows.append({
            "gene": name,
            "allowed": int(d["n_allowed_donors"]),
            "recipients": int(d["n_deleted"]),   # genes it could donate to
            "annot_productive": g["annotated_productive"] == "True",
            "used": g["used"] == "True",
            "n_transcripts": int(g["n_transcripts"]),
            "best_identity": float(g["best_identity"]),
        })
    if not rows:
        raise SystemExit(f"no deletional {args.locus} genes found")

    print(f"\n{args.locus}: {len(rows)} deletionally-rearranging V genes", file=sys.stderr)

    with open(args.out, "w") as fh:
        fh.write("test\tmethod\tn1\tn2\tstatistic\tp_value\n")

        print(" [1] donor supply vs functionality", file=sys.stderr)
        mwu([r["allowed"] for r in rows if r["annot_productive"]],
            [r["allowed"] for r in rows if not r["annot_productive"]],
            "greater", "annotated_functional_have_more_donors", fh)
        mwu([r["allowed"] for r in rows if r["used"]],
            [r["allowed"] for r in rows if not r["used"]],
            "greater", "expressed_genes_have_more_donors", fh)

        # NOT independent evidence: for a deletional gene, allowed donors and
        # potential recipients sum to a constant, so this is test 1 mirrored and
        # returns the identical U and p. Reported only because the framing
        # ("pseudogenes accumulate distally where they are still good donors")
        # is the useful one; do not count it as a second result.
        print(" [2] pseudogenes as donors -- mirror of [1], not independent",
              file=sys.stderr)
        mwu([r["recipients"] for r in rows if not r["annot_productive"]],
            [r["recipients"] for r in rows if r["annot_productive"]],
            "greater", "pseudogenes_can_serve_more_recipients_MIRROR_OF_TEST1", fh)

        print(" [3] expression vs donor supply", file=sys.stderr)
        spearman([r["allowed"] for r in rows], [r["n_transcripts"] for r in rows],
                 "transcripts_vs_allowed_donors", fh)

        print(" [4] conversion load: identity falls as donor supply rises", file=sys.stderr)
        # per-transcript, using each transcript's own identity to its assigned gene
        allowed_by_gene = {r["gene"]: r["allowed"] for r in rows}
        x, y = [], []
        for a in assigns:
            g = a["best_gene"]
            if g in allowed_by_gene:
                x.append(allowed_by_gene[g])
                y.append(float(a["identity"]))
        spearman(x, y, "transcript_identity_vs_allowed_donors", fh)

    # per-gene table for plotting / inspection
    table = args.out.replace(".tsv", "_per_gene.tsv")
    with open(table, "w") as fh:
        fh.write("gene\tallowed_donors\tpotential_recipients\tannot_productive"
                 "\tused\tn_transcripts\tbest_identity\n")
        for r in sorted(rows, key=lambda r: -r["allowed"]):
            fh.write(f"{r['gene']}\t{r['allowed']}\t{r['recipients']}"
                     f"\t{r['annot_productive']}\t{r['used']}"
                     f"\t{r['n_transcripts']}\t{r['best_identity']:.4f}\n")
    print(f"\nwrote {args.out} and {table}", file=sys.stderr)


if __name__ == "__main__":
    main()
