"""
Convert BrepConvert's CSV output into the same tract schema gc_detect_tracts.py
emits, so both methods can be plotted and scored by identical code.

One row is written per (event, donor).  Where BrepConvert lists several donors
for an event they all appear, since that ambiguity is part of its result.
Events with gene == NA (a clustered-mismatch stretch it could not assign) are
dropped, as they carry no donor to test.

The parent is taken from BrepConvert's own 'allele' column so the method is
judged on its own parent call, falling back to our assignment when absent.
"""
import argparse
import csv
import sys


def read_tsv(path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(header, line.rstrip("\n").split("\t")))
                for line in fh if line.strip()]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--brepconvert-csv", required=True)
    ap.add_argument("--gene-map", required=True)
    ap.add_argument("--transcript-map", required=True)
    ap.add_argument("--assignments", required=True)
    ap.add_argument("--donor-pool", required=True)
    ap.add_argument("--locus", required=True)
    ap.add_argument("--resolve-donors", default="all",
                    choices=["all", "possible-first", "top1"],
                    help="How to handle events listing several donors. "
                         "'all' keeps every one (each is a separate claim). "
                         "'possible-first' keeps only the topologically possible "
                         "donors when at least one exists, which is what "
                         "BrepConvert's own filterResultsByGeneOrder() does -- "
                         "note this makes the impossible-donor rate circular, so "
                         "the result is for display, NOT evidence. "
                         "'top1' keeps only the first listed donor.")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    alias2gene = {r["alias"]: r["gene"] for r in read_tsv(args.gene_map)}
    alias2tx = {r["alias"]: r["transcript"] for r in read_tsv(args.transcript_map)}
    fallback_parent = {r["transcript"]: r["best_gene"]
                       for r in read_tsv(args.assignments) if r["locus"] == args.locus}

    allowed_of, mech_of = {}, {}
    for r in read_tsv(args.donor_pool):
        allowed_of[r["rearranged_gene"]] = (
            set() if r["allowed_donors"] == "NONE"
            else set(r["allowed_donors"].split(",")))
        mech_of[r["rearranged_gene"]] = r["mechanism"]

    n_in = n_out = n_na = n_unresolvable = 0
    with open(args.brepconvert_csv) as fh, open(args.out, "w") as out:
        out.write("transcript\tparent\tdonor\tmechanism\tdonor_allowed\tstart\tend"
                  "\tspan_bp\tn_support\tp_raw\tp_corrected\tsignificant\n")
        rdr = csv.DictReader(fh)
        if not rdr.fieldnames or "SeqID" not in rdr.fieldnames:
            print("empty or unrecognised BrepConvert output", file=sys.stderr)
            return
        for row in rdr:
            n_in += 1
            gene_field = (row.get("gene") or "").strip()
            if gene_field in ("", "NA"):
                n_na += 1
                continue
            tx = alias2tx.get(row["SeqID"], row["SeqID"])
            allele = (row.get("allele") or "").strip()
            parent = alias2gene.get(allele) or fallback_parent.get(tx)
            if parent is None:
                continue
            allowed = allowed_of.get(parent, set())
            mech = mech_of.get(parent, "unknown")
            try:
                start, end = int(row["start"]), int(row["end"])
            except (TypeError, ValueError):
                continue
            donors = [alias2gene.get(a.strip(), a.strip())
                      for a in gene_field.split(";") if a.strip()]
            donors = [d for d in donors if d != parent]
            if not donors:
                continue
            if args.resolve_donors == "top1":
                donors = donors[:1]
            elif args.resolve_donors == "possible-first":
                possible = [d for d in donors if d in allowed]
                if possible:
                    donors = possible          # ambiguity resolved in favour of biology
                else:
                    n_unresolvable += 1        # every option impossible -> a real problem
            for donor in donors:
                out.write(f"{tx}\t{parent}\t{donor}\t{mech}\t{donor in allowed}"
                          f"\t{start}\t{end}\t{end - start + 1}\tNA\tNA\tNA\tTrue\n")
                n_out += 1

    print(f"{args.locus}: {n_in} BrepConvert rows -> {n_out} (transcript,donor) tracts "
          f"[resolve-donors={args.resolve_donors}] "
          f"({n_na} rows had no donor assigned)", file=sys.stderr)
    if args.resolve_donors == "possible-first":
        print(f"  events where EVERY listed donor was impossible: {n_unresolvable} "
              "-- these are the calls the topology constraint genuinely rejects",
              file=sys.stderr)


if __name__ == "__main__":
    main()
