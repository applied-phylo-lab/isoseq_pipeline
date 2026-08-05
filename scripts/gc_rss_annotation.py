"""
Derive a per-V-gene functionality call from the recombination signal sequence
(RSS) rather than from IgDetective's Productive prediction.

Why RSS is the better criterion in principle
--------------------------------------------
A V gene can only be rearranged if RAG can bind the RSS immediately 3' of its
coding end.  A gene with an intact heptamer/nonamer is usable regardless of
what an ORF-based predictor thinks; a gene without one cannot be rearranged
even if its reading frame is perfect.  So RSS integrity, not Productive, is the
mechanistically meaningful test.

Few RSS-positive genes is the expected architecture, not a failure
------------------------------------------------------------------
In species that diversify by gene conversion -- chicken, duck, and by extension
these songbirds -- only a handful of V genes need to be rearrangeable.  The
chicken IGL locus has a single functional V gene backed by ~25 pseudogene
donors, and IGH is the same shape.  So finding an RSS on ~8% of V genes is what
the biology predicts: a small rearrangeable set plus a large donor array whose
members never need an RSS because they are only ever copied FROM.

This has a sharp consequence for parent assignment.  Only RSS-positive genes
can be the rearranged parent; every other gene is a donor.  Assigning a
transcript to whichever germline gene it most resembles is therefore wrong
whenever conversion has been extensive -- a heavily converted transcript will
best-match a donor rather than its true parent.  Pass this table to
gc_call_functional_genes.py --candidate-parents to constrain the assignment.

The call is binary, deliberately
--------------------------------
    rss_present  an RSS is recorded for this gene -> it can be rearranged
    rss_absent   none recorded -> donor-only

The upstream gene lists are already screened, so a recorded RSS is taken at
face value.  This script does NOT re-score the motifs or apply its own
mismatch budget: doing so would silently overrule that screening with an
arbitrary threshold.  Hamming distances to the canonical motifs are still
written out as informational columns, but they do not affect the call.

Canonical motifs used only for those reference columns:
heptamer CACAGTG, nonamer ACAAAAACC.
"""
import argparse
import csv
import sys
from collections import Counter

CANONICAL_HEPTAMER = "CACAGTG"
CANONICAL_NONAMER = "ACAAAAACC"


def hamming(a, b):
    """Mismatches over the compared length; length difference counts as mismatch."""
    a, b = a.upper(), b.upper()
    d = sum(1 for x, y in zip(a, b) if x != y)
    return d + abs(len(a) - len(b))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gene-list-csv", required=True,
                    help="clean_birds/gene_list.csv (IGH+all loci) or the IGL variant")
    ap.add_argument("--source-pattern", required=True,
                    help="Substring identifying the species/assembly in the Source "
                         "column, e.g. RedWinged_Blackbird")
    ap.add_argument("--vgene-fastas", nargs="+", required=True,
                    help="Pipeline V gene FASTAs, used to join on (locus,pos,strand)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    sys.path.insert(0, "scripts")
    from gc_lib import read_fasta, parse_gene_names

    genes = {}
    for path in args.vgene_fastas:
        genes.update(read_fasta(path))
    info = parse_gene_names(genes)
    # join key: (locus, pos, strand)
    by_key = {(g.locus, g.pos, g.strand): name for name, g in info.items()}

    rss_by_key = {}
    with open(args.gene_list_csv, newline="") as fh:
        for r in csv.DictReader(fh):
            if args.source_pattern not in r["Source"]:
                continue
            if r["GeneType"] != "V":
                continue
            try:
                key = (r["Locus"], int(r["Pos"]), r["Strand"])
            except (KeyError, ValueError):
                continue
            rss_by_key[key] = r

    counts = Counter()
    with open(args.out, "w") as out:
        # rss_valid mirrors rss_state; kept as its own column for callers that
        # want a plain boolean rather than the three-state string.
        out.write("gene\tlocus\tpos\tstrand\tannotated_productive\theptamer\tnonamer"
                  "\tspacer_bp\theptamer_mismatch\tnonamer_mismatch\trss_state"
                  "\trss_valid\n")
        for key, name in sorted(by_key.items(), key=lambda kv: (kv[0][0], kv[0][1])):
            g = info[name]
            r = rss_by_key.get(key)
            hep = (r or {}).get("Heptamer", "").strip()
            non = (r or {}).get("Nonamer", "").strip()
            spacer = (r or {}).get("Number of bp Downstream (RSS)", "").strip() or "NA"
            if not hep and not non:
                state, valid, hm, nm = "rss_absent", "False", "NA", "NA"
            else:
                # Taken at face value: the gene lists are pre-screened, so an
                # RSS being recorded IS the criterion. The mismatch counts below
                # are reported for reference only and gate nothing.
                state, valid = "rss_present", "True"
                hm = str(hamming(hep, CANONICAL_HEPTAMER)) if hep else "NA"
                nm = str(hamming(non, CANONICAL_NONAMER)) if non else "NA"
            counts[state] += 1
            out.write(f"{name}\t{g.locus}\t{g.pos}\t{g.strand}\t{g.productive}"
                      f"\t{hep or 'NA'}\t{non or 'NA'}\t{spacer}\t{hm}\t{nm}"
                      f"\t{state}\t{valid}\n")

    n = sum(counts.values())
    print(f"{args.source_pattern}: {n} V genes joined to the gene list", file=sys.stderr)
    for state in ("rss_present", "rss_absent"):
        c = counts[state]
        print(f"  {state:11s}: {c:4d} ({c/n*100:5.1f}%)" if n else "", file=sys.stderr)
    if n and counts["rss_present"]:
        print(f"  -> {counts['rss_present']} gene(s) can be a rearranged parent; "
              f"the remaining {n - counts['rss_present']} are donor candidates. "
              "This is the expected gene-conversion architecture.", file=sys.stderr)


if __name__ == "__main__":
    main()
