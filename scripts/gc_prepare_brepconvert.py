"""
Prepare inputs for BrepConvert.

Two things need fixing before BrepConvert will accept our data.

1. Gene names.  batchConvertAnalysis() rewrites germline names with

       stringr::str_extract(names, "IG.V[0-9A-Z\\-\\*]*\\||IG.V.*$")

   which assumes IMGT naming (IGHV1-2*01|...).  Our names look like
   VGP_redwinged_blackbird_IGL.8749176.CM036732.1.V.True.- , where the
   characters after "IG" are "L." rather than "?V", so the regex returns NA
   and Biostrings dies with "'names(x)' contains NAs".  We therefore emit
   IMGT-shaped placeholders (IGLV1, IGLV2, ...) plus a mapping file so calls
   can be translated back to real gene names afterwards.

2. Transcript names.  Ids containing "/" break the temporary FASTA/BLAT
   round-trip, so they are sanitised too.

The functional/pseudogene split comes from the expression-based call in
gc_call_functional_genes.py, not from the annotation's Productive flag.
"""
import argparse
import sys

from gc_lib import read_fasta, write_fasta


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--vgene-fasta", required=True)
    p.add_argument("--functional-genes", required=True,
                   help="functional_genes.tsv from gc_call_functional_genes.py")
    p.add_argument("--assignments", required=True)
    p.add_argument("--transcripts", required=True,
                   help="combined_ig_transcripts.fasta")
    p.add_argument("--locus", required=True)
    p.add_argument("--out-functional", required=True)
    p.add_argument("--out-pseudogene", required=True)
    p.add_argument("--out-transcripts", required=True)
    p.add_argument("--out-gene-map", required=True)
    p.add_argument("--out-transcript-map", required=True)
    args = p.parse_args()

    genes = read_fasta(args.vgene_fasta)

    used = set()
    with open(args.functional_genes) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            r = dict(zip(header, line.rstrip("\n").split("\t")))
            if r["locus"] == args.locus and r["used"] == "True":
                used.add(r["gene"])
    if not used:
        raise SystemExit(
            f"No expression-confirmed functional {args.locus} genes. BrepConvert "
            "needs at least one functional allele to align against.")

    # IMGT-shaped placeholder names, ordered so numbering is stable
    gene_map, functional, pseudo = {}, {}, {}
    for i, name in enumerate(sorted(genes), start=1):
        alias = f"{args.locus}V{i}"
        gene_map[alias] = name
        (functional if name in used else pseudo)[alias] = genes[name]

    if not pseudo:
        raise SystemExit(f"No donor (pseudogene) sequences left for {args.locus}")

    write_fasta(args.out_functional, functional)
    write_fasta(args.out_pseudogene, pseudo)

    ids = set()
    with open(args.assignments) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            r = dict(zip(header, line.rstrip("\n").split("\t")))
            if r["locus"] == args.locus:
                ids.add(r["transcript"])

    tx = read_fasta(args.transcripts)
    tx_map, out_tx = {}, {}
    for name in sorted(ids):
        if name not in tx:
            continue
        alias = name.replace("/", "_")
        tx_map[alias] = name
        out_tx[alias] = tx[name]
    write_fasta(args.out_transcripts, out_tx)

    with open(args.out_gene_map, "w") as fh:
        fh.write("alias\tgene\trole\n")
        for alias, name in gene_map.items():
            fh.write(f"{alias}\t{name}\t{'functional' if name in used else 'pseudogene'}\n")
    with open(args.out_transcript_map, "w") as fh:
        fh.write("alias\ttranscript\n")
        for alias, name in tx_map.items():
            fh.write(f"{alias}\t{name}\n")

    print(f"{args.locus}: {len(functional)} functional, {len(pseudo)} donor genes, "
          f"{len(out_tx)} transcripts", file=sys.stderr)


if __name__ == "__main__":
    main()
