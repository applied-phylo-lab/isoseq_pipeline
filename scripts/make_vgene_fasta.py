"""
Convert a `combined_genes_<LOCUS>_clean.txt` germline annotation table into the
V gene FASTA format the pipeline expects.

Input is a tab-separated table with the columns

    GeneType  Contig  Pos  Strand  Sequence  Productive  Locus

Only rows whose GeneType matches --gene-type (default "V") are written out.
Headers follow the convention already used by the tufted_duck databases:

    >{prefix}_{locus}.{pos}.{contig}.{gene_type}.{productive}.{strand}

Identical sequences are collapsed to a single record (keeping the first
occurrence) because duplicate reference entries would split minimap2's top-N
alignments across indistinguishable targets.
"""
import argparse
import sys
from collections import Counter


REQUIRED_COLUMNS = ["GeneType", "Contig", "Pos", "Strand",
                    "Sequence", "Productive", "Locus"]


def read_table(path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        missing = [c for c in REQUIRED_COLUMNS if c not in header]
        if missing:
            raise SystemExit(f"{path}: missing required column(s): {', '.join(missing)}")
        idx = {c: header.index(c) for c in REQUIRED_COLUMNS}
        for lineno, line in enumerate(fh, start=2):
            line = line.rstrip("\n")
            if not line.strip():
                continue
            f = line.split("\t")
            if len(f) < len(header):
                raise SystemExit(f"{path}:{lineno}: expected {len(header)} fields, got {len(f)}")
            yield {c: f[i] for c, i in idx.items()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--table", required=True,
                        help="combined_genes_<LOCUS>_clean.txt input table")
    parser.add_argument("--prefix", required=True,
                        help="Species/assembly prefix for FASTA headers, "
                             "e.g. VGP_redwinged_blackbird")
    parser.add_argument("--locus", required=True,
                        help="Locus name used in the header, e.g. IGH")
    parser.add_argument("--gene-type", default="V",
                        help="GeneType rows to keep (default: V)")
    parser.add_argument("--min-length", type=int, default=0,
                        help="Drop sequences shorter than this many bp (default: 0)")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    rows = list(read_table(args.table))
    kept, seen_seqs, seen_ids = [], {}, Counter()
    n_wrong_type = n_dup_seq = n_too_short = n_wrong_locus = 0

    for row in rows:
        if row["GeneType"] != args.gene_type:
            n_wrong_type += 1
            continue
        if row["Locus"] != args.locus:
            n_wrong_locus += 1
            continue
        seq = row["Sequence"].strip().upper()
        if len(seq) < args.min_length:
            n_too_short += 1
            continue
        if seq in seen_seqs:
            n_dup_seq += 1
            continue

        name = (f"{args.prefix}_{args.locus}.{row['Pos']}.{row['Contig']}."
                f"{row['GeneType']}.{row['Productive']}.{row['Strand']}")
        seen_ids[name] += 1
        if seen_ids[name] > 1:
            # Same locus position annotated twice — keep both but keep IDs unique.
            name = f"{name}_{seen_ids[name]}"
        seen_seqs[seq] = name
        kept.append((name, seq))

    with open(args.output, "w") as fh:
        for name, seq in kept:
            fh.write(f">{name}\n{seq}\n")

    print(
        f"{args.table} -> {args.output}\n"
        f"  rows read            : {len(rows)}\n"
        f"  skipped (gene type)  : {n_wrong_type}\n"
        f"  skipped (locus)      : {n_wrong_locus}\n"
        f"  skipped (too short)  : {n_too_short}\n"
        f"  skipped (duplicate)  : {n_dup_seq}\n"
        f"  written              : {len(kept)}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
