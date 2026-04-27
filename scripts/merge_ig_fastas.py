"""
Merge two IG-filtered FASTA files (minimap2 filter + immunotools diversity_analyzer)
into a single deduplicated FASTA.

Deduplication is by sequence ID. If the same ID appears in both sources the
sequence from source A (minimap2) is kept and the duplicate is counted.

Outputs a stats TSV reporting how many reads came from each source exclusively
and how many overlapped, so you can see the added value of each tool.
"""
import argparse
import sys


def read_fasta(path):
    """Return ordered dict: id -> sequence string."""
    seqs = {}
    current = None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                current = line[1:].split()[0]
                seqs[current] = []
            elif current is not None:
                seqs[current].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--minimap-fasta", required=True,
                        help="ig_transcripts.fasta from the minimap2 filter step")
    parser.add_argument("--immunotools-fasta", required=True,
                        help="cleaned_sequences.fasta from diversity_analyzer")
    parser.add_argument("--out-fasta", required=True,
                        help="Combined deduplicated output FASTA")
    parser.add_argument("--out-stats", required=True,
                        help="TSV with per-source counts")
    args = parser.parse_args()

    minimap_seqs = read_fasta(args.minimap_fasta)
    immunotools_seqs = read_fasta(args.immunotools_fasta)

    minimap_ids = set(minimap_seqs)
    immunotools_ids = set(immunotools_seqs)

    both_ids       = minimap_ids & immunotools_ids
    minimap_only   = minimap_ids - immunotools_ids
    immuno_only    = immunotools_ids - minimap_ids

    # Build merged dict: minimap first (preserves its sequences for overlaps),
    # then add immunotools-only entries
    merged = dict(minimap_seqs)
    for tid in immuno_only:
        merged[tid] = immunotools_seqs[tid]

    # Write combined FASTA in a stable order: minimap first, then immuno-only
    with open(args.out_fasta, "w") as fh:
        for tid in list(minimap_seqs) + [t for t in immunotools_seqs if t in immuno_only]:
            fh.write(f">{tid}\n{merged[tid]}\n")

    # Write stats
    with open(args.out_stats, "w") as fh:
        fh.write("metric\tcount\n")
        fh.write(f"minimap2_filter_total\t{len(minimap_ids)}\n")
        fh.write(f"immunotools_total\t{len(immunotools_ids)}\n")
        fh.write(f"in_both\t{len(both_ids)}\n")
        fh.write(f"minimap2_only\t{len(minimap_only)}\n")
        fh.write(f"immunotools_only\t{len(immuno_only)}\n")
        fh.write(f"combined_total\t{len(merged)}\n")

    print(
        f"minimap2 filter : {len(minimap_ids):>6} transcripts\n"
        f"immunotools     : {len(immunotools_ids):>6} transcripts\n"
        f"overlap         : {len(both_ids):>6} transcripts\n"
        f"combined total  : {len(merged):>6} transcripts "
        f"({len(immuno_only)} added by immunotools)",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
