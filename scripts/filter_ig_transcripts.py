"""
Filter clustered IsoSeq transcripts to retain only those with significant
alignment to any germline V gene locus (IGH, IGL, IGK, ...).

PAF columns (0-indexed):
 0  query_name      6  target_name
 1  query_len       7  target_len
 2  query_start     8  target_start
 3  query_end       9  target_end
 4  strand         10  n_matches
 5  target_name    11  aln_block_len
               ...12  mapq
"""
import argparse
import sys
from collections import defaultdict
from pathlib import Path


def parse_paf(paf_path):
    """Yield dicts for each PAF alignment line."""
    with open(paf_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip().split("\t")
            yield {
                "query": f[0],
                "query_len": int(f[1]),
                "target": f[5],
                "target_len": int(f[6]),
                "n_matches": int(f[9]),
                "aln_block": int(f[10]),
                "mapq": int(f[11]),
            }


def identity(rec):
    return rec["n_matches"] / rec["aln_block"] if rec["aln_block"] > 0 else 0.0


def coverage(rec):
    """Fraction of V gene (target) covered by this alignment."""
    return rec["aln_block"] / rec["target_len"] if rec["target_len"] > 0 else 0.0


def read_fasta_ids(fasta_path):
    ids = []
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                ids.append(line[1:].split()[0])
    return ids


def read_fasta(fasta_path):
    seqs = {}
    current = None
    with open(fasta_path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                current = line[1:].split()[0]
                seqs[current] = []
            elif current:
                seqs[current].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--transcripts", required=True)
    parser.add_argument("--pafs", nargs="+", required=True)
    parser.add_argument("--min-identity", type=float, default=0.85)
    parser.add_argument("--min-coverage", type=float, default=0.50)
    parser.add_argument("--out-fasta", required=True)
    parser.add_argument("--out-ids", required=True)
    parser.add_argument("--out-stats", required=True)
    args = parser.parse_args()

    # Collect best alignment per transcript across all loci
    best = defaultdict(lambda: {"identity": 0.0, "coverage": 0.0, "locus": None, "gene": None})

    for paf_path in args.pafs:
        locus = Path(paf_path).stem  # e.g. "IGH"
        for rec in parse_paf(paf_path):
            tid = rec["query"]
            idt = identity(rec)
            cov = coverage(rec)
            if idt > best[tid]["identity"]:
                best[tid] = {"identity": idt, "coverage": cov, "locus": locus, "gene": rec["target"]}

    # Apply thresholds
    passing = {
        tid for tid, b in best.items()
        if b["identity"] >= args.min_identity and b["coverage"] >= args.min_coverage
    }

    # Read all transcript sequences
    seqs = read_fasta(args.transcripts)
    total = len(seqs)

    # Write filtered FASTA
    with open(args.out_fasta, "w") as fh:
        for tid in sorted(passing):
            if tid in seqs:
                fh.write(f">{tid}\n{seqs[tid]}\n")

    # Write ID list
    with open(args.out_ids, "w") as fh:
        for tid in sorted(passing):
            fh.write(f"{tid}\n")

    # Write per-locus stats
    locus_counts = defaultdict(int)
    for tid in passing:
        locus_counts[best[tid]["locus"]] += 1

    with open(args.out_stats, "w") as fh:
        fh.write("metric\tvalue\n")
        fh.write(f"total_transcripts\t{total}\n")
        fh.write(f"ig_transcripts\t{len(passing)}\n")
        fh.write(f"filtered_out\t{total - len(passing)}\n")
        fh.write(f"retention_rate\t{len(passing)/total:.4f}\n")
        for locus, cnt in sorted(locus_counts.items()):
            fh.write(f"best_locus_{locus}\t{cnt}\n")

    print(
        f"Retained {len(passing)}/{total} transcripts "
        f"({len(passing)/total*100:.1f}%) as IG-related",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
