"""
Detect germline V gene sequences that appear EXACTLY (100% identity, 100%
gene coverage) within IG transcripts, indicating a functionally rearranged
(VDJ-recombined) locus.

A match is "exact" when:
  - alignment block length == target (V gene) length  → no indels
  - n_matches == target length                        → no mismatches
  i.e. the complete V gene sequence is present verbatim inside the transcript.

Outputs
-------
exact_match_summary.tsv   — per-gene counts across all transcripts
exact_match_per_transcript.tsv — which transcripts contain which exact genes
"""
import argparse
import sys
from collections import defaultdict
from pathlib import Path


def parse_paf(paf_path):
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
                "locus": Path(paf_path).stem.replace("_detailed", ""),
            }


def is_exact(rec):
    tlen = rec["target_len"]
    return rec["aln_block"] == tlen and rec["n_matches"] == tlen


def read_vgene_lengths(fasta_paths):
    """Return dict gene_name → length for all V gene FASTAs."""
    lengths = {}
    for path in fasta_paths:
        current, seq_len = None, 0
        with open(path) as fh:
            for line in fh:
                line = line.rstrip()
                if line.startswith(">"):
                    if current:
                        lengths[current] = seq_len
                    current = line[1:].split()[0]
                    seq_len = 0
                else:
                    seq_len += len(line)
        if current:
            lengths[current] = seq_len
    return lengths


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pafs", nargs="+", required=True)
    parser.add_argument("--vgene-fastas", nargs="+", required=True)
    parser.add_argument("--out-summary", required=True)
    parser.add_argument("--out-per-transcript", required=True)
    args = parser.parse_args()

    gene_lengths = read_vgene_lengths(args.vgene_fastas)

    # gene → set of transcripts that contain it exactly
    gene_exact_hits = defaultdict(set)
    # transcript → set of genes found exactly within it
    transcript_exact_hits = defaultdict(set)
    # all transcripts seen
    all_transcripts = set()
    # all genes per locus
    gene_to_locus = {}

    for paf_path in args.pafs:
        locus = Path(paf_path).stem.replace("_detailed", "")
        for rec in parse_paf(paf_path):
            all_transcripts.add(rec["query"])
            gene_to_locus[rec["target"]] = locus
            if is_exact(rec):
                gene_exact_hits[rec["target"]].add(rec["query"])
                transcript_exact_hits[rec["query"]].add(rec["target"])

    n_transcripts = len(all_transcripts)

    # ── Summary: per gene ──────────────────────────────────────────────────────
    with open(args.out_summary, "w") as fh:
        fh.write(
            "gene\tlocus\tgene_length_bp\t"
            "transcripts_with_exact_match\t"
            "fraction_of_ig_transcripts\t"
            "is_functional_candidate\n"
        )
        for gene in sorted(gene_to_locus):
            locus = gene_to_locus[gene]
            n_hits = len(gene_exact_hits.get(gene, set()))
            frac = n_hits / n_transcripts if n_transcripts > 0 else 0.0
            is_functional = n_hits > 0
            fh.write(
                f"{gene}\t{locus}\t{gene_lengths.get(gene, 'NA')}\t"
                f"{n_hits}\t{frac:.6f}\t{is_functional}\n"
            )

    # ── Per transcript: which genes are exact ──────────────────────────────────
    with open(args.out_per_transcript, "w") as fh:
        fh.write("transcript\tn_exact_genes\texact_genes\tloci\n")
        for tid in sorted(all_transcripts):
            genes = sorted(transcript_exact_hits.get(tid, set()))
            loci = sorted({gene_to_locus.get(g, "?") for g in genes})
            fh.write(
                f"{tid}\t{len(genes)}\t"
                f"{','.join(genes) if genes else 'none'}\t"
                f"{','.join(loci) if loci else 'none'}\n"
            )

    n_functional_genes = sum(1 for v in gene_exact_hits.values() if v)
    n_transcripts_with_exact = sum(1 for v in transcript_exact_hits.values() if v)

    print(
        f"Exact matches: {n_functional_genes} genes found verbatim in transcripts\n"
        f"  {n_transcripts_with_exact}/{n_transcripts} IG transcripts contain at least one exact V gene",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
