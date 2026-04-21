"""
For each IG transcript, collect the top-N V gene alignments ranked by identity,
then compute detailed metrics for each.  This quantifies gene conversion
signatures: if a transcript is a mosaic of multiple V genes, you'll see the
top alignments cover overlapping or adjacent transcript regions at high
(but non-100%) identity to different donor genes.

Metrics per alignment rank:
  - gene name + locus
  - identity     : n_matches / aln_block_len
  - v_coverage   : fraction of V gene length covered (aln_block / target_len)
  - query_frac   : fraction of transcript covered  (aln_block / query_len)
  - query_start/end on transcript
  - delta_to_rank1_identity : how much identity drops from best hit
"""
import argparse
import sys
from collections import defaultdict
from pathlib import Path


def parse_paf(paf_path):
    locus = Path(paf_path).stem.replace("_detailed", "")
    with open(paf_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip().split("\t")
            n_matches = int(f[9])
            aln_block = int(f[10])
            target_len = int(f[6])
            query_len = int(f[1])
            yield {
                "query": f[0],
                "query_len": query_len,
                "query_start": int(f[2]),
                "query_end": int(f[3]),
                "strand": f[4],
                "gene": f[5],
                "target_len": target_len,
                "n_matches": n_matches,
                "aln_block": aln_block,
                "mapq": int(f[11]),
                "locus": locus,
                "identity": n_matches / aln_block if aln_block > 0 else 0.0,
                "v_coverage": aln_block / target_len if target_len > 0 else 0.0,
                "query_frac": aln_block / query_len if query_len > 0 else 0.0,
                "is_exact": (aln_block == target_len and n_matches == target_len),
            }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pafs", nargs="+", required=True)
    parser.add_argument("--top-n", type=int, default=5)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    # Collect all alignments per transcript
    transcript_alns = defaultdict(list)

    for paf_path in args.pafs:
        for rec in parse_paf(paf_path):
            transcript_alns[rec["query"]].append(rec)

    with open(args.output, "w") as fh:
        fh.write(
            "transcript\ttranscript_len\t"
            "rank\t"
            "gene\tlocus\t"
            "identity\tv_coverage\tquery_frac\t"
            "query_start\tquery_end\tstrand\t"
            "n_matches\taln_block\ttarget_len\t"
            "delta_identity_from_rank1\t"
            "is_exact\t"
            "mapq\n"
        )

        for tid in sorted(transcript_alns):
            alns = transcript_alns[tid]
            # Sort by identity descending, then by v_coverage descending as tie-break
            alns.sort(key=lambda r: (r["identity"], r["v_coverage"]), reverse=True)
            top = alns[: args.top_n]
            best_identity = top[0]["identity"] if top else 0.0

            for rank, rec in enumerate(top, start=1):
                delta = best_identity - rec["identity"]
                fh.write(
                    f"{rec['query']}\t{rec['query_len']}\t"
                    f"{rank}\t"
                    f"{rec['gene']}\t{rec['locus']}\t"
                    f"{rec['identity']:.6f}\t{rec['v_coverage']:.6f}\t{rec['query_frac']:.6f}\t"
                    f"{rec['query_start']}\t{rec['query_end']}\t{rec['strand']}\t"
                    f"{rec['n_matches']}\t{rec['aln_block']}\t{rec['target_len']}\t"
                    f"{delta:.6f}\t"
                    f"{rec['is_exact']}\t"
                    f"{rec['mapq']}\n"
                )

    n_transcripts = len(transcript_alns)
    print(f"Top-{args.top_n} alignments written for {n_transcripts} transcripts", file=sys.stderr)


if __name__ == "__main__":
    main()
