"""
Per-gene alignment position maps: matches vs mismatches along each V gene.

For each V gene with sufficient alignment coverage, shows a pileup plot where:
  - Green fill  = bases that match the reference at that position
  - Orange fill = bases that mismatch the reference at that position
  - Gray fill   = total coverage (matches + mismatches)

CS string parsing:
  :N   → N consecutive matches
  *xy  → single mismatch (ref x → query y)
  +seq → insertion in query (skipped; does not advance target position)
  -seq → deletion in query (advances target position, counted as gap/mismatch)

Inputs are the *_detailed.paf files (produced with --cs flag) for all samples.
Only genes with >= min_alignments alignments are plotted.
Output is a multi-page PDF, sorted by total alignment count (most covered first).
"""
import argparse
import re
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np


# ─── CS string parser ─────────────────────────────────────────────────────────

_CS_RE = re.compile(r'(:[0-9]+|\*[a-z][a-z]|[+\-][a-z]+)')


def parse_cs(cs_tag):
    """Parse a CS string into list of (op, value) tuples.

    op values: 'M' (match count), 'X' (mismatch), 'I' (insert seq), 'D' (del seq)
    """
    ops = []
    for token in _CS_RE.findall(cs_tag):
        if token[0] == ':':
            ops.append(('M', int(token[1:])))
        elif token[0] == '*':
            ops.append(('X', 1))
        elif token[0] == '+':
            ops.append(('I', len(token) - 1))
        elif token[0] == '-':
            ops.append(('D', len(token) - 1))
    return ops


def accumulate_positions(target_start, target_len, cs_ops,
                         n_match, n_mismatch, n_covered):
    """Update per-position arrays given one alignment's CS ops."""
    pos = target_start
    for op, val in cs_ops:
        if op == 'M':
            end = min(pos + val, target_len)
            if end > pos:
                n_match[pos:end] += 1
                n_covered[pos:end] += 1
            pos += val
        elif op == 'X':
            if pos < target_len:
                n_mismatch[pos] += 1
                n_covered[pos] += 1
            pos += 1
        elif op == 'D':
            end = min(pos + val, target_len)
            if end > pos:
                n_mismatch[pos:end] += 1
                n_covered[pos:end] += 1
            pos += val
        # insertions don't advance target position
        if pos >= target_len:
            break


# ─── PAF parser ───────────────────────────────────────────────────────────────

def parse_paf(paf_path):
    with open(paf_path) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            f = line.rstrip().split('\t')
            cs_tag = None
            for field in f[12:]:
                if field.startswith('cs:Z:'):
                    cs_tag = field[5:]
                    break
            yield {
                'query':        f[0],
                'target':       f[5],
                'target_len':   int(f[6]),
                'target_start': int(f[7]),
                'target_end':   int(f[8]),
                'n_matches':    int(f[9]),
                'aln_block':    int(f[10]),
                'cs':           cs_tag,
            }


# ─── Collect pileups across all PAF files ─────────────────────────────────────

def build_pileups(paf_paths):
    """Return dict: gene → {'len': int, 'match': array, 'mismatch': array, 'covered': array, 'n_transcripts': int}.

    Only the best-identity alignment per transcript is used — matching how the
    combined summary plot assigns each transcript to exactly one gene.
    """
    # Pass 1: for each transcript, find its best-identity gene and record
    best = {}  # transcript → {'gene', 'identity', 'rec'}
    all_recs = []
    for paf_path in paf_paths:
        for rec in parse_paf(paf_path):
            all_recs.append(rec)
            tid = rec['query']
            aln_block = rec['aln_block']
            identity = rec['n_matches'] / aln_block if aln_block > 0 else 0.0
            rec['_identity'] = identity
            if tid not in best or identity > best[tid]['identity']:
                best[tid] = {'gene': rec['target'], 'identity': identity}

    # Pass 2: accumulate pileup arrays only for each transcript's best-hit gene
    gene_len      = {}
    gene_match    = {}
    gene_mismatch = {}
    gene_covered  = {}
    gene_transcripts = defaultdict(set)

    for rec in all_recs:
        tid  = rec['query']
        gene = rec['target']
        if best[tid]['gene'] != gene:
            continue  # skip secondary hits

        tlen = rec['target_len']
        if gene not in gene_len:
            gene_len[gene]      = tlen
            gene_match[gene]    = np.zeros(tlen, dtype=np.int32)
            gene_mismatch[gene] = np.zeros(tlen, dtype=np.int32)
            gene_covered[gene]  = np.zeros(tlen, dtype=np.int32)

        gene_transcripts[gene].add(tid)

        if rec['cs'] is None:
            s, e = rec['target_start'], min(rec['target_end'], tlen)
            gene_covered[gene][s:e] += 1
            continue

        cs_ops = parse_cs(rec['cs'])
        accumulate_positions(
            rec['target_start'], tlen, cs_ops,
            gene_match[gene], gene_mismatch[gene], gene_covered[gene]
        )

    result = {}
    for gene in gene_len:
        result[gene] = {
            'len':           gene_len[gene],
            'match':         gene_match[gene],
            'mismatch':      gene_mismatch[gene],
            'covered':       gene_covered[gene],
            'n_transcripts': len(gene_transcripts[gene]),
        }
    return result


# ─── Plotting ─────────────────────────────────────────────────────────────────

GENES_PER_PAGE = 6


def plot_gene(ax, gene, data):
    tlen  = data['len']
    x     = np.arange(tlen)
    match    = data['match'].astype(float)
    mismatch = data['mismatch'].astype(float)
    covered  = data['covered'].astype(float)

    # Stack: gray background (total covered), then match (green), mismatch stacked on top
    ax.fill_between(x, 0, covered,  step='mid', color='#cccccc', label='covered')
    ax.fill_between(x, 0, match,    step='mid', color='#2ca02c', alpha=0.85, label='match')
    ax.fill_between(x, match, match + mismatch, step='mid',
                    color='#d62728', alpha=0.85, label='mismatch')

    ymax = covered.max() if covered.max() > 0 else 1
    ax.set_xlim(0, tlen)
    ax.set_ylim(0, ymax * 1.1)
    ax.set_ylabel('# transcripts', fontsize=7)
    ax.set_title(f'{gene}  (n={data["n_transcripts"]} transcripts, {tlen} bp)', fontsize=8)
    ax.tick_params(axis='both', labelsize=7)


def make_legend(fig):
    patches = [
        mpatches.Patch(color='#cccccc', label='covered (no CS detail)'),
        mpatches.Patch(color='#2ca02c', label='match'),
        mpatches.Patch(color='#d62728', label='mismatch / deletion'),
    ]
    fig.legend(handles=patches, loc='lower center', ncol=3,
               fontsize=8, bbox_to_anchor=(0.5, 0.01))


def write_pdf(pileups, output_path, min_alignments):
    # Sort genes by alignment count descending
    genes = sorted(
        [g for g, d in pileups.items() if d['n_transcripts'] >= min_alignments],
        key=lambda g: pileups[g]['n_transcripts'],
        reverse=True
    )

    if not genes:
        # Write a blank page so Snakemake output file exists
        with PdfPages(output_path) as pdf:
            fig, ax = plt.subplots(figsize=(8, 4))
            ax.text(0.5, 0.5, f'No genes with >= {min_alignments} alignments',
                    ha='center', va='center', transform=ax.transAxes, fontsize=12)
            ax.axis('off')
            pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)
        return len(genes)

    with PdfPages(output_path) as pdf:
        for page_start in range(0, len(genes), GENES_PER_PAGE):
            page_genes = genes[page_start: page_start + GENES_PER_PAGE]
            n = len(page_genes)
            fig, axes = plt.subplots(n, 1, figsize=(14, n * 2.2 + 1.0),
                                     squeeze=False)
            for i, gene in enumerate(page_genes):
                plot_gene(axes[i, 0], gene, pileups[gene])
            axes[-1, 0].set_xlabel('Position in V gene (bp)', fontsize=8)
            make_legend(fig)
            fig.suptitle('V gene alignment position maps', fontsize=10, y=1.0)
            plt.tight_layout(rect=[0, 0.06, 1, 0.98])
            pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)

    return len(genes)


# ─── CLI ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pafs', nargs='+', required=True,
                        help='Detailed PAF files (must have been produced with --cs flag)')
    parser.add_argument('--output', required=True, help='Output PDF path')
    parser.add_argument('--min-alignments', type=int, default=2,
                        help='Minimum alignments for a gene to appear in the plot (default: 2)')
    args = parser.parse_args()

    print(f'Building pileups from {len(args.pafs)} PAF file(s)…', file=sys.stderr)
    pileups = build_pileups(args.pafs)
    print(f'  {len(pileups)} unique V genes found', file=sys.stderr)

    n_plotted = write_pdf(pileups, args.output, args.min_alignments)
    print(f'  {n_plotted} genes plotted (min_transcripts={args.min_alignments})', file=sys.stderr)
    print(f'Written: {args.output}', file=sys.stderr)


if __name__ == '__main__':
    main()
