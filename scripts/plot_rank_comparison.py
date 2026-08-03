"""
Rank-1 vs rank-2 alignment comparison for the top V genes.

Produces two page types per (rank-1 gene, rank-2 partner) pair:

PAGE TYPE 1 — gene-level sequence alignment
  • Transcript match/mismatch pileup in rank-1 gene coordinates (from CS strings).
  • Needleman-Wunsch pairwise alignment of rank-1 vs rank-2 gene sequence, drawn
    as a per-position colour strip (green=identical, red=mismatch, grey=deletion,
    blue▲=insertion).  Mismatch positions are annotated with the actual base pair.

PAGE TYPE 2 — transcript-level alignment examples
  Answers "why is this gene frequent as rank-2 despite many sequence differences?"
  For up to N example transcripts, both alignments are drawn in *transcript*
  coordinates — only the region that aligns — so you can see:
    • Do rank-1 and rank-2 cover the same stretch of the transcript?
    • Where does each gene match / mismatch the transcript?
  A rank-2 gene can be frequent yet differ greatly from rank-1 if it aligns to
  a different (or shorter) region of the transcript at locally high identity.

CS string ops (minimap2 --cs):
  :N   → N matches       (advance both transcript and V gene)
  *xy  → mismatch        (advance both by 1)
  +seq → insertion in transcript  (advance transcript only)
  -seq → deletion in transcript   (advance V gene only, no transcript position consumed)
"""
import argparse
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd


# ─── Needleman-Wunsch global alignment ────────────────────────────────────────

def needleman_wunsch(seq1, seq2, match=2, mismatch=-1, gap=-2):
    """Return (aligned_seq1, aligned_seq2) as equal-length strings with '-' gaps."""
    n, m = len(seq1), len(seq2)
    dp = np.zeros((n + 1, m + 1), dtype=np.int32)
    dp[0, :] = np.arange(m + 1) * gap
    dp[:, 0] = np.arange(n + 1) * gap
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s = match if seq1[i-1] == seq2[j-1] else mismatch
            dp[i, j] = max(dp[i-1, j-1] + s,
                           dp[i-1, j]   + gap,
                           dp[i, j-1]   + gap)
    a1, a2 = [], []
    i, j = n, m
    while i > 0 or j > 0:
        if i > 0 and j > 0:
            s = match if seq1[i-1] == seq2[j-1] else mismatch
            if dp[i, j] == dp[i-1, j-1] + s:
                a1.append(seq1[i-1]); a2.append(seq2[j-1])
                i -= 1; j -= 1; continue
        if i > 0 and (j == 0 or dp[i, j] == dp[i-1, j] + gap):
            a1.append(seq1[i-1]); a2.append('-'); i -= 1
        else:
            a1.append('-'); a2.append(seq2[j-1]); j -= 1
    return ''.join(reversed(a1)), ''.join(reversed(a2))


def alignment_arrays(aligned_ref, aligned_qry):
    """Return per-ref-position arrays from a pairwise alignment string pair."""
    ref_bases, qry_bases = [], []
    match_arr, mm_arr, del_arr = [], [], []
    insertions = []   # list of (ref_pos_before, inserted_string)
    ins_buf = ''
    ref_pos = 0

    for rb, qb in zip(aligned_ref, aligned_qry):
        if rb == '-':
            ins_buf += qb
        else:
            if ins_buf:
                insertions.append((ref_pos, ins_buf))
                ins_buf = ''
            ref_bases.append(rb)
            if qb == '-':
                qry_bases.append('-')
                match_arr.append(False); mm_arr.append(False); del_arr.append(True)
            else:
                qry_bases.append(qb)
                same = rb.upper() == qb.upper()
                match_arr.append(same); mm_arr.append(not same); del_arr.append(False)
            ref_pos += 1

    if ins_buf:
        insertions.append((ref_pos, ins_buf))

    return (np.array(match_arr), np.array(mm_arr), np.array(del_arr),
            ''.join(ref_bases), qry_bases, insertions)


# ─── Sequence utilities ───────────────────────────────────────────────────────

_RC_TABLE = str.maketrans('ACGTacgtNn', 'TGCAtgcaNn')

def revcomp(seq):
    return seq.translate(_RC_TABLE)[::-1]


def dominant_strand(gene, transcript_ids, paf_index):
    """Return the most common alignment strand ('+' or '-') for gene across
    the given transcripts.  Ties go to '+'."""
    counts = {'+': 0, '-': 0}
    for tid in transcript_ids:
        rec = paf_index.get((tid, gene))
        if rec and rec.get('strand') in ('+', '-'):
            counts[rec['strand']] += 1
    return '+' if counts['+'] >= counts['-'] else '-'


# ─── FASTA reader ─────────────────────────────────────────────────────────────

def read_fasta(paths):
    seqs = {}
    for path in paths:
        current = None
        with open(path) as fh:
            for line in fh:
                line = line.rstrip()
                if line.startswith('>'):
                    current = line[1:].split()[0]
                    seqs[current] = []
                elif current:
                    seqs[current].append(line)
    return {k: ''.join(v) for k, v in seqs.items()}


# ─── CS string parsing ────────────────────────────────────────────────────────

_CS_RE = re.compile(r'(:[0-9]+|\*[a-z][a-z]|[+\-][a-z]+)')


def parse_cs(cs_tag):
    ops = []
    for token in _CS_RE.findall(cs_tag):
        if token[0] == ':':  ops.append(('M', int(token[1:])))
        elif token[0] == '*': ops.append(('X', 1))
        elif token[0] == '+': ops.append(('I', len(token) - 1))
        elif token[0] == '-': ops.append(('D', len(token) - 1))
    return ops


def cs_to_target_pileup(target_start, target_len, cs_ops, n_match, n_mismatch):
    """Accumulate match/mismatch in *target* (V gene) coordinates."""
    pos = target_start
    for op, val in cs_ops:
        if op == 'M':
            end = min(pos + val, target_len)
            if end > pos:
                n_match[pos:end] += 1
            pos += val
        elif op == 'X':
            if pos < target_len:
                n_mismatch[pos] += 1
            pos += 1
        elif op == 'D':
            pos += val
        if pos >= target_len:
            break


# STATUS values for transcript-space track arrays
_NOT_COVERED = 0
_MATCH       = 1
_MISMATCH    = 2   # substitution or insertion in transcript


def cs_to_query_track(query_start, query_end, query_len, cs_ops, strand='+'):
    """
    Return a per-transcript-position status array (length = query_len).

    Values: 0 = not covered, 1 = match to V gene, 2 = mismatch/insertion.

    For '+' strand the CS ops are processed left-to-right starting at
    query_start.  For '-' strand minimap2 emits the CS ops in *target*
    direction, which corresponds to right-to-left in forward-strand query
    coordinates: the first op maps to positions near query_end and subsequent
    ops work back towards query_start.  Ignoring this produces a mirror-image
    mismatch pattern — which is exactly the "flipped" artefact.

    CS op → query coordinate behaviour
      :N   match N bases        → advance by N (fwd) / retreat by N (rev)
      *xy  1 mismatch           → advance/retreat by 1
      +seq insertion in query   → advance/retreat by len(seq); marked mismatch
      -seq deletion in query    → target advances, query position unchanged
    """
    status = np.zeros(query_len, dtype=np.int8)

    if strand == '+':
        q = query_start
        for op, val in cs_ops:
            if op == 'M':
                end = min(q + val, query_len)
                if end > q:
                    status[q:end] = _MATCH
                q += val
            elif op == 'X':
                if q < query_len:
                    status[q] = _MISMATCH
                q += 1
            elif op == 'I':
                end = min(q + val, query_len)
                if end > q:
                    status[q:end] = _MISMATCH
                q += val
            # 'D': target advances, query position unchanged
            if q >= query_len:
                break

    else:  # strand == '-': process right-to-left in forward-strand coordinates
        q = query_end - 1   # start at the rightmost covered query position
        for op, val in cs_ops:
            if op == 'M':
                start = max(q - val + 1, 0)
                if start <= q:
                    status[start:q + 1] = _MATCH
                q -= val
            elif op == 'X':
                if 0 <= q < query_len:
                    status[q] = _MISMATCH
                q -= 1
            elif op == 'I':
                start = max(q - val + 1, 0)
                if start <= q:
                    status[start:q + 1] = _MISMATCH
                q -= val
            # 'D': target advances, query position unchanged
            if q < 0:
                break

    return status


# ─── PAF index ────────────────────────────────────────────────────────────────

def build_paf_index(paf_paths):
    """(transcript, gene) → best PAF record dict (highest n_matches).

    Stores both query (transcript) and target (V gene) coordinates so the
    index can serve both the gene-coordinate pileup and the transcript-space
    alignment examples.
    """
    index = {}
    for path in paf_paths:
        with open(path) as fh:
            for line in fh:
                if line.startswith('#') or not line.strip():
                    continue
                f = line.rstrip().split('\t')
                cs_tag = None
                for field in f[12:]:
                    if field.startswith('cs:Z:'):
                        cs_tag = field[5:]
                        break
                key = (f[0], f[5])
                nm  = int(f[9])
                if key not in index or nm > index[key]['n_matches']:
                    index[key] = {
                        'query_len':    int(f[1]),
                        'query_start':  int(f[2]),
                        'query_end':    int(f[3]),
                        'strand':       f[4],
                        'target_len':   int(f[6]),
                        'target_start': int(f[7]),
                        'target_end':   int(f[8]),
                        'n_matches':    nm,
                        'aln_block':    int(f[10]),
                        'cs':           cs_tag,
                    }
    return index


# ─── Per-gene pileup (target/gene coordinates) ───────────────────────────────

def build_gene_pileup(gene, transcript_ids, paf_index):
    tlen = None
    for tid in transcript_ids:
        rec = paf_index.get((tid, gene))
        if rec:
            tlen = rec['target_len']
            break
    if tlen is None:
        return np.zeros(1), np.zeros(1)

    n_match    = np.zeros(tlen, dtype=np.float64)
    n_mismatch = np.zeros(tlen, dtype=np.float64)
    for tid in transcript_ids:
        rec = paf_index.get((tid, gene))
        if rec is None or rec['cs'] is None:
            continue
        cs_to_target_pileup(rec['target_start'], tlen,
                            parse_cs(rec['cs']), n_match, n_mismatch)
    return n_match, n_mismatch


# ─── Top-alignment loader ─────────────────────────────────────────────────────

def load_top_alns(paths, sample_labels):
    frames = []
    for path, label in zip(paths, sample_labels):
        try:
            df = pd.read_csv(path, sep='\t')
            df['sample'] = label
            frames.append(df)
        except Exception as e:
            print(f'  Warning: {path}: {e}', file=sys.stderr)
    if not frames:
        raise RuntimeError('No top-alignment TSV files loaded.')
    return pd.concat(frames, ignore_index=True)


# ─── Colour constants ─────────────────────────────────────────────────────────

C_MATCH  = np.array([0.18, 0.63, 0.18])   # green
C_MM     = np.array([0.84, 0.15, 0.15])   # red
C_DEL    = np.array([0.80, 0.80, 0.80])   # light grey
C_INS    = '#1f77b4'

# Transcript-space colours: rank-1 green/red, rank-2 blue/orange
_TRACK_COLORS = {
    'r1': {_NOT_COVERED: [0.96, 0.96, 0.96],
           _MATCH:       [0.18, 0.63, 0.18],
           _MISMATCH:    [0.84, 0.15, 0.15]},
    'r2': {_NOT_COVERED: [0.96, 0.96, 0.96],
           _MATCH:       [0.12, 0.47, 0.71],
           _MISMATCH:    [1.00, 0.50, 0.05]},
}


# ─── Page type 1: gene-level sequence alignment ───────────────────────────────

def _draw_alignment_strip(ax, match_arr, mm_arr, del_arr,
                          ref_bases, qry_bases, insertions,
                          gene_len, partner_name, n_transcripts, seq_id_pct,
                          orientation_note=''):
    L   = len(match_arr)
    img = np.ones((1, L, 3), dtype=np.float32)
    img[0, match_arr] = C_MATCH
    img[0, mm_arr]    = C_MM
    img[0, del_arr]   = C_DEL
    ax.imshow(img, aspect='auto', extent=[0, L, 0, 1],
              interpolation='nearest', zorder=1)

    for ref_pos, _ in insertions:
        ax.annotate('▲', xy=(ref_pos, 1), xytext=(ref_pos, 1.25),
                    fontsize=5, ha='center', va='bottom', color=C_INS,
                    annotation_clip=False)

    for pos in np.where(mm_arr)[0]:
        ax.text(pos + 0.5, -0.35, f'{ref_bases[pos]}\n{qry_bases[pos]}',
                fontsize=4.5, ha='center', va='top',
                fontfamily='monospace', color='#8b0000', clip_on=False)

    ax.set_xlim(0, L)
    ax.set_ylim(-0.6, 1.3)
    ax.set_yticks([])
    ax.tick_params(axis='x', labelsize=6)
    n_mm  = int(mm_arr.sum())
    n_del = int(del_arr.sum())
    ax.set_title(
        f'{partner_name}  '
        f'[seq identity {seq_id_pct:.1f}%,  {n_mm} mismatches,  {n_del} deletions,  '
        f'n={n_transcripts} transcripts as rank-2]{orientation_note}',
        fontsize=7.5, pad=8)
    return [mpatches.Patch(color=C_MATCH, label='identical'),
            mpatches.Patch(color=C_MM,   label='mismatch'),
            mpatches.Patch(color=C_DEL,  label='deletion in rank-2'),
            mpatches.Patch(color=C_INS,  label='insertion in rank-2 (▲)')]


def plot_gene_page(pdf, rank1_gene, locus, sequences,
                   all_r1_tids, partners, paf_index):
    n_partners = len(partners)
    fig = plt.figure(figsize=(14, 2.8 + n_partners * 1.6 + 0.6))
    gs  = gridspec.GridSpec(1 + n_partners, 1, figure=fig,
                            height_ratios=[2.8] + [1.6] * n_partners,
                            hspace=0.65, top=0.93, bottom=0.06)

    ref_seq  = sequences.get(rank1_gene, '')
    gene_len = len(ref_seq)

    # Transcript pileup
    ax0 = fig.add_subplot(gs[0])
    match, mismatch = build_gene_pileup(rank1_gene, all_r1_tids, paf_index)
    x = np.arange(len(match))
    ax0.fill_between(x, 0, match + mismatch, step='mid', color='#cccccc')
    ax0.fill_between(x, 0, match,            step='mid', color='#2ca02c', alpha=0.9, label='match')
    ax0.fill_between(x, match, match + mismatch, step='mid', color='#d62728', alpha=0.9, label='mismatch')
    ax0.set_xlim(0, gene_len)
    ax0.set_ylabel('# transcripts', fontsize=7)
    ax0.tick_params(labelsize=6)
    ax0.set_title(
        f'[{locus}]  {rank1_gene}  —  transcript pileup  '
        f'(n={len(all_r1_tids)} transcripts where this is rank-1)',
        fontsize=9, fontweight='bold')
    ax0.legend(fontsize=7, loc='upper right', ncol=3)

    r1_strand = dominant_strand(rank1_gene, all_r1_tids, paf_index)

    # Alignment strips
    for row, (partner, tids) in enumerate(partners, start=1):
        ax = fig.add_subplot(gs[row])
        qry_seq = sequences.get(partner, '')
        if not ref_seq or not qry_seq:
            ax.text(0.5, 0.5, f'Sequence not found for {partner}',
                    ha='center', va='center', transform=ax.transAxes)
            continue

        # Determine orientation: if the two genes align on opposite strands to
        # transcripts, their raw sequences are in opposite orientations — RC the
        # partner before NW so we compare what actually faces the transcript.
        r2_strand = dominant_strand(partner, tids, paf_index)
        if r1_strand != r2_strand:
            # RC whichever gene is on the minus strand
            qry_seq_aln = revcomp(qry_seq) if r2_strand == '-' else qry_seq
            ref_seq_aln = revcomp(ref_seq) if r1_strand == '-' else ref_seq
            orientation_note = (f'\n  ⚠ opposite strands (rank-1:{r1_strand} '
                                f'rank-2:{r2_strand}) — rank-2 reverse-complemented '
                                f'for comparison')
        else:
            ref_seq_aln, qry_seq_aln = ref_seq, qry_seq
            orientation_note = f'  (both align {r1_strand} strand)'

        a_ref, a_qry = needleman_wunsch(ref_seq_aln, qry_seq_aln)
        m_arr, mm_arr, d_arr, ref_b, qry_b, ins = alignment_arrays(a_ref, a_qry)
        seq_id = m_arr.sum() / len(m_arr) * 100 if len(m_arr) > 0 else 0.0
        patches = _draw_alignment_strip(ax, m_arr, mm_arr, d_arr,
                                        ref_b, qry_b, ins,
                                        gene_len, partner, len(tids), seq_id,
                                        orientation_note=orientation_note)
        if row == 1:
            ax.legend(handles=patches, fontsize=6.5, loc='upper right', ncol=4)
        if row == n_partners:
            ax.set_xlabel(f'Position in {rank1_gene}  (bp)', fontsize=8)

    fig.suptitle(
        f'Gene-level sequence alignment: {rank1_gene} (rank-1) vs top rank-2 partners',
        fontsize=8.5, y=0.99)
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)


# ─── Page type 2: transcript-level alignment examples ─────────────────────────

def _short_name(gene):
    """Return last meaningful part of a long gene ID for axis labels."""
    parts = gene.split('.')
    # e.g. tufted_duck_IGH.229270.WNMM01000072.1.V.True.-
    # return something like "IGH·229270"
    for i, p in enumerate(parts):
        if p.startswith('WNMM') or p.startswith('NW_') or p.startswith('NC_'):
            return '.'.join(parts[:i])[-30:]
    return gene[-30:]


def plot_transcript_examples(pdf, rank1_gene, partner_gene, locus,
                              candidate_tids, paf_index, aln_df,
                              n_examples=5):
    """
    One PDF page: up to n_examples transcripts showing both alignments
    drawn in transcript (query) coordinates.
    """
    # Keep only transcripts that have PAF records for both genes
    usable = [t for t in candidate_tids
              if (t, rank1_gene) in paf_index and (t, partner_gene) in paf_index]
    if not usable:
        return

    # Sort by rank-2 identity descending so we show the most interesting cases first
    r2_id_lookup = (aln_df[aln_df['gene'] == partner_gene]
                    .set_index('transcript')['identity']
                    .to_dict())
    usable.sort(key=lambda t: r2_id_lookup.get(t, 0.0), reverse=True)
    examples = usable[:n_examples]

    n = len(examples)
    row_h = 1.8
    fig   = plt.figure(figsize=(14, n * row_h + 1.4))
    gs    = gridspec.GridSpec(n, 1, figure=fig,
                              hspace=0.75, top=0.91, bottom=0.06)

    r1_short = _short_name(rank1_gene)
    r2_short = _short_name(partner_gene)

    for row, tid in enumerate(examples):
        ax  = fig.add_subplot(gs[row])
        rec1 = paf_index[(tid, rank1_gene)]
        rec2 = paf_index[(tid, partner_gene)]
        qlen = max(rec1['query_len'], rec2['query_len'])

        # Build per-position status tracks, respecting alignment strand
        track1 = np.zeros(qlen, dtype=np.int8)
        track2 = np.zeros(qlen, dtype=np.int8)
        if rec1['cs']:
            track1 = cs_to_query_track(rec1['query_start'], rec1['query_end'],
                                        qlen, parse_cs(rec1['cs']),
                                        strand=rec1['strand'])
        if rec2['cs']:
            track2 = cs_to_query_track(rec2['query_start'], rec2['query_end'],
                                        qlen, parse_cs(rec2['cs']),
                                        strand=rec2['strand'])

        # Crop to the region covered by at least one alignment + small padding
        pad  = 10
        x_lo = max(0,    min(rec1['query_start'], rec2['query_start']) - pad)
        x_hi = min(qlen, max(rec1['query_end'],   rec2['query_end'])   + pad)
        w    = x_hi - x_lo

        # Build 2-row RGB image  (row 0 = rank-1, row 1 = rank-2)
        img = np.ones((2, w, 3), dtype=np.float32)
        for i in range(w):
            qi = i + x_lo
            img[0, i] = _TRACK_COLORS['r1'][track1[qi]]
            img[1, i] = _TRACK_COLORS['r2'][track2[qi]]

        ax.imshow(img, aspect='auto', extent=[x_lo, x_hi, 0, 2],
                  interpolation='nearest', zorder=1)
        ax.axhline(1.0, color='white', linewidth=0.8, zorder=2)

        # Dashed lines at alignment boundaries
        for start, end, col in [
            (rec1['query_start'], rec1['query_end'], '#2ca02c'),
            (rec2['query_start'], rec2['query_end'], '#1f77b4'),
        ]:
            for xv in (start, end):
                ax.axvline(xv, color=col, linewidth=0.8,
                           linestyle='--', alpha=0.8, zorder=3)

        # Y-axis labels
        ax.set_yticks([0.5, 1.5])
        ax.set_yticklabels([r2_short, r1_short], fontsize=5.5)
        ax.set_ylim(0, 2)
        ax.set_xlim(x_lo, x_hi)
        ax.tick_params(axis='x', labelsize=6)

        # Per-transcript stats from the TSV
        def _stats(gene):
            row_df = aln_df[(aln_df['transcript'] == tid) &
                            (aln_df['gene'] == gene)]
            if row_df.empty:
                return '—'
            r = row_df.iloc[0]
            return f'id={r["identity"]:.3f}  cov={r["v_coverage"]:.2f}'

        s1, s2 = rec1['strand'], rec2['strand']
        strand_note = (f'  ⚠ opposite strands (rank-1:{s1} rank-2:{s2})'
                       if s1 != s2 else f'  (both {s1})')
        ax.set_title(
            f'{tid}  (len={qlen} bp){strand_note}  |  '
            f'rank-1 [{r1_short}]: {_stats(rank1_gene)}  |  '
            f'rank-2 [{r2_short}]: {_stats(partner_gene)}',
            fontsize=6.5, pad=3)

        if row == n - 1:
            ax.set_xlabel('Transcript position (bp)', fontsize=7)

    # Shared legend
    legend_patches = [
        mpatches.Patch(color=_TRACK_COLORS['r1'][_MATCH],    label='rank-1 match'),
        mpatches.Patch(color=_TRACK_COLORS['r1'][_MISMATCH], label='rank-1 mismatch'),
        mpatches.Patch(color=_TRACK_COLORS['r2'][_MATCH],    label='rank-2 match'),
        mpatches.Patch(color=_TRACK_COLORS['r2'][_MISMATCH], label='rank-2 mismatch'),
        mpatches.Patch(color=_TRACK_COLORS['r1'][_NOT_COVERED], label='not covered'),
    ]
    fig.legend(handles=legend_patches, loc='upper right',
               fontsize=7, ncol=5, bbox_to_anchor=(0.99, 0.99))

    fig.suptitle(
        f'Transcript alignments — {locus}  |  '
        f'rank-1: {rank1_gene}  vs  rank-2: {partner_gene}\n'
        f'Dashed lines mark alignment boundaries.  '
        f'Showing {n} example transcripts (highest rank-2 identity first).',
        fontsize=8, y=0.995)
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)


# ─── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--top-alns',     nargs='+', required=True)
    parser.add_argument('--pafs',         nargs='+', required=True)
    parser.add_argument('--vgene-fastas', nargs='+', required=True)
    parser.add_argument('--samples',      nargs='+', required=True)
    parser.add_argument('--output',       required=True)
    parser.add_argument('--top-genes',    type=int, default=5)
    parser.add_argument('--top-partners', type=int, default=3)
    parser.add_argument('--n-examples',   type=int, default=5,
                        help='Example transcripts per rank-2 partner (default: 5)')
    args = parser.parse_args()

    print('Reading V gene sequences…', file=sys.stderr)
    sequences = read_fasta(args.vgene_fastas)
    print(f'  {len(sequences)} sequences loaded', file=sys.stderr)

    print('Loading top-alignment tables…', file=sys.stderr)
    aln = load_top_alns(args.top_alns, args.samples)

    print('Indexing PAF files…', file=sys.stderr)
    paf_index = build_paf_index(args.pafs)
    print(f'  {len(paf_index):,} (transcript, gene) records indexed', file=sys.stderr)

    best_per_locus = (
        aln.sort_values('identity', ascending=False)
           .groupby(['sample', 'transcript', 'locus'], as_index=False)
           .first()
    )

    with PdfPages(args.output) as pdf:
        for locus in sorted(aln['locus'].unique()):
            r1 = best_per_locus[best_per_locus['locus'] == locus]
            top_genes = (r1.groupby('gene')['transcript']
                           .nunique()
                           .sort_values(ascending=False)
                           .head(args.top_genes)
                           .index.tolist())

            print(f'\n{locus}: {len(top_genes)} rank-1 genes', file=sys.stderr)

            for rank1_gene in top_genes:
                all_r1_tids = set(r1[r1['gene'] == rank1_gene]['transcript'])

                # Best non-rank1 hit per transcript for this locus
                r2 = (aln[aln['transcript'].isin(all_r1_tids) &
                          (aln['locus'] == locus) &
                          (aln['gene'] != rank1_gene)]
                      .sort_values('identity', ascending=False)
                      .groupby(['sample', 'transcript'], as_index=False)
                      .first())

                partners = [(p, list(r2[r2['gene'] == p]['transcript'].unique()))
                            for p in (r2.groupby('gene')['transcript']
                                        .nunique()
                                        .sort_values(ascending=False)
                                        .head(args.top_partners)
                                        .index)]

                if not partners:
                    continue

                print(f'  {rank1_gene}: {len(all_r1_tids)} transcripts, '
                      f'{len(partners)} partners', file=sys.stderr)

                # Page 1: gene-level NW alignment
                plot_gene_page(pdf, rank1_gene, locus, sequences,
                               list(all_r1_tids), partners, paf_index)

                # Page 2+: transcript-level examples, one page per partner
                for partner, tids in partners:
                    plot_transcript_examples(pdf, rank1_gene, partner, locus,
                                             tids, paf_index, aln,
                                             n_examples=args.n_examples)

    print(f'\nWritten: {args.output}', file=sys.stderr)


if __name__ == '__main__':
    main()
