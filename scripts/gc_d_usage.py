"""
Work out which IGHD genes are used, from the transcripts.

Why this is harder than V or J
------------------------------
A D gene enters the transcript already damaged.  Exonuclease chews both ends
during V-D and D-J joining and N nucleotides are added on either side, so what
survives is an internal fragment of a gene that was only 20-66 bp to begin
with.  Median usable D remnants in vertebrate repertoires are often under
10 bp.

That creates a statistical problem rather than an alignment problem: a 6 bp
match to one of ~20 D genes happens by chance all the time.  Assigning the
best-scoring D and stopping would produce a confident-looking D usage profile
out of pure noise.

So every call here is scored against a null
-------------------------------------------
For each transcript the junction region (between the end of V and the start of
J) is extracted, and the longest exact match to any D gene is found, in both
orientations.  The same search is then repeated on shuffled versions of that
same junction -- preserving its length and base composition -- to get the match
length expected by chance for that particular junction.  A D call is reported
as confident only when the observed match is longer than the null essentially
always (empirical p below --p-threshold).

Ambiguity is reported rather than hidden: where several D genes tie at the best
match length, all are listed, because with genes this short ties are common and
picking one arbitrarily invents precision.
"""
import argparse
import random
import statistics
import sys
from collections import Counter, defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from gc_lib import read_fasta, parse_paf
from gc_palette import save_figure, LOCUS, SEGMENT, GREY, GREY_DARK, YES


def rc(s):
    return s.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def longest_shared(a, b):
    """Length of the longest substring common to a and b, and where it sits in a."""
    if not a or not b:
        return 0, -1
    best, best_i = 0, -1
    prev = [0] * (len(b) + 1)
    for i in range(1, len(a) + 1):
        cur = [0] * (len(b) + 1)
        for j in range(1, len(b) + 1):
            if a[i - 1] == b[j - 1]:
                cur[j] = prev[j - 1] + 1
                if cur[j] > best:
                    best, best_i = cur[j], i - cur[j]
        prev = cur
    return best, best_i


def best_d(junction, dgenes):
    """Longest match to any D, in either orientation. Returns (len, [names], pos)."""
    best_len, hits, pos = 0, [], -1
    for name, seq in dgenes.items():
        for oriented in (seq, rc(seq)):
            L, i = longest_shared(junction, oriented)
            if L > best_len:
                best_len, hits, pos = L, [name], i
            elif L == best_len and L > 0 and name not in hits:
                hits.append(name)
    return best_len, hits, pos


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--transcripts", required=True)
    ap.add_argument("--paf", required=True, help="transcripts vs V genes")
    ap.add_argument("--assignments", required=True)
    ap.add_argument("--dgene-fasta", required=True)
    ap.add_argument("--jgene-fasta", required=True)
    ap.add_argument("--locus", default="IGH")
    ap.add_argument("--min-match", type=int, default=6,
                    help="shortest D remnant considered at all (default: 6)")
    ap.add_argument("--p-threshold", type=float, default=0.05)
    ap.add_argument("--permutations", type=int, default=200)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out-per-transcript", required=True)
    ap.add_argument("--out-usage", required=True)
    ap.add_argument("--out-figure", required=True)
    args = ap.parse_args()

    rng = random.Random(args.seed)
    tx = read_fasta(args.transcripts)
    dgenes = read_fasta(args.dgene_fasta)
    jseq = next(iter(read_fasta(args.jgene_fasta).values()))

    parent_of = {}
    with open(args.assignments) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            r = dict(zip(hdr, line.rstrip("\n").split("\t")))
            if r["locus"] == args.locus:
                parent_of[r["transcript"]] = r["best_gene"]

    # end of the V alignment in transcript coordinates
    v_end = {}
    for rec in parse_paf(args.paf):
        if parent_of.get(rec.query) != rec.target:
            continue
        if rec.qend > v_end.get(rec.query, -1):
            v_end[rec.query] = rec.qend

    def find_j(seq):
        """Start of the J in the transcript: exact, else best 12-mer anchor."""
        i = seq.find(jseq)
        if i != -1:
            return i
        for k in range(0, len(jseq) - 12):
            i = seq.find(jseq[k:k + 12])
            if i != -1:
                return max(0, i - k)
        return -1

    rows, junction_lens, obs_lens, null_lens = [], [], [], []
    for tid, ve in sorted(v_end.items()):
        s = tx.get(tid)
        if s is None:
            continue
        js = find_j(s)
        if js < 0 or js <= ve:
            continue
        junction = s[ve:js]
        junction_lens.append(len(junction))
        if len(junction) < args.min_match:
            rows.append((tid, len(junction), 0, "", "no_junction", 1.0))
            continue
        L, hits, _ = best_d(junction, dgenes)
        # null: same junction, shuffled, same D set
        pool = list(junction)
        null = []
        for _ in range(args.permutations):
            rng.shuffle(pool)
            nl, _h, _p = best_d("".join(pool), dgenes)
            null.append(nl)
        p = (sum(1 for x in null if x >= L) + 1) / (len(null) + 1)
        call = ("confident" if (L >= args.min_match and p <= args.p_threshold)
                else ("weak" if L >= args.min_match else "no_call"))
        rows.append((tid, len(junction), L, ";".join(sorted(hits)), call, p))
        obs_lens.append(L)
        null_lens.append(statistics.mean(null) if null else 0)

    with open(args.out_per_transcript, "w") as fh:
        fh.write("transcript\tjunction_bp\tmatch_bp\td_genes\tcall\tp_value\n")
        for tid, jl, L, hits, call, p in rows:
            fh.write(f"{tid}\t{jl}\t{L}\t{hits or 'NA'}\t{call}\t{p:.4g}\n")

    usage = Counter()
    for tid, jl, L, hits, call, p in rows:
        if call == "confident" and hits:
            names = hits.split(";")
            for n in names:
                usage[n] += 1 / len(names)     # split credit across ties

    with open(args.out_usage, "w") as fh:
        fh.write("d_gene\tpos\tweighted_transcripts\n")
        for name, c in sorted(usage.items(), key=lambda kv: -kv[1]):
            fh.write(f"{name}\t{name.split('.')[1]}\t{c:.2f}\n")

    counts = Counter(r[4] for r in rows)
    n_conf = counts["confident"]

    fig, axes = plt.subplots(1, 4, figsize=(17, 4.0))
    C = LOCUS.get(args.locus, SEGMENT["D"])

    ax = axes[0]
    ax.hist(junction_lens, bins=range(0, max(junction_lens or [1]) + 3, 2),
            color=GREY_DARK, edgecolor="black", lw=.4)
    ax.set_xlabel("junction length V-end → J-start (bp)")
    ax.set_ylabel("transcripts")
    ax.set_title(f"A  junction length\nmedian {statistics.median(junction_lens or [0]):.0f} bp",
                 fontsize=10, fontweight="bold", loc="left")

    ax = axes[1]
    if obs_lens:
        mx = max(obs_lens + null_lens)
        bins = np.arange(0, mx + 2, 1)
        ax.hist(null_lens, bins=bins, color=GREY, edgecolor="black", lw=.4,
                label="expected by chance")
        ax.hist(obs_lens, bins=bins, color=C, alpha=.85, edgecolor="black", lw=.5,
                label="observed")
        ax.legend(fontsize=8)
    ax.set_xlabel("longest D match (bp)")
    ax.set_ylabel("transcripts")
    ax.set_title("B  is the match longer than chance?", fontsize=10,
                 fontweight="bold", loc="left")

    ax = axes[2]
    order = ["confident", "weak", "no_call", "no_junction"]
    vals = [counts.get(k, 0) for k in order]
    ax.bar(order, vals, color=[YES, GREY_DARK, GREY, GREY], edgecolor="black", lw=.6)
    for i, v in enumerate(vals):
        ax.text(i, v, f"\n{v}", ha="center", va="bottom", fontsize=9)
    ax.set_ylabel("transcripts")
    ax.set_ylim(0, max(vals + [1]) * 1.3)
    ax.tick_params(axis="x", labelsize=8)
    ax.set_title(f"C  call confidence\n{n_conf} of {len(rows)} assignable",
                 fontsize=10, fontweight="bold", loc="left")

    ax = axes[3]
    if usage:
        top = sorted(usage.items(), key=lambda kv: -kv[1])[:15][::-1]
        ax.barh(range(len(top)), [v for _, v in top], color=C,
                edgecolor="black", lw=.4)
        ax.set_yticks(range(len(top)))
        ax.set_yticklabels([k.split(".")[1] for k, _ in top], fontsize=7)
        ax.set_xlabel("transcripts (ties split)")
    else:
        ax.text(.5, .5, "no confident D calls", ha="center", transform=ax.transAxes)
    ax.set_title("D  D gene usage", fontsize=10, fontweight="bold", loc="left")

    for a in axes:
        for s in ("top", "right"):
            a.spines[s].set_visible(False)
    fig.suptitle(f"{args.locus} — D gene usage from the transcripts",
                 fontsize=12, fontweight="bold", y=1.04)
    fig.tight_layout()
    save_figure(fig, args.out_figure)

    print(f"{args.locus}: {len(rows)} transcripts with a junction; "
          f"{n_conf} confident D calls, {counts.get('weak',0)} weak, "
          f"{counts.get('no_call',0)} none\n"
          f"  junction median {statistics.median(junction_lens or [0]):.0f} bp; "
          f"observed match median {statistics.median(obs_lens or [0]):.0f} bp "
          f"vs chance {statistics.median(null_lens or [0]):.1f} bp",
          file=sys.stderr)


if __name__ == "__main__":
    main()
