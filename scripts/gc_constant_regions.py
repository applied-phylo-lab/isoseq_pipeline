"""
Pull the constant (C) region out of each IG transcript and group the classes.

Where the C region starts
-------------------------
A full-length IsoSeq IG transcript reads V - D - J - C.  The V end is known from
the V gene alignment, but the V/C distance is not fixed: J is short and D adds a
variable amount, so cutting at "V end + a constant" lands somewhere different in
every transcript.  Instead the J sequence is located directly by local alignment
and the C region is taken as everything 3' of the J match.  Transcripts where J
cannot be found are reported rather than silently trimmed at a guess.

Grouping without a reference
----------------------------
There is no bird constant-region database here, so classes are recovered from
the data: C regions are clustered by k-mer similarity, which separates isotypes
easily because different constant genes are far less similar to each other than
alleles of one are.  Clusters are then characterised by what can be read off the
sequence itself -- length, translation, and the number of immunoglobulin
domains, counted from the spacing of the conserved cysteine pairs that define
them.

Naming the classes needs one more thing
---------------------------------------
Cluster identity (which one is IgM, IgY, IgA) cannot be settled by clustering
alone.  Two handles are written out for that: the per-cluster CH domain count,
and -- with --assembly -- the genomic position of each cluster downstream of the
V array, since the constant genes sit in a diagnostic order along the locus.
"""
import argparse
import csv
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict

from Bio import Align
from Bio.Seq import Seq

from gc_lib import read_fasta, parse_paf

RC = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def revcomp(s):
    return s[::-1].translate(RC)


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def kmers(s, k=12):
    return {s[i:i + k] for i in range(len(s) - k + 1)}


def jaccard(a, b):
    return len(a & b) / max(1, len(a | b))


def find_j_end(seq, jseq, aligner, min_score):
    """3' end of the best J match, or None. Both orientations are tried."""
    best = None
    for cand, rc in ((seq, False), (revcomp(seq), True)):
        try:
            aln = aligner.align(cand, jseq)[0]
        except (ValueError, IndexError):
            continue
        if aln.score < min_score:
            continue
        end = max(b for a, b in aln.aligned[0]) if len(aln.aligned[0]) else None
        if end is None:
            continue
        if best is None or aln.score > best[0]:
            best = (aln.score, end, rc)
    return best


def ig_domains(aa):
    """Count immunoglobulin domains from conserved cysteine spacing.

    An Ig fold is pinned by two cysteines ~55-75 residues apart with a tryptophan
    ~12-18 residues after the first. Counting those is crude but it is a property
    of the sequence rather than of a reference, which is the point here.
    """
    cys = [i for i, c in enumerate(aa) if c == "C"]
    n, used = 0, set()
    for i in cys:
        if i in used:
            continue
        for j in cys:
            if 55 <= j - i <= 80 and j not in used:
                w = aa[i + 10:i + 22]
                if "W" in w:
                    n += 1
                    used |= {i, j}
                    break
    return n


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--transcripts", required=True)
    ap.add_argument("--paf", required=True, help="detailed V gene PAF for this locus")
    ap.add_argument("--jgene", required=True)
    ap.add_argument("--locus", required=True)
    ap.add_argument("--min-c-len", type=int, default=200)
    ap.add_argument("--cluster-threshold", type=float, default=0.30)
    ap.add_argument("--assembly", help="optional: locate each cluster in the genome")
    ap.add_argument("--out-fasta", required=True)
    ap.add_argument("--out-table", required=True)
    args = ap.parse_args()

    tx = read_fasta(args.transcripts)
    jseq = next(iter(read_fasta(args.jgene).values())).upper()

    # local alignment: J is a short island inside a long transcript
    al = Align.PairwiseAligner(mode="local", match_score=2, mismatch_score=-1,
                              open_gap_score=-6, extend_gap_score=-1)
    min_score = 1.1 * len(jseq)          # ~55% of a perfect local match

    vend = {}
    for rec in parse_paf(args.paf):
        if rec.query not in tx:
            continue
        cur = vend.get(rec.query)
        if cur is None or rec.nmatch > cur[1]:
            vend[rec.query] = (rec.qend, rec.nmatch)

    cregions, no_j, too_short = {}, [], []
    for q in sorted(vend):
        seq = tx[q].upper()
        hit = find_j_end(seq, jseq, al, min_score)
        if hit is None:
            no_j.append(q)
            continue
        _score, jend, rc = hit
        c = (revcomp(seq) if rc else seq)[jend:]
        if len(c) < args.min_c_len:
            too_short.append(q)
            continue
        cregions[q] = c

    print(f"{args.locus}: {len(vend)} transcripts with a V alignment", file=sys.stderr)
    print(f"  J found and C region >= {args.min_c_len} bp : {len(cregions)}", file=sys.stderr)
    print(f"  J not found                                : {len(no_j)}", file=sys.stderr)
    print(f"  J found but C too short                    : {len(too_short)}", file=sys.stderr)
    if not cregions:
        raise SystemExit("no constant regions recovered")

    ks = {q: kmers(c) for q, c in cregions.items()}
    clusters = []
    for q in sorted(cregions, key=lambda x: -len(cregions[x])):
        for cl in clusters:
            if jaccard(ks[q], ks[cl[0]]) >= args.cluster_threshold:
                cl.append(q)
                break
        else:
            clusters.append([q])
    clusters.sort(key=len, reverse=True)

    with open(args.out_fasta, "w") as fh:
        for i, cl in enumerate(clusters, 1):
            for q in cl:
                fh.write(f">{q} class={i} len={len(cregions[q])}\n{cregions[q]}\n")

    rows = []
    for i, cl in enumerate(clusters, 1):
        rep = max(cl, key=lambda q: len(cregions[q]))
        seq = cregions[rep]
        # frame with the fewest stops; a C region is one long ORF
        best = min(range(3), key=lambda f: str(
            Seq(seq[f:len(seq) - ((len(seq) - f) % 3)]).translate()).count("*"))
        aa = str(Seq(seq[best:len(seq) - ((len(seq) - best) % 3)]).translate())
        lens = sorted(len(cregions[q]) for q in cl)
        rows.append({
            "locus": args.locus, "class": i, "n_transcripts": len(cl),
            "pct": f"{len(cl) / len(cregions):.1%}",
            "median_len": lens[len(lens) // 2], "max_len": lens[-1],
            "ig_domains_in_longest": ig_domains(aa),
            "representative": rep, "protein": aa,
        })

    # between-class identity: isotypes are far apart, alleles of one are not
    print("\n  class  n      %   median bp  max bp  Ig domains", file=sys.stderr)
    for r in rows:
        print(f"  {r['class']:>5}  {r['n_transcripts']:<5} {r['pct']:>5}  "
              f"{r['median_len']:>9}  {r['max_len']:>6}  {r['ig_domains_in_longest']:>10}",
              file=sys.stderr)
        print(f"        {r['protein'][:96]}", file=sys.stderr)
    if len(rows) > 1:
        print("\n  pairwise k-mer similarity between classes:", file=sys.stderr)
        for a in range(len(rows)):
            for b in range(a + 1, len(rows)):
                ja = jaccard(kmers(cregions[rows[a]["representative"]]),
                             kmers(cregions[rows[b]["representative"]]))
                print(f"    class {rows[a]['class']} vs {rows[b]['class']}: {ja:.3f}",
                      file=sys.stderr)

    if args.assembly:
        # genomic order of the constant genes is diagnostic in birds, so each
        # class representative is placed on the assembly
        with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as fh:
            for r in rows:
                fh.write(f">class{r['class']}\n{cregions[r['representative']]}\n")
            qpath = fh.name
        try:
            out = subprocess.run(
                ["/programs/bin/blast+/blastn", "-query", qpath,
                 "-subject", args.assembly, "-outfmt",
                 "6 qseqid sseqid sstart send pident length bitscore",
                 "-max_target_seqs", "5"],
                capture_output=True, text=True, timeout=1800).stdout
            print("\n  genomic placement (best hits):", file=sys.stderr)
            seen = set()
            for ln in out.splitlines():
                f = ln.split("\t")
                if f[0] in seen or float(f[5]) < 200:
                    continue
                seen.add(f[0])
                print(f"    {f[0]:>7}  {f[1]}:{f[2]}-{f[3]}  "
                      f"{f[4]}% over {f[5]} bp", file=sys.stderr)
                for r in rows:
                    if f"class{r['class']}" == f[0]:
                        r["contig"], r["gstart"], r["gend"] = f[1], f[2], f[3]
        except Exception as exc:
            print(f"  assembly placement failed: {exc}", file=sys.stderr)

    cols = ["locus", "class", "n_transcripts", "pct", "median_len", "max_len",
            "ig_domains_in_longest", "contig", "gstart", "gend",
            "representative", "protein"]
    with open(args.out_table, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "NA")) for c in cols) + "\n")


if __name__ == "__main__":
    main()
