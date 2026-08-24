"""
Collapse V gene annotations that are the same DNA entered twice.

The problem
-----------
An annotator that scans both strands emits a hit on each.  The result is pairs of
"V genes" whose genomic intervals OVERLAP and which, over the shared block, are
byte-identical on the forward strand -- one locus, written down twice, once as a
plus-strand gene and once as a minus-strand gene.  In bAgePho2 IGH there are 59
such pairs among 162 entries; IGL has one.

This is NOT the tandem duplication that fills a V array.  Tandem duplicates sit
at DIFFERENT coordinates and are 90-99% similar.  These sit at the SAME
coordinates and are 100% identical to themselves.

Which member is real
--------------------
Only one member of a pair can carry the V reading frame: the reverse complement
of a V exon is not a V exon.  So the test is whether the entry, AS STORED (i.e.
in its own declared coding orientation), translates to something with V domain
hallmarks -- the FR2 tryptophan W[VILM]RQ, the [LIV]EW[VILMA] that runs into
CDR2, and the FR3 cysteine Y[YFH]C.

On bAgePho2 IGH that test is unanimous: in all 59 pairs the losing member scores
ZERO motifs while the winner scores 1-3.  It also agrees with the two independent
annotations -- the winner is the one carrying an RSS in 18 of 18 informative
pairs, and the one carrying transcripts in 23 of 24.

What collapsing changes
-----------------------
Denominators, and nothing else.  "25 of 162 IGH V genes carry an RSS" counts 59
loci twice and becomes 25 of 103.  No tract call moves: no donor->parent pair in
either locus is an overlap pair.  Nor could one -- both members point at the same
DNA, so as a conversion donor they are one option, not two.
"""
import argparse
import csv
import re
import sys
from collections import Counter

from Bio.Seq import Seq

from gc_lib import read_fasta

RC = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")

FR2 = re.compile(r"W[VILM]RQ")
LEW = re.compile(r"[LIV]EW[VILMA]")
CYS2 = re.compile(r"Y[YFH]C")


def revcomp(s):
    return s[::-1].translate(RC)


def v_score(nt):
    """(motif hits, -stop codons) over the 3 forward frames of the entry AS STORED.

    Forward frames only, deliberately: the entry is written in its own declared
    coding orientation, so if the annotation is right the reading frame is in
    there. Searching the reverse complement too would just rediscover the
    partner's frame and score both members of every pair identically.
    """
    best = (0, -99)
    for f in range(3):
        end = len(nt) - ((len(nt) - f) % 3)
        aa = str(Seq(nt[f:end]).translate())
        hits = bool(FR2.search(aa)) + bool(LEW.search(aa)) + bool(CYS2.search(aa))
        cand = (hits, -aa.count("*"))
        if cand > best:
            best = cand
    return best


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vgene-fasta", required=True)
    ap.add_argument("--rss-annotation", required=True)
    ap.add_argument("--usage-assignments", required=True)
    ap.add_argument("--locus", required=True)
    ap.add_argument("--out-fasta", required=True)
    ap.add_argument("--out-report", required=True)
    args = ap.parse_args()

    rss = {r["gene"]: r for r in read_tsv(args.rss_annotation)
           if r["locus"] == args.locus}
    seqs = read_fasta(args.vgene_fasta)
    usage = Counter(r["best_gene"] for r in read_tsv(args.usage_assignments)
                    if r["locus"] == args.locus)

    def genomic(n):
        return seqs[n] if rss[n]["strand"] == "+" else revcomp(seqs[n])

    entries = sorted((int(rss[n]["pos"]), int(rss[n]["pos"]) + len(seqs[n]), n)
                     for n in seqs if n in rss)

    # Group into connected components of overlapping intervals. Components rather
    # than pairs because nothing guarantees an overlap set is exactly two -- a
    # three-way pile-up has to collapse to one locus, not to a pair plus a
    # leftover.
    comps = []
    for s, e, n in entries:
        if comps and s < comps[-1][1]:
            comps[-1][1] = max(comps[-1][1], e)
            comps[-1][2].append(n)
        else:
            comps.append([s, e, [n]])

    def rank(n):
        hits, negstops = v_score(seqs[n])
        return (hits, negstops,
                rss[n]["rss_state"] == "rss_present",
                usage.get(n, 0),
                len(seqs[n]), n)

    keep, dropped, rows = [], [], []
    for s, e, members in comps:
        if len(members) == 1:
            keep.append(members[0])
            continue
        ordered = sorted(members, key=rank, reverse=True)
        win, losers = ordered[0], ordered[1:]
        keep.append(win)
        dropped.extend(losers)
        wh, ws = v_score(seqs[win])
        for lo in losers:
            lh, ls = v_score(seqs[lo])
            o = max(int(rss[win]["pos"]), int(rss[lo]["pos"]))
            c = min(int(rss[win]["pos"]) + len(seqs[win]),
                    int(rss[lo]["pos"]) + len(seqs[lo]))
            same = (genomic(win)[o - int(rss[win]["pos"]):c - int(rss[win]["pos"])]
                    == genomic(lo)[o - int(rss[lo]["pos"]):c - int(rss[lo]["pos"])])
            rows.append({
                "kept": win, "kept_strand": rss[win]["strand"],
                "kept_motifs": wh, "kept_rss": rss[win]["rss_state"],
                "kept_transcripts": usage.get(win, 0),
                "dropped": lo, "dropped_strand": rss[lo]["strand"],
                "dropped_motifs": lh, "dropped_rss": rss[lo]["rss_state"],
                "dropped_transcripts": usage.get(lo, 0),
                "overlap_bp": c - o, "identical_over_overlap": same,
                # flag the calls the evidence does not make unanimously
                "unanimous": wh > lh and (rss[lo]["rss_state"] != "rss_present")
                             and usage.get(lo, 0) == 0,
            })

    keep_set = set(keep)
    with open(args.out_fasta, "w") as fh:
        for n in sorted(seqs, key=lambda g: int(rss[g]["pos"]) if g in rss else 0):
            if n in keep_set or n not in rss:
                fh.write(f">{n}\n{seqs[n]}\n")

    with open(args.out_report, "w") as fh:
        cols = list(rows[0].keys()) if rows else ["kept", "dropped"]
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    n_amb = sum(1 for r in rows if not r["unanimous"])
    print(f"{args.locus}: {len(entries)} entries -> {len(keep_set)} kept "
          f"({len(dropped)} dropped as same-DNA duplicates)", file=sys.stderr)
    print(f"  collapsed groups: {sum(1 for c in comps if len(c[2]) > 1)}", file=sys.stderr)
    print(f"  groups where the three tests are NOT unanimous: {n_amb}", file=sys.stderr)
    for r in rows:
        if not r["unanimous"]:
            print(f"    kept {r['kept'].split('.')[1]}({r['kept_strand']}) "
                  f"motifs={r['kept_motifs']} rss={r['kept_rss']} "
                  f"tx={r['kept_transcripts']}  |  dropped "
                  f"{r['dropped'].split('.')[1]}({r['dropped_strand']}) "
                  f"motifs={r['dropped_motifs']} rss={r['dropped_rss']} "
                  f"tx={r['dropped_transcripts']}", file=sys.stderr)
    kept_rss = sum(1 for n in keep_set if rss[n]["rss_state"] == "rss_present")
    print(f"  RSS-bearing among kept: {kept_rss}/{len(keep_set)}", file=sys.stderr)


if __name__ == "__main__":
    main()
