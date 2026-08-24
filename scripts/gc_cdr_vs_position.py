"""
Are conversion tracts in CDRs, or just in the middle of the gene?

The confound
------------
CDR1 and CDR2 sit in the middle of a V exon, and conversion tracts are also in
the middle -- nothing is called in the first ~40 bp or the last ~60 bp.  So "most
tracts overlap a CDR" is true under either hypothesis and cannot distinguish
them.  Three nulls are used, each breaking a different thing:

  UNIFORM     tracts placed anywhere in the assessable region. Tests only whether
              tracts are non-random at all; it does NOT control for position, so
              a significant result here means little on its own.

  POSITIONAL  tracts placed only within the positional range where tracts are
              actually observed, pooled across genes. This holds "tracts are in
              the middle" fixed and asks whether CDRs are still favoured inside
              that middle. This is the test that matters.

  SWAP        each tract keeps its exact coordinates but is reassigned to a
              different parent gene. Absolute position is preserved perfectly;
              only the alignment with THIS gene's CDR is broken. The parents'
              CDR2 boundaries range over ~40 bp, so this is a real shuffle.

A tract is scored by the fraction of its positions that fall inside a CDR, so a
tract straddling a boundary contributes partially rather than being forced into
one bin.

Caveat carried in the output
----------------------------
Only the framework BOUNDARIES are measured from protein motifs; CDR extents are
canonical lengths (CDR1 24 bp before FR2, CDR2 30 bp after it). The result is
therefore reported across a range of assumed CDR2 lengths, because a test whose
answer depends on an assumed constant should show that dependence.
"""
import argparse
import csv
import random
import statistics as st
import sys
from collections import defaultdict

from gc_conversion_peaks import landmarks
from gc_lib import read_fasta


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def short(g):
    p = g.split(".")
    return p[1] if len(p) > 1 else g


def cdr_positions(seq, cdr1_len, cdr2_len):
    """Set of positions inside CDR1 or CDR2, from the measured FR2 boundaries."""
    lm = dict((lab, x) for x, lab in landmarks(seq)[0])
    out = set()
    a, b = lm.get("CDR1|FR2"), lm.get("FR2|CDR2")
    if a is not None:
        out |= set(range(max(0, a - cdr1_len), a))
    if b is not None:
        out |= set(range(b, min(len(seq), b + cdr2_len)))
    return out


def frac_in_cdr(spans, cdr):
    tot = hit = 0
    for s, e in spans:
        for i in range(s, e + 1):
            tot += 1
            hit += i in cdr
    return hit / tot if tot else 0.0


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--locus", action="append", required=True,
                    metavar="NAME=fasta,tracts")
    ap.add_argument("--perms", type=int, default=20000)
    ap.add_argument("--cdr1-len", type=int, default=24)
    ap.add_argument("--out-table", required=True)
    args = ap.parse_args()
    random.seed(0)

    rows = []
    for spec in args.locus:
        name, rest = spec.split("=", 1)
        fa, tr = rest.split(",")
        seqs = read_fasta(fa)
        tracts = [t for t in read_tsv(tr) if t.get("significant") == "True"]
        # DISTINCT tracts only: 123 IGL calls collapse to 21 real events, and
        # counting calls would let one expanded clone carry the whole test.
        uniq = sorted({(t["parent"], int(t["start"]), int(t["end"])) for t in tracts})
        byp = defaultdict(list)
        for p, s, e in uniq:
            byp[p].append((s, e))
        parents = sorted(byp)
        print(f"\n{'='*70}\n{name}: {len(uniq)} distinct tracts over {len(parents)} parents",
              file=sys.stderr)

        for cdr2_len in (21, 30, 39):
            cdr = {p: cdr_positions(seqs[p], args.cdr1_len, cdr2_len) for p in parents}
            lim = {p: len(seqs[p]) - 20 for p in parents}
            # scored per parent, so every position is judged against ITS OWN
            # gene's CDRs rather than a pooled coordinate set
            tot = hit = 0
            for p in parents:
                for s, e in byp[p]:
                    for i in range(s, e + 1):
                        tot += 1
                        hit += i in cdr[p]
            obs = hit / tot

            # what fraction of each gene, and of the tract-bearing middle, is CDR
            cdr_share_gene = st.mean(len(cdr[p]) / lim[p] for p in parents)
            lo = min(s for p in parents for s, _e in byp[p])
            hi = max(e for p in parents for _s, e in byp[p])
            cdr_share_mid = st.mean(
                len({i for i in cdr[p] if lo <= i <= hi}) / max(1, min(hi, lim[p]) - lo + 1)
                for p in parents)

            nulls = {}
            # UNIFORM
            cnt = 0
            for _ in range(args.perms):
                t = h = 0
                for p in parents:
                    for s, e in byp[p]:
                        w = e - s + 1
                        st0 = random.randint(0, max(0, lim[p] - w))
                        for i in range(st0, st0 + w):
                            t += 1
                            h += i in cdr[p]
                cnt += (h / t) >= obs
            nulls["uniform"] = (cnt + 1) / (args.perms + 1)

            # POSITIONAL: same pooled positional range as observed
            cnt = 0
            for _ in range(args.perms):
                t = h = 0
                for p in parents:
                    for s, e in byp[p]:
                        w = e - s + 1
                        top = min(hi, lim[p] - 1) - w
                        st0 = random.randint(lo, max(lo, top))
                        for i in range(st0, st0 + w):
                            t += 1
                            h += i in cdr[p]
                cnt += (h / t) >= obs
            nulls["positional"] = (cnt + 1) / (args.perms + 1)

            # SWAP: identical coordinates, different gene
            if len(parents) > 1:
                cnt = 0
                for _ in range(args.perms):
                    t = h = 0
                    for p in parents:
                        for s, e in byp[p]:
                            q = random.choice([x for x in parents if x != p])
                            if e >= lim[q]:
                                q = p
                            for i in range(s, e + 1):
                                t += 1
                                h += i in cdr[q]
                    cnt += (h / t) >= obs
                nulls["swap"] = (cnt + 1) / (args.perms + 1)
            else:
                nulls["swap"] = float("nan")

            print(f"  CDR2 assumed {cdr2_len} bp: observed {obs:.1%} of tract "
                  f"positions in a CDR", file=sys.stderr)
            print(f"     CDR is {cdr_share_gene:.1%} of the assessable gene, "
                  f"{cdr_share_mid:.1%} of the tract-bearing middle ({lo}-{hi})",
                  file=sys.stderr)
            for k, v in nulls.items():
                print(f"     p({k:<10}) = {v:.4f}", file=sys.stderr)
            rows.append((name, cdr2_len, len(uniq), obs, cdr_share_gene,
                         cdr_share_mid, nulls["uniform"], nulls["positional"],
                         nulls["swap"]))

    with open(args.out_table, "w") as fh:
        fh.write("locus\tcdr2_len_assumed\tn_distinct_tracts\tobserved_frac_in_cdr\t"
                 "cdr_frac_of_gene\tcdr_frac_of_tract_range\tp_uniform\t"
                 "p_positional\tp_swap\n")
        for r in rows:
            fh.write("\t".join(f"{x:.4f}" if isinstance(x, float) else str(x)
                               for x in r) + "\n")


if __name__ == "__main__":
    main()
