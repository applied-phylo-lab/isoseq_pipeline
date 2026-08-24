"""
Score candidate parents on FR1 alone.

Why restrict to FR1
-------------------
The mosaic-parsimony test (gc_parent_parsimony.py) cannot separate three things:
a wrong RSS restriction, an incomplete RSS annotation, and the fact that a
converted transcript ends up closer to its DONOR than to its parent.  Splitting
that test by whether a tract was detected shows the third effect is real --
converted transcripts pick an RSS parent far less often than unconverted ones.

FR1 sidesteps it rather than trying to price it.  Conversion in these data lands
in CDR1, CDR2 and FR3; FR1 is the conserved 5' end where a transcript still
follows the gene it rearranged from.  Scoring candidate parents on FR1 ONLY
removes the converted sequence from the measurement instead of charging for it.

The cost is power.  FR1 is ~40-80 bp and the candidate genes are ~86% identical,
so only a handful of positions discriminate.  That is the trade being made and
the script reports it: the number of positions that actually differ between the
best RSS and best non-RSS candidate is printed alongside every result, because a
tie on two informative positions is not evidence of anything.

FR1 is taken from the measured framework landmark (see gc_conversion_peaks.py):
FR1 runs from the start of the gene to the beginning of CDR1, which sits 24 bp
before the FR2 start motif.  Each candidate is scored over ITS OWN FR1, which is
the homologous region, and scores are expressed as a mismatch RATE so that genes
with slightly different FR1 lengths remain comparable.
"""
import argparse
import csv
import sys
from collections import defaultdict

import numpy as np

from gc_conversion_peaks import landmarks
from gc_lib import read_fasta, parse_paf, projected_query


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def short(g):
    p = g.split(".")
    return p[1] if len(p) > 1 else g


def scoring_positions(seq, mode, cdr1_len=24, cdr2_len=30):
    """Positions of the gene used to score a candidate parent.

    fr1        the conserved 5' end only: [0, CDR1 start)
    framework  every framework region: FR1 + FR2 + FR3, i.e. the whole gene
               except CDR1, CDR2 and the 3' junction

    `framework` exists because FR1 alone has no power in IGH -- the candidates
    are identical across it, so the median comparison has ZERO informative
    positions. Widening to all framework regions is not an arbitrary rescue: the
    CDR-versus-position test showed conversion concentrates in CDR1/CDR2, so
    their complement is the conversion-poor part of the gene and is the natural
    larger version of the same idea. It is not conversion-FREE -- IGL has a
    tract window at 185-202, inside FR3 -- so `framework` trades a little
    contamination for a lot of power, and both modes are reported.
    """
    lm = dict((lab, x) for x, lab in landmarks(seq)[0])
    a, b, c = lm.get("CDR1|FR2"), lm.get("FR2|CDR2"), lm.get("FR3|CDR3")
    if a is None:
        return None
    fr1 = list(range(0, max(0, a - cdr1_len)))
    if mode == "fr1":
        return fr1
    if mode == "fr1cdr1":
        # FR1 plus CDR1. CDR1 is a hypervariable region, so it is rich in the
        # positions that separate near-identical candidates -- exactly what IGH
        # lacks in FR1 alone. It is only safe to add where CDR1 is not itself a
        # conversion target, and that differs by locus: CDR1 holds 7% of IGH
        # tract bp (one gene) but 23% of IGL's, so this mode helps IGH and
        # contaminates IGL.
        return fr1 + list(range(max(0, a - cdr1_len), a))
    pos = list(fr1)
    if b is not None:
        pos += list(range(a, b))                     # FR2
    if mode == "fr1fr2":
        return sorted(set(pos))
    if b is not None:
        end = c if c is not None else len(seq) - 20  # FR3, stopping before CDR3
        pos += list(range(min(b + cdr2_len, end), end))
    return sorted(set(pos))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--paf", required=True)
    ap.add_argument("--vgene-fasta", required=True)
    ap.add_argument("--rss-annotation", required=True)
    ap.add_argument("--tracts", help="to split converted vs unconverted transcripts")
    ap.add_argument("--clusters",
                    help="gene_clusters.tsv; score FAMILIES rather than genes. "
                         "A family counts as RSS-bearing if ANY member carries "
                         "one, which is the point: near-identical paralogues "
                         "cannot be told apart by 70 transcripts, but if the "
                         "family a transcript matches contains an RSS gene then "
                         "SOME member of it was rearrangeable, and the annotation "
                         "may simply have put the RSS on the wrong copy.")
    ap.add_argument("--locus", required=True)
    ap.add_argument("--mode",
                    choices=("fr1", "fr1cdr1", "fr1fr2", "framework"), default="fr1")
    ap.add_argument("--min-fr1-covered", type=int, default=30)
    ap.add_argument("--perms", type=int, default=5000)
    ap.add_argument("--out-table", required=True)
    ap.add_argument("--out-summary", required=True)
    args = ap.parse_args()
    rng = np.random.default_rng(0)

    seqs = read_fasta(args.vgene_fasta)
    rss = {r["gene"]: r["rss_state"] == "rss_present"
           for r in read_tsv(args.rss_annotation) if r["locus"] == args.locus}
    names = sorted([g for g in seqs if g in rss], key=lambda g: int(g.split(".")[1]))
    ends = {g: scoring_positions(seqs[g], args.mode) for g in names}
    usable = [g for g in names if ends[g] and len(ends[g]) >= args.min_fr1_covered]
    n_rss = sum(rss[g] for g in usable)
    print(f"{args.locus}: {len(usable)}/{len(names)} genes have a locatable FR1 "
          f"({n_rss} with an RSS)", file=sys.stderr)
    print(f"  mode={args.mode}: median {int(np.median([len(ends[g]) for g in usable]))} "
          f"scored bp per gene", file=sys.stderr)

    # optional family collapse: keep, per family, the member each transcript
    # matches best, and label the family by whether ANY member has an RSS
    fam = None
    if args.clusters:
        fam = {}
        for r in read_tsv(args.clusters):
            fam[r["gene"]] = r["cluster"]
        famrss = defaultdict(bool)
        for g, c in fam.items():
            famrss[c] = famrss[c] or rss.get(g, False)
        usable_f = sorted({fam[g] for g in usable if g in fam})
        n_rss_f = sum(famrss[c] for c in usable_f)
        print(f"  families: {len(usable_f)} ({n_rss_f} contain an RSS gene)",
              file=sys.stderr)

    conv = set()
    if args.tracts:
        conv = {t["transcript"] for t in read_tsv(args.tracts)
                if t.get("significant") == "True"}

    # mismatch rate over FR1 for every (transcript, candidate parent)
    score = defaultdict(dict)
    for rec in parse_paf(args.paf):
        if rec.target not in seqs or rec.target not in ends or ends[rec.target] is None:
            continue
        pr = projected_query(rec, seqs[rec.target])
        if pr is None:
            continue
        cov = mm = 0
        for i in ends[rec.target]:
            b = pr[i]
            if b in (None, "-"):
                continue
            cov += 1
            mm += b != seqs[rec.target][i]
        if cov < args.min_fr1_covered:
            continue
        prev = score[rec.query].get(rec.target)
        if prev is None or cov > prev[1]:
            score[rec.query][rec.target] = (mm, cov, mm / cov)

    if fam is not None:
        # collapse each transcript's per-gene scores to per-family bests
        coll = defaultdict(dict)
        for q, per in score.items():
            for g, v in per.items():
                c = fam.get(g)
                if c is None:
                    continue
                if c not in coll[q] or v[2] < coll[q][c][2]:
                    coll[q][c] = v
        score = coll
        rss = dict(famrss)
        usable = usable_f
        n_rss = n_rss_f
        ends = None

    out = open(args.out_table, "w")
    out.write("transcript\tparent\thas_rss\tfr1_mismatches\tfr1_covered\tfr1_rate\n")
    for q, per in score.items():
        for g, (mm, cov, rate) in per.items():
            out.write(f"{q}\t{g}\t{rss[g]}\t{mm}\t{cov}\t{rate:.4f}\n")
    out.close()

    sm = open(args.out_summary, "w")

    def say(s):
        print(s, file=sys.stderr)
        sm.write(s + "\n")

    base = n_rss / len(usable)
    say(f"locus\t{args.locus}\ngenes_with_FR1\t{len(usable)}\trss\t{n_rss}\t"
        f"chance\t{base:.3f}")
    say("")
    say("subset\tn_tx\tmedian_rate_RSS\tmedian_rate_nonRSS\tRSS_cheapest\t"
        "enrichment\tmedian_informative\tp_permuted")

    def analyse(label, sel):
        deltas, wins, tot, rr, nn, infos = [], 0, 0, [], [], []
        keep = []
        for q, per in score.items():
            if not sel(q):
                continue
            a = [(v[2], g) for g, v in per.items() if rss[g]]
            b = [(v[2], g) for g, v in per.items() if not rss[g]]
            if not a or not b:
                continue
            ba, bb = min(a), min(b)
            rr.append(ba[0])
            nn.append(bb[0])
            deltas.append(bb[0] - ba[0])
            # how many FR1 positions actually distinguish the two candidates
            if ends is not None:
                ga, gb = ba[1], bb[1]
                shared = set(ends[ga]) & set(ends[gb])
                infos.append(sum(1 for i in shared if seqs[ga][i] != seqs[gb][i]))
            else:
                infos.append(float("nan"))
            best = min(per.items(), key=lambda kv: kv[1][2])[0]
            wins += rss[best]
            tot += 1
            keep.append(q)
        if not tot:
            return
        obs = wins / tot
        null = []
        for _ in range(args.perms):
            lab = dict(zip(usable, rng.permutation([rss[g] for g in usable])))
            w = t = 0
            for q in keep:
                per = score[q]
                if not per:
                    continue
                best = min(per.items(), key=lambda kv: kv[1][2])[0]
                w += lab.get(best, False)
                t += 1
            if t:
                null.append(w / t)
        p = (np.sum(np.array(null) >= obs) + 1) / (len(null) + 1) if null else float("nan")
        say(f"{label}\t{tot}\t{np.median(rr):.3f}\t{np.median(nn):.3f}\t"
            f"{wins}/{tot} = {obs:.1%}\t{obs / base:.2f}x\t"
            f"{np.median(infos):.0f}\t{p:.4f}"
            if not np.isnan(np.median(infos)) else
            f"{label}\t{tot}\t{np.median(rr):.3f}\t{np.median(nn):.3f}\t"
            f"{wins}/{tot} = {obs:.1%}\t{obs / base:.2f}x\tNA\t{p:.4f}")

    analyse("all", lambda q: True)
    if conv:
        analyse("with_tract", lambda q: q in conv)
        analyse("no_tract", lambda q: q not in conv)

    # per-gene: how often each gene is the best FR1 match
    say("")
    say("# per-gene, best FR1 match")
    say("gene\thas_rss\tn_tx_best\tmedian_fr1_rate")
    win = defaultdict(list)
    for q, per in score.items():
        if per:
            g = min(per.items(), key=lambda kv: kv[1][2])[0]
            win[g].append(per[g][2])
    for g in sorted(win, key=lambda g: -len(win[g]))[:15]:
        say(f"{short(g) if ends is not None else 'cluster' + str(g)}\t"
            f"{rss[g]}\t{len(win[g])}\t{np.median(win[g]):.3f}")
    sm.close()


if __name__ == "__main__":
    main()
