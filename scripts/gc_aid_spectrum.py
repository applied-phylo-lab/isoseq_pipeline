"""
Separate somatic hypermutation from gene conversion using the mutation spectrum
alone -- no donor database required.

The logic
---------
AID initiates both processes, but only SHM leaves AID's fingerprint on the
result:

  * SHM         the mutation happens AT the AID-targeted cytosine (or its
                error-prone repair).  Mutations therefore pile up on the
                WRCY / RGYW hotspot motifs, favour transitions (C->T, G->A),
                and concentrate on C:G pairs.
  * conversion  AID makes the same lesion, but repair copies a DONOR.  The
                differences that result sit wherever the donor happens to
                differ from the parent, which has nothing to do with where AID
                bound.  Tract-internal differences should therefore show no
                hotspot enrichment, no transition bias, and no C:G preference.

So a called conversion tract carrying an AID-like spectrum is clustered SHM
wearing a disguise, and this can be checked without knowing, or trusting, which
donor was assigned.  That independence is the point: it still works when the
germline annotation is incomplete.

Motifs (the mutated base is capitalised)
    WRCY   [AT][AG] C [CT]        W=A/T  R=A/G  Y=C/T
    RGYW   [AG] G [CT][AT]        (reverse complement of WRCY)
    coldspot  SYC / GRS           S=G/C

The null
--------
Hotspot fraction depends on the base composition of the V gene, so raw
percentages mean nothing on their own.  For each transcript its observed
mutations are redistributed at random over the positions actually covered in
that same gene, keeping the mutation count fixed, and the hotspot fraction is
recomputed.  Repeating that gives the fraction expected with no AID targeting,
and an empirical p-value.
"""
import argparse
import random
import statistics
import sys
from collections import Counter, defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np

from gc_lib import read_fasta, parse_paf, projected_query

from gc_palette import save_figure, CLASS_A as OUTSIDE_C, CLASS_B as INSIDE_C, GREY_DARK as NULL_C

PURINES = set("AG")


def is_wrcy(seq, i):
    """Mutated C at i, preceded by W R and followed by Y."""
    if i < 2 or i + 1 >= len(seq) or seq[i] != "C":
        return False
    return seq[i - 2] in "AT" and seq[i - 1] in "AG" and seq[i + 1] in "CT"


def is_rgyw(seq, i):
    """Mutated G at i, preceded by R and followed by Y W."""
    if i < 1 or i + 2 >= len(seq) or seq[i] != "G":
        return False
    return seq[i - 1] in "AG" and seq[i + 1] in "CT" and seq[i + 2] in "AT"


def is_hotspot(seq, i):
    return is_wrcy(seq, i) or is_rgyw(seq, i)


def is_coldspot(seq, i):
    """SYC (mutated C) or its complement GRS (mutated G)."""
    if i >= 2 and seq[i] == "C" and seq[i - 2] in "GC" and seq[i - 1] in "CT":
        return True
    if i + 2 < len(seq) and seq[i] == "G" and seq[i + 1] in "AG" and seq[i + 2] in "GC":
        return True
    return False


def transition(ref, alt):
    return {ref, alt} in ({"A", "G"}, {"C", "T"})


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--paf", required=True)
    ap.add_argument("--vgene-fasta", required=True)
    ap.add_argument("--assignments", required=True)
    ap.add_argument("--tracts", help="tracts TSV; without it every difference is "
                                     "pooled and only the overall spectrum is reported")
    ap.add_argument("--locus", required=True)
    ap.add_argument("--junction-margin", type=int, default=20,
                    help="ignore this many 3' bases, where V(D)J junction "
                         "formation creates differences that are not mutations")
    ap.add_argument("--permutations", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out-figure", required=True)
    ap.add_argument("--out-stats", required=True)
    args = ap.parse_args()

    rng = random.Random(args.seed)
    vgenes = read_fasta(args.vgene_fasta)

    parent_of = {}
    with open(args.assignments) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            r = dict(zip(hdr, line.rstrip("\n").split("\t")))
            if r["locus"] == args.locus and r["best_gene"] in vgenes:
                parent_of[r["transcript"]] = r["best_gene"]

    # tract intervals per transcript, in parent coordinates
    tracts = defaultdict(list)
    if args.tracts:
        with open(args.tracts) as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            for line in fh:
                r = dict(zip(hdr, line.rstrip("\n").split("\t")))
                if r.get("significant") == "True":
                    tracts[r["transcript"]].append((int(r["start"]), int(r["end"])))

    # best projection per transcript against its own parent
    proj_of = {}
    for rec in parse_paf(args.paf):
        p = parent_of.get(rec.query)
        if p is None or rec.target != p:
            continue
        pr = projected_query(rec, vgenes[p])
        if pr is None:
            continue
        prev = proj_of.get(rec.query)
        if prev is None or sum(x is not None for x in pr) > sum(x is not None for x in prev):
            proj_of[rec.query] = pr

    muts = {"inside": [], "outside": []}       # (gene, pos, ref, alt)
    covered_by_gene = defaultdict(set)
    obs_counts = {"inside": defaultdict(int), "outside": defaultdict(int)}

    for tid, proj in proj_of.items():
        gene = parent_of[tid]
        gseq = vgenes[gene]
        limit = max(0, len(gseq) - args.junction_margin)
        ivs = tracts.get(tid, [])
        for i in range(limit):
            b = proj[i]
            if b is None or b == "-" or b not in "ACGT":
                continue
            covered_by_gene[gene].add(i)
            if b == gseq[i]:
                continue
            cls = "inside" if any(s <= i <= e for s, e in ivs) else "outside"
            muts[cls].append((gene, i, gseq[i], b))
            obs_counts[cls][tid] += 1

    def spectrum(rows):
        if not rows:
            return None
        hot = sum(1 for g, i, r, a in rows if is_hotspot(vgenes[g], i))
        cold = sum(1 for g, i, r, a in rows if is_coldspot(vgenes[g], i))
        ts = sum(1 for g, i, r, a in rows if transition(r, a))
        cg = sum(1 for g, i, r, a in rows if r in "CG")
        n = len(rows)
        return {"n": n, "hotspot": hot / n, "coldspot": cold / n,
                "ti_tv": ts / max(1, n - ts), "cg_frac": cg / n}

    def null_motifs(rows, k):
        """
        Redistribute the same number of mutations over covered positions.

        Returns hotspot AND coldspot fractions from the same draws, so the two
        are directly comparable: the coldspot arm is what makes the test
        two-sided in substance rather than just in wording.  AID does not merely
        prefer WRCY/RGYW, it actively avoids SYC/GRS, so a genuine AID signature
        must show enrichment in one and depletion in the other.  A process that
        simply concentrates mutations anywhere would move both the same way.
        """
        by_gene = Counter(g for g, *_ in rows)
        hots, colds, cgs = [], [], []
        for _ in range(k):
            hot = cold = cg = tot = 0
            for gene, cnt in by_gene.items():
                pool = sorted(covered_by_gene[gene])
                if not pool:
                    continue
                pick = rng.sample(pool, min(cnt, len(pool)))
                hot += sum(1 for i in pick if is_hotspot(vgenes[gene], i))
                cold += sum(1 for i in pick if is_coldspot(vgenes[gene], i))
                cg += sum(1 for i in pick if vgenes[gene][i] in "CG")
                tot += len(pick)
            if tot:
                hots.append(hot / tot)
                colds.append(cold / tot)
                cgs.append(cg / tot)
        return hots, colds, cgs

    results, nulls, nulls_cold, nulls_cg = {}, {}, {}, {}
    for cls in ("outside", "inside"):
        s = spectrum(muts[cls])
        if s is None:
            continue
        results[cls] = s
        (nulls[cls], nulls_cold[cls],
         nulls_cg[cls]) = null_motifs(muts[cls], args.permutations)

    with open(args.out_stats, "w") as fh:
        fh.write("class\tmetric\tvalue\n")
        for cls, s in results.items():
            for motif, nl in (("hotspot", nulls[cls]), ("coldspot", nulls_cold[cls])):
                exp = statistics.mean(nl) if nl else float("nan")
                obs = s[motif]
                p_up = p_dn = p_two = float("nan")
                if nl:
                    n = len(nl)
                    p_up = (sum(1 for x in nl if x >= obs) + 1) / (n + 1)
                    p_dn = (sum(1 for x in nl if x <= obs) + 1) / (n + 1)
                    p_two = min(1.0, 2 * min(p_up, p_dn))
                if motif == "hotspot":
                    fh.write(f"{cls}\tn_differences\t{s['n']}\n")
                fh.write(f"{cls}\t{motif}_fraction\t{obs:.4f}\n")
                fh.write(f"{cls}\t{motif}_expected\t{exp:.4f}\n")
                fh.write(f"{cls}\t{motif}_enrichment\t{obs/exp:.3f}\n")
                # TWO-SIDED. A one-sided test in a fixed direction reports p=1
                # when the effect is strongly in the other direction, which reads
                # as "no effect" and is exactly backwards. The direction is
                # carried by the enrichment value, so the p-value only has to say
                # whether the departure from the null is real.
                fh.write(f"{cls}\t{motif}_p_value\t{p_two:.4g}\n")
                fh.write(f"{cls}\t{motif}_p_enriched_onesided\t{p_up:.4g}\n")
                fh.write(f"{cls}\t{motif}_p_depleted_onesided\t{p_dn:.4g}\n")
            fh.write(f"{cls}\ttransition_transversion\t{s['ti_tv']:.3f}\n")
            ng = nulls_cg[cls]
            exp_cg = statistics.mean(ng) if ng else float("nan")
            fh.write(f"{cls}\tCG_targeted_fraction\t{s['cg_frac']:.4f}\n")
            fh.write(f"{cls}\tCG_expected\t{exp_cg:.4f}\n")
            fh.write(f"{cls}\tCG_enrichment\t{s['cg_frac']/exp_cg:.3f}\n")

    # ── figure ───────────────────────────────────────────────────────────────
    have_both = len(results) == 2
    fig, axes = plt.subplots(1, 3 if not have_both else 4,
                             figsize=(15 if have_both else 11.5, 4.0))
    order = [c for c in ("outside", "inside") if c in results]
    cols = {"outside": OUTSIDE_C, "inside": INSIDE_C}
    lbl = {"outside": "outside tracts\n(expected: SHM)",
           "inside": "inside tracts\n(expected: conversion)"}

    # Panel A shows all FOUR cases -- {hotspot, coldspot} x {outside, inside} --
    # as fold-change against the permutation null, because that is what makes the
    # hypothesis falsifiable.  The prediction is directional and specific:
    #
    #   outside tracts = SHM        hotspot UP    coldspot DOWN
    #   inside tracts  = conversion hotspot ~1    coldspot ~1
    #
    # Plotting only hotspots (the old panel) could not distinguish "AID acted
    # here" from "mutations are clustered here for any reason at all". The
    # coldspot arm is the discriminator: only AID pushes the two in opposite
    # directions. Drawn on a log axis so enrichment and depletion are
    # symmetrical about the no-effect line rather than one being visually
    # compressed into it.
    ax = axes[0]
    cases = [(c, m) for c in order for m in ("hotspot", "coldspot")]
    xc = np.arange(len(cases))
    vals, ps, colours = [], [], []
    for c, m in cases:
        nl = nulls[c] if m == "hotspot" else nulls_cold[c]
        exp = statistics.mean(nl) if nl else float("nan")
        obs = results[c][m]
        vals.append(obs / exp if exp else float("nan"))
        # two-sided: the bar already shows which way the effect goes
        n = len(nl)
        up = (sum(1 for v in nl if v >= obs) + 1) / (n + 1)
        dn = (sum(1 for v in nl if v <= obs) + 1) / (n + 1)
        ps.append(min(1.0, 2 * min(up, dn)))
        colours.append(cols[c])
    bars = ax.bar(xc, vals, 0.62, color=colours, edgecolor="black", lw=.7)
    for bar, (_, m) in zip(bars, cases):
        if m == "coldspot":
            bar.set_hatch("//")
    ax.axhline(1.0, color="black", lw=1.2)
    ax.set_yscale("log")
    ax.set_yticks([0.25, 0.5, 1, 2, 4])
    ax.set_yticklabels(["0.25×", "0.5×", "1× (no effect)", "2×", "4×"], fontsize=8)
    for xi, (v, p) in enumerate(zip(vals, ps)):
        # Uniformly just above the bar top: for a depleted bar that lands in
        # the empty band under the 1x line, which keeps it off the bar itself.
        arrow = "↑" if v > 1 else "↓"
        ax.text(xi, v * 1.07, f"{arrow}×{v:.2f}\np={p:.3g}", ha="center",
                va="bottom", fontsize=7.2, fontweight="bold")
    ax.set_xticks(xc)
    ax.set_xticklabels([f"{m}\n{'outside' if c == 'outside' else 'inside'}"
                        for c, m in cases], fontsize=8)
    lo, hi = min(vals), max(vals)
    ax.set_ylim(min(lo * 0.45, 0.3), max(hi * 2.6, 3.0))
    ax.set_ylabel("observed / expected", fontsize=9)
    handles = [Patch(facecolor=cols[c], edgecolor="black", label=lbl[c].split("\n")[0])
               for c in order]
    handles += [Patch(facecolor="white", edgecolor="black", label="hotspot WRCY/RGYW"),
                Patch(facecolor="white", edgecolor="black", hatch="//",
                      label="coldspot SYC/GRS")]
    ax.legend(handles=handles, fontsize=6.8, loc="upper right", ncol=2)
    ax.set_title("A  AID signature — all four cases\n"
                 "p two-sided; arrow gives the direction",
                 fontsize=9.5, fontweight="bold", loc="left")

    x = np.arange(len(order))

    ax = axes[1]
    ax.bar(x, [results[c]["ti_tv"] for c in order], .5,
           color=[cols[c] for c in order], edgecolor="black", lw=.6)
    ax.axhline(0.5, color="black", ls="--", lw=1)
    ax.text(len(order) - .5, .52, "no bias (0.5)", fontsize=7.5, ha="right")
    ax.set_xticks(x)
    ax.set_xticklabels([lbl[c] for c in order], fontsize=8)
    ax.set_ylabel("transitions / transversions")
    ax.set_title("B  transition bias\nconfounded inside tracts — see caption",
                 fontsize=9.5, fontweight="bold", loc="left")

    # The baseline here is NOT 0.5.  It is the C+G content of the positions
    # actually covered in these genes, which is what a process with no C:G
    # preference would hit -- taken from the same permutation as panel A.
    ax = axes[2]
    ax.bar(x - .2, [results[c]["cg_frac"] for c in order], .4,
           color=[cols[c] for c in order], edgecolor="black", lw=.6,
           label="observed")
    ax.bar(x + .2, [statistics.mean(nulls_cg[c]) for c in order], .4,
           color=NULL_C, edgecolor="black", lw=.6, label="expected (composition)")
    for i, c in enumerate(order):
        e = statistics.mean(nulls_cg[c])
        ax.text(i, max(results[c]["cg_frac"], e) * 1.03,
                f"×{results[c]['cg_frac']/e:.2f}", ha="center", fontsize=8,
                fontweight="bold")
    ax.set_ylim(0, max(max(results[c]["cg_frac"] for c in order),
                       max(statistics.mean(nulls_cg[c]) for c in order)) * 1.30)
    ax.set_xticks(x)
    ax.set_xticklabels([lbl[c] for c in order], fontsize=8)
    ax.set_ylabel("fraction at a C or G")
    ax.legend(fontsize=7, loc="upper right")
    ax.set_title("C  C:G targeting\nbaseline = actual C+G content, not 0.5",
                 fontsize=9.5, fontweight="bold", loc="left")

    if have_both:
        ax = axes[3]
        types = ["C>T", "G>A", "A>G", "T>C", "C>A", "C>G", "G>C", "G>T",
                 "A>C", "A>T", "T>A", "T>G"]
        for k, c in enumerate(order):
            cnt = Counter(f"{r}>{a}" for g, i, r, a in muts[c])
            tot = sum(cnt.values())
            ax.bar(np.arange(len(types)) + (k - .5) * .4,
                   [cnt.get(t, 0) / tot for t in types], .4,
                   color=cols[c], edgecolor="black", lw=.3, label=lbl[c])
        ax.set_xticks(range(len(types)))
        ax.set_xticklabels(types, rotation=90, fontsize=7)
        ax.set_ylabel("fraction of differences")
        ax.legend(fontsize=7)
        ax.set_title("D  substitution spectrum", fontsize=10,
                     fontweight="bold", loc="left")

    for a in axes:
        for s_ in ("top", "right"):
            a.spines[s_].set_visible(False)
    fig.suptitle(f"{args.locus} — AID mutation-spectrum test "
                 f"(donor-database independent)", fontsize=12,
                 fontweight="bold", y=1.05)
    fig.tight_layout()
    save_figure(fig, args.out_figure)

    for cls in order:
        s, nl = results[cls], nulls[cls]
        exp = statistics.mean(nl)
        p = (sum(1 for v in nl if v >= s["hotspot"]) + 1) / (len(nl) + 1)
        print(f"{args.locus} {cls:8s}: n={s['n']:6d}  hotspot={s['hotspot']:.3f} "
              f"(exp {exp:.3f}, ×{s['hotspot']/exp:.2f}, p={p:.3g})  "
              f"Ti/Tv={s['ti_tv']:.2f}  C:G={s['cg_frac']:.3f}", file=sys.stderr)


if __name__ == "__main__":
    main()
