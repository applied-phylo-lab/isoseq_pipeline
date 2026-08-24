"""
Summary of the evidence that only RSS-bearing V genes act as parents.

What the three panels answer
----------------------------
A  Does restricting the scored region help? Enrichment of "the best-matching
   candidate carries an RSS" over the chance rate, for each scored region.
   Widening the window from FR1 to the full framework does NOT help: it
   reintroduces the converted sequence the restriction was meant to exclude.

B  Is the signal an artefact of conversion? The same enrichment split by whether
   the transcript carries a detected tract. If conversion drives the result, the
   two bars diverge; if the scored region is genuinely conversion-poor, they
   agree. This is the panel that validates the choice of region.

C  Which genes are actually being picked, and do they carry an RSS? Genes whose
   FR1 the transcripts match best. A non-RSS gene high in this ranking with a
   near-zero mismatch rate is a missing-RSS candidate, not a counterexample.
"""
import argparse
import csv
import sys
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from gc_palette import save_figure, LOCUS, SEGMENT, GREY, GREY_DARK, INK, NO, YES


def parse_summary(path):
    """(chance, {subset: (n, rss_frac, enrichment, informative)})"""
    chance, rows = None, {}
    with open(path) as fh:
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if f and f[0] == "genes_with_FR1":
                chance = float(f[-1])
            if len(f) >= 8 and f[0] in ("all", "with_tract", "no_tract"):
                frac = float(f[4].split("=")[1].strip().rstrip("%")) / 100
                try:
                    info = float(f[6])
                except ValueError:
                    info = float("nan")     # family level has no gene coordinates
                rows[f[0]] = (int(f[1]), frac, float(f[5].rstrip("x")), info,
                              float(f[7]) if len(f) > 7 else float("nan"))
    return chance, rows


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--summary", nargs="+", required=True,
                    metavar="LOCUS:REGION:path")
    ap.add_argument("--fr1-table", nargs="+", default=[], metavar="LOCUS:path")
    ap.add_argument("--clusters", nargs="+", default=[], metavar="LOCUS:path")
    ap.add_argument("--rss-annotation", required=True)
    ap.add_argument("--out-figure", required=True)
    args = ap.parse_args()

    data = defaultdict(dict)
    chance = {}
    for spec in args.summary:
        loc, region, path = spec.split(":", 2)
        ch, rows = parse_summary(path)
        chance[loc] = ch
        data[loc][region] = rows

    fig, axes = plt.subplots(2, 2, figsize=(14.0, 9.4))
    axes = axes.ravel()
    regions = ["FR1", "FR1+CDR1", "FR1+FR2", "framework", "family"]
    loci = [l for l in ("IGH", "IGL") if l in data]

    # ── A: enrichment by scored region ────────────────────────────────────────
    ax = axes[0]
    w = 0.36
    for k, loc in enumerate(loci):
        vals = [data[loc].get(r, {}).get("all", (0, 0, np.nan))[2] for r in regions]
        x = np.arange(len(regions)) + (k - (len(loci) - 1) / 2) * w
        ax.bar(x, vals, w, color=LOCUS.get(loc, INK), edgecolor="black", lw=0.6,
               label=loc)
        for xi, v in zip(x, vals):
            if not np.isnan(v):
                ax.text(xi, v + 0.3, f"{v:.1f}×", ha="center", va="bottom",
                        fontsize=7.5, fontweight="bold", color=LOCUS.get(loc, INK))
    ax.axhline(1.0, color=INK, lw=1.2, ls="--")
    ax.annotate("chance", (len(regions) - 0.45, 1.0), xytext=(0, 3),
                textcoords="offset points", fontsize=7.5, color=INK, ha="right")
    ax.set_xticks(range(len(regions)))
    ax.set_xticklabels(regions, fontsize=8.5)
    ax.set_ylabel("enrichment: best match carries an RSS", fontsize=9)
    ax.set_title(r"$\bf{A}$  IGL wants the narrowest region; IGH wants families",
                 fontsize=10.5, loc="left", x=0.0)
    ax.legend(fontsize=8.5, frameon=False)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

    # ── B: conversion control ─────────────────────────────────────────────────
    ax = axes[1]
    labels, vals, cols, hat = [], [], [], []
    for loc in loci:
        best = "FR1" if loc == "IGL" else "family"
        for reg, h in ((best, ""), ("framework", "//")):
            for sub, tag in (("with_tract", "conv"), ("no_tract", "unconv")):
                r = data[loc].get(reg, {}).get(sub)
                if not r:
                    continue
                labels.append(f"{loc}\n{reg}\n{tag}")
                vals.append(r[2]); cols.append(LOCUS.get(loc, INK)); hat.append(h)
    xs = np.arange(len(vals))
    for b, h in zip(ax.bar(xs, vals, 0.75, color=cols, edgecolor="black", lw=0.6), hat):
        b.set_hatch(h)
    ax.axhline(1.0, color=INK, lw=1.2, ls="--")
    ax.set_xticks(xs); ax.set_xticklabels(labels, fontsize=6.5)
    ax.set_ylabel("enrichment", fontsize=9)
    ax.set_title(r"$\bf{B}$  in the chosen region, converted and unconverted agree",
                 fontsize=10.5, loc="left", x=0.0)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

    # ── C: which genes the FR1 evidence points to ─────────────────────────────
    ax = axes[2]
    rss = {}
    with open(args.rss_annotation) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            rss[r["gene"]] = r["rss_state"] == "rss_present"
    picked = []
    for spec in args.fr1_table:
        loc, path = spec.split(":", 1)
        best = {}
        with open(path) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                q, g, rate = r["transcript"], r["parent"], float(r["fr1_rate"])
                if q not in best or rate < best[q][1]:
                    best[q] = (g, rate)
        cnt = defaultdict(list)
        for g, rate in best.values():
            cnt[g].append(rate)
        for g, rates in cnt.items():
            picked.append((loc, g, len(rates), float(np.median(rates)), rss.get(g, False)))
    picked.sort(key=lambda t: -t[2])
    picked = picked[:12]
    y = np.arange(len(picked)); mx = max(p[2] for p in picked)
    ax.barh(y, [p[2] for p in picked],
            color=[LOCUS.get(p[0], INK) for p in picked], edgecolor="black", lw=0.5)
    for i, p in enumerate(picked):
        if p[4]:
            ax.plot([-mx * 0.035], [i], marker="o", ms=6, color=SEGMENT["V"],
                    markeredgecolor="black", markeredgewidth=0.5, clip_on=False)
        ax.text(p[2] + mx * 0.02, i, f" {p[3]:.3f}", va="center", fontsize=7,
                color=GREY_DARK)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{p[0]} {p[1].split('.')[1]}" for p in picked], fontsize=7.5)
    ax.invert_yaxis()
    ax.set_xlabel("transcripts whose FR1 matches this gene best", fontsize=9)
    ax.set_title(r"$\bf{C}$  genes the FR1 evidence points to"
                 "\n(• has an RSS; number = median mismatch rate)",
                 fontsize=10.5, loc="left", x=0.0)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

    # ── D: mixed families ─────────────────────────────────────────────────────
    ax = axes[3]
    mixed = []
    for spec in args.clusters:
        loc, path = spec.split(":", 1)
        cl = defaultdict(list)
        for r in csv.DictReader(open(path), delimiter="\t"):
            cl[r["cluster"]].append(r)
        for c, v in cl.items():
            r_ = [x for x in v if x["rss_state"] == "rss_present"]
            n_ = [x for x in v if x["rss_state"] != "rss_present"]
            if r_ and n_:
                mixed.append((loc, c, len(r_), len(n_),
                              sum(int(x["n_transcripts"]) for x in n_)))
    mixed.sort(key=lambda t: -(t[2] + t[3]))
    mixed = mixed[:10]
    y = np.arange(len(mixed))
    ax.barh(y, [m[2] for m in mixed], color=SEGMENT["V"], edgecolor="black",
            lw=0.5, label="members WITH an RSS")
    ax.barh(y, [m[3] for m in mixed], left=[m[2] for m in mixed], color=GREY,
            edgecolor="black", lw=0.5, label="members without")
    for i, m in enumerate(mixed):
        if m[4]:
            ax.text(m[2] + m[3] + 0.12, i, f" {m[4]} tx on non-RSS members",
                    va="center", fontsize=6.8, color=NO)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{m[0]} fam {m[1]}" for m in mixed], fontsize=7.5)
    ax.invert_yaxis()
    ax.set_xlabel("genes in the family", fontsize=9)
    ax.set_title(r"$\bf{D}$  mixed families: near-identical paralogues where"
                 "\nonly some copies carry an annotated RSS",
                 fontsize=10.5, loc="left", x=0.0)
    ax.legend(fontsize=7.5, frameon=False, loc="lower right")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

    fig.tight_layout()
    save_figure(fig, args.out_figure)
    print("wrote", args.out_figure, file=sys.stderr)


if __name__ == "__main__":
    main()
