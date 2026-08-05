"""
One-page overview of an IsoSeq IG run.

Panels
------
  A  read/transcript funnel, from raw HiFi reads down to IG transcripts
  B  IGH vs IGL split
  C  identity of each transcript to its assigned germline V gene
  D  top-N IGH V gene usage
  E  top-N IGL V gene usage
  F  minimap2 filter vs immunotools: what each method recovers, and why
     immunotools drops what it drops
  G  IGH locus map
  H  IGL locus map

Locus map encoding
------------------
  tick direction : + strand up, - strand down
  tick height    : gene has an RSS (tall) or not (short)
  tick colour    : transcripts calling that gene their parent

A note on the IG percentage
---------------------------
When the pipeline screens FLNC reads for IG before clustering, the transcripts
reaching the clustering step are already IG-enriched, so
"IG transcripts / total transcripts" approaches 1 and says nothing about the
library.  The meaningful figure is IG reads as a fraction of all FLNC reads.
Panel A shows the whole funnel so the two cannot be confused.
"""
import argparse
from collections import Counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap, PowerNorm
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, Patch

from gc_palette import save_figure, LOCUS as LOCUS_C, GREY, GREY_DARK, SEGMENT, ramp

ZERO_C = GREY             # gene with no transcripts
FUNNEL_C = GREY_DARK      # neutral: the funnel is not a locus


def locus_cmap(locus):
    """Sequential ramp for a locus; see gc_palette.ramp for why it is not white-based."""
    return ramp(LOCUS_C[locus])


def make_norm(vmax):
    """Compressed so a locus with a much lower maximum stays legible."""
    return PowerNorm(gamma=0.5, vmin=0, vmax=max(1, vmax))


def read_tsv(path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(header, line.rstrip("\n").split("\t")))
                for line in fh if line.strip()]


def read_kv(path):
    d = {}
    for r in read_tsv(path):
        k = list(r)
        d[r[k[0]]] = r[k[1]]
    return d


def short(name):
    try:
        return name.split(".")[1]
    except IndexError:
        return name


def panel_funnel(ax, steps):
    ax.axis("off")
    n0 = steps[0][1]
    y = len(steps)
    for i, (label, n, note) in enumerate(steps):
        frac = n / n0 if n0 else 0
        w = max(0.035, frac ** 0.28)
        shade = 1 - 0.55 * (1 - i / max(1, len(steps) - 1))
        ax.add_patch(FancyBboxPatch((0.5 - w / 2, y - i - 0.78), w, 0.6,
                                    boxstyle="round,pad=0.006",
                                    facecolor=mcolors.to_rgba(FUNNEL_C, shade),
                                    edgecolor="black", lw=.7))
        txt = f"{n:,}" + (f"   ({note})" if note else "")
        if w > 0.42:
            ax.text(0.5, y - i - 0.40, label, ha="center", va="center",
                    fontsize=9.5, fontweight="bold")
            ax.text(0.5, y - i - 0.62, txt, ha="center", va="center", fontsize=8.5)
        else:
            ax.text(0.5 + w / 2 + 0.03, y - i - 0.40, label, ha="left",
                    va="center", fontsize=9.5, fontweight="bold")
            ax.text(0.5 + w / 2 + 0.03, y - i - 0.62, txt, ha="left",
                    va="center", fontsize=8.5)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, y)
    ax.set_title("A  from raw reads to IG transcripts", fontsize=11,
                 fontweight="bold", loc="left")


def panel_usage(ax, usage, locus, top_n, vmax, letter, source=""):
    rows = [r for r in sorted(usage, key=lambda r: -int(r["total_transcripts"]))
            if int(r["total_transcripts"]) > 0][:top_n]
    if not rows:
        ax.axis("off")
        ax.text(.5, .5, f"no {locus} usage", ha="center", transform=ax.transAxes)
        return
    names = [short(r["gene"]) for r in rows][::-1]
    vals = [int(r["total_transcripts"]) for r in rows][::-1]
    cmap, norm = locus_cmap(locus), make_norm(vmax)
    ax.barh(range(len(vals)), vals, color=[cmap(norm(v)) for v in vals],
            edgecolor="black", lw=.4)
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(names, fontsize=7)
    ax.set_xlabel("transcripts", fontsize=9)
    ax.set_title(f"{letter}  top {len(rows)} {locus} V genes by usage"
                 + (f"\n{source}" if source else ""),
                 fontsize=11, fontweight="bold", loc="left")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for i, v in enumerate(vals):
        ax.text(v, i, f" {v}", va="center", fontsize=6.5)


def panel_methods(ax, merge, reasons, letter):
    """
    Each method's TOTAL split into shared and unique, rather than three separate
    bars. Drawing "both" as its own bar invites reading it as a third method;
    stacking makes it clear that the 504 shared transcripts are counted inside
    both totals, and that immunotools' set sits almost entirely inside
    minimap2's.
    """
    both = int(merge["in_both"])
    mm_only = int(merge["minimap2_only"])
    it_only = int(merge["immunotools_only"])
    mm_tot, it_tot = both + mm_only, both + it_only

    from gc_palette import BOTH as SHARED_C
    UNIQ_MM, UNIQ_IT = LOCUS_C["IGH"], LOCUS_C["IGL"]
    ax.barh([1], [both], color=SHARED_C, edgecolor="black", lw=.5)
    ax.barh([1], [mm_only], left=[both], color=UNIQ_MM, edgecolor="black", lw=.5)
    ax.barh([0], [both], color=SHARED_C, edgecolor="black", lw=.5)
    ax.barh([0], [it_only], left=[both], color=UNIQ_IT, edgecolor="black", lw=.5)

    ax.text(both / 2, 1, f"{both}", ha="center", va="center", fontsize=8.5,
            color="white", fontweight="bold")
    ax.text(both + mm_only / 2, 1, f"+{mm_only}", ha="center", va="center",
            fontsize=8.5, fontweight="bold")
    ax.text(both / 2, 0, f"{both}", ha="center", va="center", fontsize=8.5,
            color="white", fontweight="bold")
    # the immunotools-unique slice is only 2 transcripts wide, so its "+2" and
    # its total would print on top of each other; combine them into one label
    ax.text(mm_tot, 1, f"   = {mm_tot}", ha="left", va="center", fontsize=9,
            fontweight="bold")
    ax.text(it_tot, 0, f"   +{it_only}  = {it_tot}", ha="left", va="center",
            fontsize=9, fontweight="bold")

    ax.set_yticks([1, 0])
    ax.set_yticklabels(["minimap2\nfilter", "immunotools\nvj_finder"], fontsize=8)
    ax.set_xlabel("IG transcripts", fontsize=9)
    ax.set_xlim(0, mm_tot * 1.30)
    ax.set_ylim(-7.6, 1.9)
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    ax.set_title(f"{letter}  how the two IG filters compare", fontsize=11,
                 fontweight="bold", loc="left")
    ax.legend(handles=[Patch(facecolor=SHARED_C, label=f"found by both ({both})"),
                       Patch(facecolor=UNIQ_MM, label=f"minimap2 only ({mm_only})"),
                       Patch(facecolor=UNIQ_IT, label=f"immunotools only ({it_only})")],
              loc="lower right", fontsize=7.5, framealpha=.9)

    pct = both / it_tot * 100 if it_tot else 0
    note = (f"immunotools' set is {pct:.1f}% contained within minimap2's\n"
            f"({both} of its {it_tot}); union = {both + mm_only + it_only}.\n\n"
            "minimap2 filter — align to germline V only; keep if\n"
            "identity >=0.70 and >=150 bp of V covered.\n"
            "immunotools/vj_finder — needs a V AND a J hit,\n"
            "V >=250 bp, J >=30 bp, plus uncovered-end limits.\n"
            "Requiring a J is why it recovers fewer transcripts,\n"
            "and the IGH J available here is a Gallus stand-in.")
    ax.text(0.0, -1.9, note, transform=ax.get_yaxis_transform(),
            fontsize=7.4, va="top", ha="left", clip_on=False)
    if reasons:
        txt = ("why immunotools rejected transcripts:\n" +
               "\n".join(f"  {k.replace('_', ' ').lower()}: {v}"
                          for k, v in sorted(reasons.items(),
                                             key=lambda kv: -kv[1])[:4]))
        ax.text(0.0, -5.3, txt, transform=ax.get_yaxis_transform(),
                fontsize=7.4, va="top", ha="left", family="monospace",
                clip_on=False)


def panel_map(ax, genes, locus, j_pos, vmax, letter):
    pos = [g["pos"] for g in genes]
    lo, hi = min(pos + [j_pos]), max(pos + [j_pos])
    pad = (hi - lo) * 0.02
    cmap, norm = locus_cmap(locus), make_norm(vmax)

    TALL, SHORT = 0.86, 0.34
    ax.axhline(0, color="black", lw=2.0, zorder=2)
    for g in genes:
        n = g["n"]
        col = cmap(norm(n)) if n > 0 else ZERO_C
        h = TALL if g["rss"] else SHORT
        sign = 1 if g["strand"] == "+" else -1
        ax.plot([g["pos"], g["pos"]], [0, sign * h], color=col,
                lw=2.6 if g["rss"] else 1.8, solid_capstyle="butt",
                zorder=4 if n > 0 else 3)
    ax.plot([j_pos], [0], marker="D", ms=9, color=SEGMENT["J"], zorder=6)
    ax.annotate("J", (j_pos, 0), textcoords="offset points", xytext=(0, -21),
                ha="center", fontsize=9, color=SEGMENT["J"], fontweight="bold")

    ax.set_xlim(lo - pad, hi + pad)
    ax.set_ylim(-1.12, 1.12)
    ax.set_yticks([])
    for s in ("left", "right", "top"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_position(("outward", 8))
    ax.ticklabel_format(axis="x", style="plain", useOffset=False)
    ax.set_xticklabels([f"{int(t):,}" for t in ax.get_xticks()], fontsize=8)
    ax.set_xlabel("contig position (bp)", fontsize=9, loc="left")
    n_rss = sum(1 for g in genes if g["rss"])
    n_plus = sum(1 for g in genes if g["strand"] == "+")
    ax.set_title(
        f"{letter}  {locus} locus — {len(genes)} V genes · {n_rss} with an RSS · "
        f"{n_plus} on + strand\n"
        f"up/down = strand · tall/short = RSS present/absent · colour = transcripts as parent",
        fontsize=10.5, fontweight="bold", loc="left")

    # below the axis rather than inside it: the tick region is full of genes
    cax = ax.inset_axes([0.76, -0.40, 0.22, 0.07])
    cb = plt.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=cax,
                      orientation="horizontal")
    cb.set_label(f"{locus} transcripts as parent", fontsize=7)
    cb.ax.tick_params(labelsize=6)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--skera-summary")
    ap.add_argument("--prefilter-stats", required=True)
    ap.add_argument("--filter-stats", required=True)
    ap.add_argument("--merge-stats")
    ap.add_argument("--immunotools-filtering", nargs="*", default=[],
                    help="filtering_info.txt files from diversity_analyzer")
    ap.add_argument("--igh-usage", help="legacy gene_usage.tsv (raw PAF identity)")
    ap.add_argument("--igl-usage", help="legacy gene_usage.tsv (raw PAF identity)")
    ap.add_argument("--usage-assignments",
                    help="UNCONSTRAINED transcript assignment TSV from "
                         "gc_call_functional_genes.py. Preferred over the legacy "
                         "gene_usage tables: those take the top hit by raw PAF "
                         "identity, with no junction masking and no coverage "
                         "floor, so they can name a different top gene than the "
                         "rest of the analysis. Scoring usage the same way "
                         "everywhere is the point.")
    ap.add_argument("--rss-annotation", required=True)
    ap.add_argument("--assignments", required=True)
    ap.add_argument("--igh-j-pos", type=int, required=True)
    ap.add_argument("--igl-j-pos", type=int, required=True)
    ap.add_argument("--top-n", type=int, default=20)
    ap.add_argument("--sample", default="")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    pre, filt = read_kv(args.prefilter_stats), read_kv(args.filter_stats)
    skera = {}
    if args.skera_summary:
        for line in open(args.skera_summary):
            if "," in line:
                k, v = line.rstrip("\n").rsplit(",", 1)
                skera[k] = v

    flnc, ig_flnc = int(pre["flnc_reads"]), int(pre["ig_flnc_reads"])
    total_tx, ig_tx = int(filt["total_transcripts"]), int(filt["ig_transcripts"])

    steps = []
    if skera:
        steps.append(("HiFi reads", int(skera["Input Reads"]), "raw"))
        steps.append(("segmented reads (skera)",
                      int(skera["Segmented Reads (S-Reads)"]),
                      f"×{float(skera['Mean Array Size (Concatenation Factor)']):.1f} array"))
    steps += [("FLNC reads", flnc, ""),
              ("IG-screened FLNC", ig_flnc, f"{ig_flnc/flnc*100:.4f}% of FLNC"),
              ("clustered transcripts", total_tx, ""),
              ("IG transcripts", ig_tx, f"{ig_tx/total_tx*100:.1f}% of clustered")]

    if args.usage_assignments:
        rows = read_tsv(args.usage_assignments)
        counts = {}
        for r in rows:
            counts.setdefault(r["locus"], Counter())[r["best_gene"]] += 1
        def as_usage(locus):
            return [{"gene": g, "total_transcripts": str(n)}
                    for g, n in counts.get(locus, Counter()).most_common()]
        igh_u, igl_u = as_usage("IGH"), as_usage("IGL")
        usage_src = "best-matching V gene (junction-aware)"
    else:
        if not (args.igh_usage and args.igl_usage):
            raise SystemExit("give --usage-assignments, or both --igh-usage and --igl-usage")
        igh_u, igl_u = read_tsv(args.igh_usage), read_tsv(args.igl_usage)
        usage_src = "top PAF hit (legacy scoring)"
    rss = {r["gene"]: r for r in read_tsv(args.rss_annotation)}
    assigns = read_tsv(args.assignments)
    merge = read_kv(args.merge_stats) if args.merge_stats else None

    reasons = {}
    for f in args.immunotools_filtering:
        try:
            for line in open(f):
                p = line.rstrip("\n").split("\t")
                if len(p) >= 2:
                    reasons[p[1]] = reasons.get(p[1], 0) + 1
        except OSError:
            pass

    usage_n = {r["gene"]: int(r["total_transcripts"]) for r in igh_u + igl_u}
    vmax = {loc: max([usage_n.get(g, 0) for g, r in rss.items()
                      if r["locus"] == loc] + [1]) for loc in ("IGH", "IGL")}

    def genes_for(locus):
        return sorted(({"pos": int(r["pos"]), "strand": r["strand"],
                        "rss": r["rss_state"] == "rss_present",
                        "n": usage_n.get(g, 0)}
                       for g, r in rss.items() if r["locus"] == locus),
                      key=lambda d: d["pos"])

    fig = plt.figure(figsize=(16, 19))
    gs = fig.add_gridspec(4, 3, height_ratios=[2.4, 2.9, 1.35, 1.35],
                          hspace=0.62, wspace=0.34)

    panel_funnel(fig.add_subplot(gs[0, 0]), steps)

    axb = fig.add_subplot(gs[0, 1])
    igh_n, igl_n = int(filt.get("best_locus_IGH", 0)), int(filt.get("best_locus_IGL", 0))
    axb.bar(["IGH", "IGL"], [igh_n, igl_n],
            color=[LOCUS_C["IGH"], LOCUS_C["IGL"]], edgecolor="black", lw=.6)
    for i, v in enumerate([igh_n, igl_n]):
        axb.text(i, v, f"\n{v}\n({v/(igh_n+igl_n)*100:.0f}%)", ha="center",
                 va="bottom", fontsize=9)
    axb.set_ylabel("IG transcripts", fontsize=9)
    axb.set_ylim(0, max(igh_n, igl_n) * 1.35)
    for s in ("top", "right"):
        axb.spines[s].set_visible(False)
    axb.set_title("B  locus split", fontsize=11, fontweight="bold", loc="left")

    axc = fig.add_subplot(gs[0, 2])
    # IGL has ~7x more transcripts, so IGH is drawn LAST and sits on top;
    # otherwise the smaller distribution is buried.
    for locus, alpha, z in (("IGL", .62, 2), ("IGH", .92, 3)):
        vals = [float(a["identity"]) for a in assigns if a["locus"] == locus]
        if vals:
            axc.hist(vals, bins=np.arange(.85, 1.005, .01), alpha=alpha,
                     color=LOCUS_C[locus], edgecolor="black", lw=.5,
                     zorder=z, label=f"{locus} (n={len(vals)})")
    handles, labels = axc.get_legend_handles_labels()
    order = [labels.index(l) for l in sorted(labels)]      # IGH first in legend
    axc.legend([handles[i] for i in order], [labels[i] for i in order], fontsize=8)
    axc.set_xlabel("identity to assigned germline V", fontsize=9)
    axc.set_ylabel("transcripts", fontsize=9)
    for s in ("top", "right"):
        axc.spines[s].set_visible(False)
    axc.set_title("C  divergence from germline\n"
                  "(SHM + gene conversion + cross-individual)",
                  fontsize=10.5, fontweight="bold", loc="left")

    panel_usage(fig.add_subplot(gs[1, 0]), igh_u, "IGH", args.top_n, vmax["IGH"], "D",
                usage_src)
    panel_usage(fig.add_subplot(gs[1, 1]), igl_u, "IGL", args.top_n, vmax["IGL"], "E",
                usage_src)
    axf = fig.add_subplot(gs[1, 2])
    if merge:
        panel_methods(axf, merge, reasons, "F")
    else:
        axf.axis("off")

    panel_map(fig.add_subplot(gs[2, :]), genes_for("IGH"), "IGH",
              args.igh_j_pos, vmax["IGH"], "G")
    panel_map(fig.add_subplot(gs[3, :]), genes_for("IGL"), "IGL",
              args.igl_j_pos, vmax["IGL"], "H")

    handles = [
        Line2D([], [], color="black", lw=2.6, label="tall = has RSS (can be rearranged)"),
        Line2D([], [], color="black", lw=1.8, label="short = no RSS (donor only)"),
        Patch(facecolor=ZERO_C, label="grey = no transcripts"),
        Line2D([], [], color="none", label="up = + strand   ·   down = − strand"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=9,
               frameon=True, bbox_to_anchor=(0.5, 0.008))
    fig.suptitle(f"IsoSeq IG overview — {args.sample}", fontsize=15,
                 fontweight="bold", y=0.983)
    save_figure(fig, args.out)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
