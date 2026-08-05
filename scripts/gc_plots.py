"""
Figures for the gene conversion analysis.

Produces, per locus:

  <locus>_donor_network.pdf   genes laid out along the locus with an arc for
                              every donor -> recipient relationship, labelled
                              with the number of supporting tracts.  Arcs that
                              respect the recombination topology are drawn
                              above the axis, impossible ones below, so the
                              false-positive load is visible at a glance.
  <locus>_geneconv_report.pdf multi-panel summary: donor/recipient matrix,
                              locus architecture, topology sweep, tract
                              positions and evidence strength.
"""
import argparse
import math
from collections import Counter, defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import ConnectionStyle, FancyArrowPatch, Patch, Rectangle

from gc_palette import save_figure, YES, NO, GREY, GREY_DARK, INK, SEGMENT, locus_ramp, HEATMAP

ALLOWED_C = YES           # topologically possible
IMPOSSIBLE_C = NO         # donor was deleted by the rearrangement
FUNC_C = INK
PSEUDO_C = GREY

# Gene colours encode RSS state, in the V-segment navy so they never compete
# with the teal/rose used for arrow topology.
RSS_COLORS = {
    "rss_present": SEGMENT["V"],  # RSS recorded: can be rearranged
    "rss_absent": GREY,           # none recorded: donor-only pseudogene
}


def locus_ramp_from(colour):
    from gc_palette import ramp
    return ramp(colour)


def read_tsv(path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(header, line.rstrip("\n").split("\t")))
                for line in fh if line.strip()]


def short(name):
    """VGP_..._IGL.8749176.CM036732.1.V.True.-  ->  8749176"""
    try:
        return name.split(".")[1]
    except IndexError:
        return name


# ─── figure 1: donor -> recipient arcs along the locus ────────────────────────

def donor_network(tracts, pool, genes, locus, j_pos, out, rss=None,
                  min_count=1, label_top=14):
    order = sorted(pool, key=lambda g: int(pool[g]["pos"]))
    xof = {g: i for i, g in enumerate(order)}
    positions = [int(pool[g]["pos"]) for g in order]

    pairs = Counter()
    allowed_flag = {}
    for t in tracts:
        if t["significant"] != "True":
            continue
        p, d = t["parent"], t["donor"]
        if p not in xof or d not in xof:
            continue
        pairs[(d, p)] += 1
        allowed_flag[(d, p)] = t["donor_allowed"] == "True"

    fig, ax = plt.subplots(figsize=(15, 8.5))

    # Gene track. Colour encodes RSS state -- whether the gene can be rearranged
    # at all -- and marker size encodes whether it is expressed. Those are two
    # independent facts, so they get two independent visual channels and two
    # legend blocks.
    def gene_colour(g):
        if rss is None:
            return FUNC_C if genes.get(g, {}).get("annotated_productive") == "True" else PSEUDO_C
        return RSS_COLORS.get(rss.get(g, {}).get("rss_state", "rss_absent"), PSEUDO_C)

    # Expression is read from the UNCONSTRAINED assignment on purpose, so the
    # size channel means the same thing in every figure. Under a constrained
    # parent model only the parent can receive transcripts, which would make
    # size a restatement of the model rather than a property of the data.
    # Marker SHAPE encodes transcriptional orientation relative to J, which is
    # what sets the recombination mechanism: same orientation -> deletion (the
    # intervening DNA is excised, so donors between V and J are lost); opposite
    # -> inversion (everything is retained).
    # J sits at one end of the array; a deletional gene is drawn pointing toward
    # it, an inversional gene away from it.
    j_dir = "left" if j_pos < min(positions) else "right"
    for g in order:
        x = xof[g]
        ntx = int(genes.get(g, {}).get("n_transcripts", 0) or 0)
        expressed = ntx > 0
        mech = pool[g].get("mechanism", "")
        # point the marker toward J when co-oriented with it (deletional)
        if mech == "deletion":
            mk = "<" if j_dir == "left" else ">"
        elif mech == "inversion":
            mk = ">" if j_dir == "left" else "<"
        else:
            mk = "o"
        ax.plot([x], [0], marker=mk, ms=13 if expressed else 7.5,
                color=gene_colour(g),
                markeredgecolor="black" if expressed else "none",
                markeredgewidth=1.4, zorder=5)
    ax.axhline(0, color="black", lw=1.2, zorder=1)

    # J marker at whichever end it sits
    j_left = j_pos < min(positions)
    jx = -1.4 if j_left else len(order) + 0.4
    ax.plot([jx], [0], marker="D", ms=12, color=SEGMENT["J"], zorder=6)
    ax.annotate(f"J\n{j_pos:,}", (jx, 0), textcoords="offset points",
                xytext=(0, -34), ha="center", fontsize=10, color=SEGMENT["J"],
                fontweight="bold")

    maxc = max(pairs.values()) if pairs else 1
    top = {k for k, _ in pairs.most_common(label_top)}

    n_allowed = n_impossible = 0
    labels = []
    for (d, p), c in sorted(pairs.items(), key=lambda kv: kv[1]):
        if c < min_count:
            continue
        ok = allowed_flag[(d, p)]
        n_allowed += c if ok else 0
        n_impossible += 0 if ok else c
        x0, x1 = xof[d], xof[p]
        sign = 1 if ok else -1
        lw = 0.7 + 3.6 * (c / maxc)
        # arc3's curvature is relative to the direction of travel, so the sign
        # has to be flipped for right-to-left arcs or "possible" and
        # "impossible" end up on the same side of the axis.
        travel = 1 if x1 > x0 else -1
        rad = -0.42 * sign * travel
        arrow = FancyArrowPatch(
            (x0, 0), (x1, 0),
            connectionstyle=f"arc3,rad={rad}",
            arrowstyle="-|>", mutation_scale=11 + 7 * (c / maxc),
            lw=lw, color=ALLOWED_C if ok else IMPOSSIBLE_C,
            alpha=0.75, zorder=3,
        )
        ax.add_patch(arrow)
        if (d, p) in top:
            labels.append((x0, x1, rad, c, ALLOWED_C if ok else IMPOSSIBLE_C))

    # With ~130 genes every label cannot be legible; show every Nth and keep
    # minor ticks for the rest so positions stay readable.
    step = 1 if len(order) <= 40 else (2 if len(order) <= 70 else 4)
    ax.set_xticks(range(0, len(order), step))
    ax.set_xticklabels([short(order[i]) for i in range(0, len(order), step)],
                       rotation=90, fontsize=8 if step == 1 else 6)
    ax.set_xticks(range(len(order)), minor=True)
    ax.set_yticks([])
    for s in ("left", "right", "top"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.set_xlim(-2.4, len(order) + 1.4)
    lim = 1.55
    ax.set_ylim(-lim, lim)

    # ── arc labels ────────────────────────────────────────────────────────────
    # These have to come after the limits are fixed.  FancyArrowPatch builds its
    # connection in DISPLAY space, not data space, so the apex of an arc depends
    # on the axes aspect and cannot be predicted from the data coordinates alone.
    # Ask matplotlib for the same Arc3 it drew, then pull the curve back into
    # data coordinates: Arc3 yields a quadratic Bezier [P0, control, P2], which
    # is exact to evaluate.  (The transform is dpi-invariant, so this survives
    # savefig at any resolution.)
    to_disp, to_data = ax.transData, ax.transData.inverted()
    placed = []
    # Try the apex first, then walk symmetrically outwards along the same curve.
    # Sliding a number along its own arc keeps it unambiguously attached to that
    # arc; lifting it off would put it in the empty band between arcs, which is
    # what made the old labels look unmoored.
    offsets = (0.5, 0.40, 0.60, 0.31, 0.69, 0.23, 0.77, 0.16, 0.84)
    for x0, x1, rad, c, col in sorted(labels, key=lambda t: -t[3]):
        bez = to_data.transform(
            ConnectionStyle.Arc3(rad=rad)(
                to_disp.transform((x0, 0.0)),
                to_disp.transform((x1, 0.0)),
            ).vertices
        )
        (bx0, by0), (cx, cy), (bx1, by1) = bez[0], bez[1], bez[2]
        lx = ly = None
        for t in offsets:
            u = 1.0 - t
            px_ = u * u * bx0 + 2 * u * t * cx + t * t * bx1
            py_ = u * u * by0 + 2 * u * t * cy + t * t * by1
            if lx is None:
                lx, ly = px_, py_          # fall back to the apex
            if all(abs(px_ - qx) > 0.9 or abs(py_ - qy) > 0.1
                   for qx, qy in placed):
                lx, ly = px_, py_
                break
        placed.append((lx, ly))
        ax.annotate(str(c), (lx, ly), ha="center", va="center", fontsize=8.5,
                    color=col, fontweight="bold",
                    bbox=dict(boxstyle="round,pad=0.16", fc="white",
                              ec="none", alpha=0.88), zorder=7)

    ax.text(0.5, 0.965,
            "donor → recipient, arc respects recombination topology",
            transform=ax.transAxes, ha="center", fontsize=10, color=ALLOWED_C)
    ax.text(0.5, 0.035,
            "donor was DELETED by the rearrangement — impossible",
            transform=ax.transAxes, ha="center", fontsize=10, color=IMPOSSIBLE_C)

    tot = n_allowed + n_impossible
    frac = n_impossible / tot if tot else float("nan")
    # Only call them "significant" if the tracts actually carry a p-value.
    # BrepConvert-derived tracts do not, so labelling them significant would
    # imply a filter that was never applied.
    scored = any(t.get("p_corrected", "NA") not in ("NA", "") for t in tracts)
    what = "significant tracts" if scored else "donor-calls, no significance filter"
    ax.set_title(
        f"{locus} — gene conversion donor→recipient relationships "
        f"({tot} {what})\n"
        f"{n_impossible}/{tot} = {frac:.0%} use a donor that recombination "
        f"had already deleted",
        fontsize=13, fontweight="bold", pad=16)

    # Line2D handles rather than Patch, so the legend shows the real marker
    # sizes -- size is a data channel here and a coloured square cannot show it.
    def dot(color, ms, edge, label):
        return Line2D([], [], marker="o", linestyle="none", markersize=ms,
                      color=color, markeredgecolor=edge,
                      markeredgewidth=1.4 if edge != "none" else 0, label=label)

    if rss is None:
        colour_handles = [dot(FUNC_C, 9, "none", "annotated productive"),
                          dot(PSEUDO_C, 9, "none", "annotated pseudogene")]
        colour_title = "gene colour — annotation"
    else:
        colour_handles = [
            dot(RSS_COLORS["rss_present"], 9, "none", "has RSS — can be rearranged"),
            dot(RSS_COLORS["rss_absent"], 9, "none", "no RSS — donor only"),
        ]
        colour_title = "gene colour — RSS state"

    # One legend, placed BELOW the axes rather than inside them. Arcs can reach
    # anywhere in the plotting area, so any in-axes placement eventually covers
    # data; putting it outside is the only placement that is safe for every
    # dataset. Section headers are blank handles so the groups stay readable
    # in a single flat row.
    def spacer(label):
        return Line2D([], [], linestyle="none", marker="", label=label)

    def tri(mk, label):
        return Line2D([], [], marker=mk, linestyle="none", markersize=10,
                      color="#777777", markeredgecolor="none", label=label)

    handles = [spacer("$\\bf{gene\\ colour}$")] + colour_handles + [
        spacer("$\\bf{orientation\\ vs\\ J}$"),
        tri("<" if j_dir == "left" else ">", "same as J → deletion"),
        tri(">" if j_dir == "left" else "<", "opposite → inversion"),
        spacer("$\\bf{gene\\ size}$"),
        # neutral fill, not white: a white dot on a white legend is invisible
        dot("#dddddd", 11, "black", "≥1 transcript best-matches it"),
        dot("#dddddd", 6.5, "#777777", "no transcript matches it"),
        spacer("$\\bf{arrow\\ colour}$"),
        Patch(facecolor=ALLOWED_C, label="donor available"),
        Patch(facecolor=IMPOSSIBLE_C, label="donor deleted (impossible)"),
    ]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.30),
              ncol=4, fontsize=9, framealpha=0.95, borderpad=0.9,
              columnspacing=1.8, handletextpad=0.7)

    ax.set_xlabel("V gene (contig position, ordered; J-proximal end nearest J marker)",
                  fontsize=10)
    save_figure(fig, out)
    plt.close(fig)
    return tot, n_impossible


# ─── figure 2: multi-panel report ─────────────────────────────────────────────

def report(tracts, pool, genes, arch, topo, locus, out):
    order = sorted(pool, key=lambda g: int(pool[g]["pos"]))
    idx = {g: i for i, g in enumerate(order)}
    n = len(order)

    fig = plt.figure(figsize=(16, 11))
    gs = fig.add_gridspec(2, 3, hspace=0.38, wspace=0.28)

    # (a) donor x recipient matrix
    ax = fig.add_subplot(gs[0, 0])
    M = np.zeros((n, n))
    for t in tracts:
        if t["significant"] != "True":
            continue
        if t["parent"] in idx and t["donor"] in idx:
            M[idx[t["donor"]], idx[t["parent"]]] += 1
    # viridis reversed: yellow = low, dark purple = high
    im = ax.imshow(M, cmap=HEATMAP, origin="lower", aspect="auto")
    # hatch the impossible half-plane
    for p in order:
        allowed = set() if pool[p]["allowed_donors"] == "NONE" else set(
            pool[p]["allowed_donors"].split(","))
        if pool[p]["mechanism"] != "deletion":
            continue
        for d in order:
            if d != p and d not in allowed:
                ax.add_patch(Rectangle((idx[p] - .5, idx[d] - .5), 1, 1,
                                       fill=False, hatch="///", lw=0,
                                       edgecolor=IMPOSSIBLE_C, alpha=.45))
    ax.set_xlabel("recipient (parent)")
    ax.set_ylabel("donor")
    ax.set_title("A  significant tracts per donor–recipient pair\n"
                 "hatching = topologically impossible", fontsize=10)
    ax.set_xticks(range(0, n, max(1, n // 12)))
    ax.set_xticklabels([short(order[i]) for i in range(0, n, max(1, n // 12))],
                       rotation=90, fontsize=6)
    ax.set_yticks(range(0, n, max(1, n // 12)))
    ax.set_yticklabels([short(order[i]) for i in range(0, n, max(1, n // 12))],
                       fontsize=6)
    fig.colorbar(im, ax=ax, shrink=.8, label="tracts")

    # (b) locus architecture: donor supply of functional vs pseudogenes
    ax = fig.add_subplot(gs[0, 1])
    if arch:
        fx = [int(r["allowed_donors"]) for r in arch if r["annot_productive"] == "True"]
        px = [int(r["allowed_donors"]) for r in arch if r["annot_productive"] != "True"]
        rng = np.random.default_rng(0)
        for vals, ypos, col, lab in ((fx, 1, FUNC_C, "annotated\nfunctional"),
                                     (px, 0, PSEUDO_C, "annotated\npseudogene")):
            if not vals:
                continue
            ax.scatter(vals, ypos + rng.uniform(-.13, .13, len(vals)),
                       color=col, s=55, edgecolor="black", lw=.5, zorder=3, label=lab)
            ax.plot([np.mean(vals)] * 2, [ypos - .28, ypos + .28],
                    color="black", lw=2.2, zorder=4)
        ax.set_yticks([0, 1])
        ax.set_yticklabels(["pseudogene", "functional"])
        ax.set_ylim(-.55, 1.55)
        ax.set_xlabel("allowed donors (higher = more J-proximal)")
        mf = np.mean(fx) if fx else float("nan")
        mp = np.mean(px) if px else float("nan")
        ax.set_title("B  functional genes sit where donors are plentiful\n"
                     f"mean {mf:.1f} vs {mp:.1f}  (bars = group means)", fontsize=10)
        ax.grid(axis="x", alpha=.25)

    # (c) power-matched topology test for THIS method's own calls.
    # Parents whose donor pool has nothing deleted cannot produce an impossible
    # call, so including them dilutes the test to meaninglessness -- they are
    # excluded here and the exclusion is stated on the panel.
    ax = fig.add_subplot(gs[0, 2])
    p_imp = {}
    for g, r in pool.items():
        if r["mechanism"] != "deletion":
            continue
        a, d = int(r["n_allowed_donors"]), int(r["n_deleted"])
        if a + d > 0:
            p_imp[g] = d / (a + d)
    powered = {g for g, v in p_imp.items() if v > 0}
    sel = [t for t in tracts if t["significant"] == "True"
           and t["mechanism"] == "deletion" and t["parent"] in powered]
    if sel:
        obs = sum(1 for t in sel if t["donor_allowed"] != "True")
        probs = [p_imp[t["parent"]] for t in sel]
        exp = sum(probs)
        sd = math.sqrt(sum(q * (1 - q) for q in probs))
        m = len(sel)
        pv = float("nan")
        if sd > 0:
            from scipy import stats as _st
            pv = _st.norm.cdf((obs + 0.5 - exp) / sd)
        ax.bar([0, 1], [obs / m, exp / m],
               color=[IMPOSSIBLE_C, "#888888"], width=.6,
               edgecolor="black", lw=.6)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["observed", "expected\nif random"])
        ax.set_ylabel("fraction of calls with an impossible donor")
        ax.set_ylim(0, max(.6, min(1.0, max(obs / m, exp / m) * 1.35)))
        for xi, v in enumerate([obs / m, exp / m]):
            ax.text(xi, v + .015, f"{v:.1%}", ha="center", fontsize=10,
                    fontweight="bold")
        verdict = ("beats chance" if (pv == pv and pv < 0.05 and obs / m < exp / m)
                   else "no better than chance")
        # Each listed donor counts as one call. Where a method reports an
        # ambiguous donor set this is the strict reading; scoring an event as
        # correct when ANY listed donor is possible is more generous and can
        # give a very different answer, so the view is stated explicitly.
        ax.set_title(f"C  topology control, power-matched\n"
                     f"n={m} donor-calls on {len(powered)} informative parents — "
                     f"{verdict} (p={pv:.2g})\n"
                     f"each listed donor counted separately (strict reading)",
                     fontsize=9)
    else:
        ax.text(.5, .5, "no calls on informative parents",
                ha="center", va="center", transform=ax.transAxes)
        ax.set_title("C  topology control, power-matched", fontsize=10)

    # (d) tract start positions along V
    ax = fig.add_subplot(gs[1, 0])
    sig = [t for t in tracts if t["significant"] == "True"]
    if sig:
        starts = [int(t["start"]) for t in sig]
        ok = [int(t["start"]) for t in sig if t["donor_allowed"] == "True"]
        bad = [int(t["start"]) for t in sig if t["donor_allowed"] != "True"]
        bins = np.arange(0, max(starts) + 12, 10)
        ax.hist([ok, bad], bins=bins, stacked=True,
                color=[ALLOWED_C, IMPOSSIBLE_C], label=["possible", "impossible"])
        ax.set_xlabel("tract start (bp along V gene)")
        ax.set_ylabel("tracts")
        ax.legend(fontsize=8)
        ax.set_title("D  where tracts fall along V", fontsize=10)

    # (e) evidence strength -- n_support is NA for BrepConvert-derived tracts,
    # which report a span rather than a count of diagnostic positions
    ax = fig.add_subplot(gs[1, 1])
    have_support = [t for t in sig if t.get("n_support", "NA") not in ("NA", "")]
    if have_support:
        sig_e = have_support
        ok = [int(t["n_support"]) for t in sig_e if t["donor_allowed"] == "True"]
        bad = [int(t["n_support"]) for t in sig_e if t["donor_allowed"] != "True"]
        mx = max([*ok, *bad]) if (ok or bad) else 1
        bins = np.arange(2.5, mx + 1.5, 1)
        ax.hist([ok, bad], bins=bins, stacked=True,
                color=[ALLOWED_C, IMPOSSIBLE_C], label=["possible", "impossible"])
        ax.set_xlabel("donor-diagnostic supporting positions per tract")
        ax.set_ylabel("tracts")
        ax.legend(fontsize=8)
        ax.set_title("E  impossible donors are as well supported\n"
                     "as possible ones — evidence strength does not separate them",
                     fontsize=10)
    elif sig:
        ok = [int(t["span_bp"]) for t in sig if t["donor_allowed"] == "True"]
        bad = [int(t["span_bp"]) for t in sig if t["donor_allowed"] != "True"]
        mx = max([*ok, *bad]) if (ok or bad) else 1
        ax.hist([ok, bad], bins=np.arange(0, min(mx, 80) + 3, 2), stacked=True,
                color=[ALLOWED_C, IMPOSSIBLE_C], label=["possible", "impossible"])
        ax.set_xlabel("tract span (bp)")
        ax.set_ylabel("tracts")
        ax.legend(fontsize=8)
        ax.set_title("E  tract length distribution", fontsize=10)

    # (f) per-gene usage vs donor supply
    ax = fig.add_subplot(gs[1, 2])
    if arch:
        x = [int(r["allowed_donors"]) for r in arch]
        y = [int(r["n_transcripts"]) for r in arch]
        c = [FUNC_C if r["annot_productive"] == "True" else PSEUDO_C for r in arch]
        ax.scatter(x, y, c=c, s=48, edgecolor="black", lw=.5)
        for r in arch:
            if int(r["n_transcripts"]) > 0:
                ax.annotate(short(r["gene"]),
                            (int(r["allowed_donors"]), int(r["n_transcripts"])),
                            fontsize=6, xytext=(3, 3), textcoords="offset points")
        ax.set_xlabel("allowed donors")
        ax.set_ylabel("transcripts assigned")
        ax.set_title("F  expression vs donor supply", fontsize=10)

    fig.suptitle(f"{locus} — gene conversion analysis", fontsize=15, fontweight="bold")
    save_figure(fig, out)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tracts", required=True)
    ap.add_argument("--donor-pool", required=True)
    ap.add_argument("--functional-genes", required=True)
    ap.add_argument("--architecture", help="*_locus_architecture_per_gene.tsv")
    ap.add_argument("--topology", help="*_topology_test.tsv")
    ap.add_argument("--locus", required=True)
    ap.add_argument("--j-pos", type=int, required=True)
    ap.add_argument("--rss-annotation",
                    help="rss_annotation.tsv from gc_rss_annotation.py. When "
                         "given, genes are coloured by RSS state instead of by "
                         "the annotation's Productive flag.")
    ap.add_argument("--out-network", required=True)
    ap.add_argument("--out-report", required=True)
    ap.add_argument("--out-matrix", help="Donor x parent matrix figure")
    args = ap.parse_args()

    tracts = read_tsv(args.tracts)
    pool = {r["rearranged_gene"]: r for r in read_tsv(args.donor_pool)}
    genes = {r["gene"]: r for r in read_tsv(args.functional_genes)
             if r["locus"] == args.locus}
    rss = None
    if args.rss_annotation:
        rss = {r["gene"]: r for r in read_tsv(args.rss_annotation)
               if r["locus"] == args.locus}
    arch = read_tsv(args.architecture) if args.architecture else []
    topo = read_tsv(args.topology) if args.topology else []

    tot, bad = donor_network(tracts, pool, genes, args.locus, args.j_pos,
                             args.out_network, rss=rss)
    report(tracts, pool, genes, arch, topo, args.locus, args.out_report)
    if args.out_matrix:
        if donor_matrix(tracts, pool, genes, args.locus, args.j_pos,
                        args.out_matrix, rss=rss):
            print(f"{args.locus}: wrote {args.out_matrix}")
        else:
            print(f"{args.locus}: no calls, matrix not written")
    print(f"{args.locus}: wrote {args.out_network} and {args.out_report} "
          f"({tot} significant tracts, {bad} impossible)")




# ─── figure 3: donor x recipient matrix, split by recombination mechanism ─────

def donor_matrix(tracts, pool, genes, locus, j_pos, out, rss=None):
    """
    One panel showing everything the arc plot shows, but readable when the array
    is large.

    Columns are the V genes actually used as a parent (usually few); rows are
    every V gene, ordered along the locus.  Cell colour is the number of tracts,
    on TWO scales: green where the donor survived recombination, red where it
    was deleted and the call is impossible.  So legality and abundance are the
    same visual channel and cannot be read apart by accident.

    The column annotation strip carries the recombination mechanism, which is
    the thing the arc plot could not show: an inversional parent retains the
    whole array, so no cell in its column can ever be red, while a deletional
    parent loses everything between itself and J.  The dashed step line marks
    each deletional parent's own position, i.e. the boundary beyond which its
    donors no longer exist.
    """
    order = sorted(pool, key=lambda g: int(pool[g]["pos"]))
    ypos = {g: i for i, g in enumerate(order)}

    counts, legal = Counter(), {}
    for t in tracts:
        if t["significant"] != "True":
            continue
        p, d = t["parent"], t["donor"]
        if p not in ypos or d not in ypos:
            continue
        counts[(d, p)] += 1
        legal[(d, p)] = t["donor_allowed"] == "True"

    # Columns are every CANDIDATE parent, not merely those that produced a call.
    # A candidate that yielded nothing is real information -- it says the gene
    # can rearrange but no conversion was detected on it -- and dropping such
    # columns makes "produced nothing" indistinguishable from "not a candidate".
    candidates = {p for _, p in counts}
    if rss is not None:
        candidates |= {g for g, r in rss.items()
                       if r.get("rss_state") == "rss_present" and g in pool}
    parents = sorted(candidates, key=lambda g: int(pool[g]["pos"]))
    if not parents:
        return False
    n_empty = sum(1 for p in parents if not any(pp == p for _, pp in counts))
    nrow, ncol = len(order), len(parents)

    pos_m = np.full((nrow, ncol), np.nan)
    neg_m = np.full((nrow, ncol), np.nan)
    for (d, p), c in counts.items():
        j = parents.index(p)
        (pos_m if legal[(d, p)] else neg_m)[ypos[d], j] = c

    fig, (axm, axc) = plt.subplots(
        2, 1, figsize=(max(5.5, 1.15 * ncol + 4.5), max(7.0, 0.085 * nrow + 3.4)),
        gridspec_kw={"height_ratios": [nrow, 2.2], "hspace": 0.02}, sharex=True)

    vmax = max(counts.values())
    im_ok = axm.imshow(pos_m, cmap=locus_ramp_from(YES), aspect="auto", origin="lower",
                       vmin=0, vmax=vmax)
    im_no = axm.imshow(neg_m, cmap=locus_ramp_from(NO), aspect="auto", origin="lower",
                       vmin=0, vmax=vmax)

    # boundary: for a deletional parent, donors beyond its own position are gone
    for j, p in enumerate(parents):
        if pool[p]["mechanism"] != "deletion":
            continue
        y = ypos[p]
        axm.plot([j - .5, j + .5], [y, y], color="black", lw=1.6,
                 linestyle="--", zorder=6)

    ystep = 1 if nrow <= 30 else (2 if nrow <= 60 else 4)
    yticks = list(range(0, nrow, ystep))
    axm.set_yticks(yticks)
    axm.set_yticklabels([short(order[i]) for i in yticks],
                        fontsize=7 if nrow <= 40 else 5.5)
    # RSS state is carried by the TICK LABEL COLOUR rather than an extra marker
    # column, which otherwise collides with the labels themselves.
    if rss is not None:
        for i, lab in zip(yticks, axm.get_yticklabels()):
            state = rss.get(order[i], {}).get("rss_state", "rss_absent")
            if state == "rss_present":
                lab.set_color(RSS_COLORS["rss_present"])
                lab.set_fontweight("bold")
    axm.set_ylabel("donor V gene (ordered along locus)", fontsize=10)
    axm.set_title(
        f"{locus} — donors used per rearranged parent\n"
        f"green = donor survived · red = donor deleted (impossible) · "
        f"dashed = deletion boundary\n"
        f"columns = all {len(parents)} candidate parents"
        + (f" ({n_empty} produced no calls)" if n_empty else ""),
        fontsize=12, fontweight="bold", pad=12)

    # mechanism strip
    MECH_C = {"deletion": SEGMENT["D"], "inversion": "#F8CD9C"}
    for j, p in enumerate(parents):
        axc.add_patch(Rectangle((j - .5, 0), 1, 1,
                                color=MECH_C.get(pool[p]["mechanism"], "#cccccc")))
        axc.text(j, 0.5, pool[p]["mechanism"][:3].upper(), ha="center",
                 va="center", fontsize=7.5, color="white", fontweight="bold")
        if not any(pp == p for _, pp in counts):
            # candidate parent with zero calls: mark it so the empty column reads
            # as a result rather than as missing data
            axm.text(j, nrow * 0.5, "no calls", rotation=90, ha="center",
                     va="center", fontsize=7.5, color="#999999", style="italic")
    axc.set_xlim(-.5, ncol - .5)
    axc.set_ylim(0, 1)
    axc.set_yticks([])
    axc.set_ylabel("mech.", fontsize=8, rotation=0, ha="right", va="center")
    axc.set_xticks(range(ncol))
    axc.set_xticklabels([short(p) for p in parents], rotation=90, fontsize=8)
    axc.set_xlabel("rearranged parent V gene", fontsize=10)
    for s in ("top", "right", "left"):
        axc.spines[s].set_visible(False)

    cb1 = fig.colorbar(im_ok, ax=[axm, axc], label="tracts — donor available",
                       fraction=0.04, pad=0.02, aspect=28)
    cb2 = fig.colorbar(im_no, ax=[axm, axc], label="tracts — impossible",
                       fraction=0.04, pad=0.06, aspect=28)
    for cb in (cb1, cb2):
        cb.ax.tick_params(labelsize=7)
        cb.set_label(cb.ax.get_ylabel(), fontsize=8)

    handles = [Patch(facecolor=MECH_C["deletion"],
                     label="deletion — donors between V and J are lost"),
               Patch(facecolor=MECH_C["inversion"],
                     label="inversion — whole array retained")]
    if rss is not None:
        handles.append(Line2D([], [], marker="s", linestyle="none", markersize=7,
                              color=RSS_COLORS["rss_present"],
                              label="donor label in blue = has RSS"))
    axc.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -3.6),
               ncol=1, fontsize=8.5, framealpha=0.95)

    save_figure(fig, out)
    plt.close(fig)
    return True


if __name__ == "__main__":
    main()
