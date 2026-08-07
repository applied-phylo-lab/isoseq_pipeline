"""
The summary figure: everything the analysis actually established, on one page.

Panels, in the order the argument runs:

  A/B locus architecture  -- both loci. Almost no V gene carries an RSS, so almost
                             none can be rearranged. That is the setup: diversity
                             cannot come from combinatorial V use, so it has to
                             come from somewhere else.
  C  sequence evidence    -- the raw alignments. A contiguous block where the
                             transcript abandons its parent and follows a donor.
                             This comes BEFORE the two panels that analyse tracts,
                             because both of them speak of differences inside
                             versus outside a tract, and the reader has to have
                             seen one before those words mean anything.
  D  donor network        -- IGL. Every pseudogene in the array donates into the
                             single functional gene, and not one call uses a donor
                             that recombination had already deleted.
  D  AID spectrum         -- AID initiates BOTH processes; what differs is how the
                             lesion is repaired. Under SHM the mutation happens AT
                             the AID-targeted base, so differences pile onto its
                             motifs. Under gene conversion the lesion is repaired
                             by copying a donor, so the resulting differences sit
                             wherever the donor happened to differ -- unrelated to
                             where AID bound. Outside-tract differences therefore
                             carry the targeting signature and inside-tract ones do
                             not, which separates the two processes without
                             trusting any donor assignment.

A and B describe the system; C shows the events; D and E are the evidence
that those events are conversion rather than hypermutation.
"""
import argparse
import statistics
import sys
from collections import Counter, defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import FuncFormatter
from matplotlib.lines import Line2D
from matplotlib.patches import ConnectionStyle, FancyArrowPatch, Patch, Rectangle
from Bio import Align

from gc_lib import read_fasta, parse_paf, projected_query
from gc_palette import (save_figure, LOCUS, SEGMENT, GREY, GREY_DARK, INK,
                        YES, NO, ramp)

DONOR_C = YES
OTHER_C = NO
UNINF_C = GREY


def L(letter):
    """Panel letter in bold, via mathtext, so the rest of the title stays normal.

    set_title takes a single fontweight, so a mixed-weight title needs either two
    overlaid text objects (fragile once titles wrap to two lines) or a mathtext
    span, which is what this is.
    """
    return r"$\bf{" + letter + r"}$"


def read_tsv(path):
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(hdr, ln.rstrip("\n").split("\t")))
                for ln in fh if ln.strip()]


def short(name):
    try:
        return name.split(".")[1]
    except IndexError:
        return name


def kv(path):
    d = defaultdict(dict)
    for r in read_tsv(path):
        d[r["class"]][r["metric"]] = r["value"]
    return d


# ─── A: locus architecture ───────────────────────────────────────────────────

def arch_extent(genes, usage, vmax, silent=0.048):
    """Vertical range a locus actually needs, per strand.

    IGL has one minus-strand gene, so a symmetric axis would spend half the panel
    on empty space. Returning the true extent lets the caller size the panels so
    that a given transcript count is the same number of millimetres in BOTH --
    which is the only way the two are actually comparable.
    """
    def h(n):
        return np.log1p(n) / np.log1p(vmax) if n else silent
    up = max([h(usage.get(g["gene"], 0)) for g in genes
              if g["strand"] == "+"] + [silent])
    dn = max([h(usage.get(g["gene"], 0)) for g in genes
              if g["strand"] != "+"] + [silent])
    return -(dn + 0.20), up + 0.22


def panel_architecture(ax, genes, locus, j_pos, usage, letter, vmax, ylim):
    genes = sorted(genes, key=lambda g: g["pos"])
    xs = [g["pos"] for g in genes]
    lo, hi = min(xs), max(xs)
    span = max(hi - lo, 1)
    col = LOCUS.get(locus, INK)
    # Height is the transcript count on a LOG axis, normalised by a vmax SHARED
    # between the two panels. Sharing it is what makes A and B comparable: the
    # same stem height means the same number of transcripts in both, so IGH's
    # ceiling of 7 is visibly nowhere near IGL's 233 rather than being rescaled
    # to fill its own panel. The top genes carry printed counts as well, since
    # log deliberately compresses the very difference it is being asked to show.
    # Silent genes sit ON the zero line, nudged just far enough off it to show
    # which strand they are on. The offset has to stay well under height(1)
    # (~0.13 on the shared scale) or a gene with no transcripts is drawn level
    # with the "1" tick and reads as having one.
    SILENT = 0.048

    def height(n):
        return np.log1p(n) / np.log1p(vmax)

    for g in genes:
        x = g["pos"]
        rss = g["rss"]
        n = usage.get(g["gene"], 0)
        sign = 1 if g["strand"] == "+" else -1
        if n:
            h = height(n)
            ax.plot([x, x], [0, sign * h], lw=1.9, color=col,
                    solid_capstyle="butt", zorder=3)
        else:
            # A silent gene is drawn as a MARKER offset to its own strand's side,
            # not as a stem. Strand is a property of the gene whether or not it is
            # expressed, so it has to stay visible -- but on a sqrt axis a stem
            # short enough to mean "zero" is indistinguishable from one meaning
            # "1-2 transcripts". A marker is sized in points rather than data
            # units, so it stays legible at any scaling while the offset still
            # carries the strand.
            h = SILENT
            if not rss:
                ax.plot([x], [sign * h], marker="o", ms=3.4, color="white",
                        markeredgecolor=GREY_DARK, markeredgewidth=0.9, zorder=4)
        if rss:
            ax.plot([x], [sign * h], marker="o", ms=5.0,
                    color=SEGMENT["V"], markeredgecolor="black",
                    markeredgewidth=0.6, zorder=4)

    # direct labels on the biggest genes -- the height carries the comparison,
    # the number carries the value
    top = sorted(genes, key=lambda g: -usage.get(g["gene"], 0))[:3]
    for g in top:
        n = usage.get(g["gene"], 0)
        if n < 3:
            continue
        sign = 1 if g["strand"] == "+" else -1
        ax.annotate(f"{n}", (g["pos"], sign * height(n)),
                    textcoords="offset points", xytext=(0, 6 * sign),
                    ha="center", va="bottom" if sign > 0 else "top",
                    fontsize=7, fontweight="bold", color=col, zorder=6,
                    bbox=dict(boxstyle="round,pad=0.12", fc="white",
                              ec="none", alpha=0.9))

    # Height now carries a number, so it needs a real axis. Stems run up for +
    # strand and down for -, so the scale is mirrored: the same tick value means
    # the same transcript count on either side of the line.
    # Ticks only on the upper half: the scale is mirrored, so labelling both
    # sides just prints every number twice, and in IGH (small range) the two
    # copies collide. Direction is strand, which the axis label states.
    # The SAME tick values in both panels, so a reader can look across. Ticks are
    # drawn on the minus side too wherever the axis reaches them -- previously the
    # lower half had gridlines but no numbers, which left it unreadable.
    refs = [r for r in (1, 5, 10, 50, 200) if r <= vmax]
    ticks, labels = [], []
    for r in refs:
        yr = height(r)
        if yr <= ylim[1]:
            ax.axhline(yr, color=GREY, lw=0.5, ls=":", zorder=1)
            ticks.append(yr)
            labels.append(str(r))
        if -yr >= ylim[0]:
            ax.axhline(-yr, color=GREY, lw=0.5, ls=":", zorder=1)
            ticks.append(-yr)
            labels.append(str(r))
    order = np.argsort(ticks)
    ax.set_yticks(np.array(ticks)[order])
    ax.set_yticklabels(np.array(labels)[order], fontsize=6.5)
    ax.set_ylabel("transcripts", fontsize=7)
    ax.tick_params(axis="y", length=2, pad=1.5)

    ax.axhline(0, color="black", lw=1.0, zorder=2)

    # J can sit far outside the V array (IGH: V genes end at 77 kb, J at 219 kb).
    # Drawing it to scale would compress every gene into the left third, so the
    # axis stops at the array and J is placed just beyond a break marker with its
    # real coordinate spelled out.
    pad = span * 0.03
    right = hi + pad
    # J can lie off either end: above the array (blackbird IGL) or below it
    # (tufted duck IGH, whose locus reads right-to-left). Placing the marker
    # unconditionally on the right would put it on the wrong side of the genes
    # and invert the geometry the panel is supposed to show.
    left = lo - pad
    if j_pos > hi + span * 0.15:
        jx = hi + span * 0.13
        brk = hi + span * 0.065
        right = hi + span * 0.20
    elif j_pos < lo - span * 0.15:
        jx = lo - span * 0.13
        brk = lo - span * 0.065
        left = lo - span * 0.20
        right = hi + pad
    else:
        jx = brk = None
        right = max(hi, j_pos) + pad
        left = min(lo, j_pos) - pad
    if brk is not None:
        ax.text(brk, 0, "//", ha="center", va="center",
                fontsize=11, color="black", zorder=6,
                bbox=dict(boxstyle="square,pad=0.06", fc="white", ec="none"))
    if jx is None:
        jx = j_pos
    ax.plot([jx], [0], marker="D", ms=9, color=SEGMENT["J"], zorder=5,
            clip_on=False)
    j_unit = "Mb" if abs(j_pos) >= 1e6 else "kb"
    j_div, j_dec = (1e6, 3) if j_unit == "Mb" else (1e3, 1)
    ax.annotate(f"J\n{j_pos / j_div:.{j_dec}f} {j_unit}", (jx, 0),
                textcoords="offset points", xytext=(0, -13), ha="center",
                va="top", fontsize=7.5, color=SEGMENT["J"], fontweight="bold")

    ax.set_xlim(left, right)
    ax.set_ylim(*ylim)
    for s in ("right", "top"):
        ax.spines[s].set_visible(False)
    ax.spines["left"].set_bounds(min(ticks), max(ticks))
    ax.spines["bottom"].set_position(("outward", 6))
    # Raw base pairs read badly at these magnitudes -- IGL's ticks came out as
    # 6.315 with a detached "1e6" in the corner. Pick the unit that gives short
    # numbers for THIS contig: kb while the array is under a megabase, Mb once it
    # is past one. The two panels are different contigs, so there is nothing to
    # compare between their absolute coordinates and no reason to force a shared
    # unit on them.
    if max(abs(min(lo, jx)), abs(right)) >= 1e6:
        unit, div, dec, jdec = "Mb", 1e6, 3, 3
    else:
        unit, div, dec, jdec = "kb", 1e3, 0, 1
    ax.xaxis.set_major_formatter(
        FuncFormatter(lambda v, _pos: f"{v / div:.{dec}f}"))
    ax.tick_params(axis="x", labelsize=7)
    ax.set_xlabel(f"contig position ({unit})", fontsize=8, loc="left")

    n_rss = sum(1 for g in genes if g["rss"])
    n_exp = sum(1 for g in genes if usage.get(g["gene"], 0) > 0)
    ax.set_title(f"{L(letter)}  {locus} — {len(genes)} V genes · "
                 f"{n_rss} with an RSS · {n_exp} expressed",
                 fontsize=10, loc="left", x=-0.048, color="black")


# ─── B: donor network ────────────────────────────────────────────────────────

def panel_network(ax, tracts, pool, usage, rssmap, locus, letter, title_x=0.0):
    order = sorted(pool, key=lambda g: int(pool[g]["pos"]))
    xof = {g: i for i, g in enumerate(order)}

    regions = defaultdict(list)
    for t in tracts:
        if t.get("significant") != "True":
            continue
        p, d = t["parent"], t["donor"]
        if p not in xof or d not in xof:
            continue
        regions[(t["transcript"], t["start"], t["end"], p)].append(
            (d, t["donor_allowed"] == "True"))

    pairs, impossible = Counter(), 0
    for (_, _, _, p), cands in regions.items():
        legal = [d for d, ok in cands if ok]
        chosen = legal[0] if legal else cands[0][0]
        if not legal:
            impossible += 1
        pairs[(chosen, p)] += 1

    maxc = max(pairs.values()) if pairs else 1
    for (d, p), c in sorted(pairs.items(), key=lambda kv: kv[1]):
        x0, x1 = xof[d], xof[p]
        travel = 1 if x1 > x0 else -1
        ax.add_patch(FancyArrowPatch(
            (x0, 0), (x1, 0), connectionstyle=f"arc3,rad={-0.40 * travel}",
            arrowstyle="-|>", mutation_scale=9 + 6 * (c / maxc),
            lw=0.6 + 3.0 * (c / maxc), color=DONOR_C, alpha=0.72, zorder=3))

    # Every gene gets the same large ringed marker. Size previously encoded
    # expression, which made the donor-only genes -- the whole subject of this
    # panel -- the hardest things on it to see. The only distinction kept is the
    # one the panel is about: which gene can be rearranged at all.
    for g in order:
        rss = rssmap.get(g, {}).get("rss_state") == "rss_present"
        ax.plot([xof[g]], [0], marker="o", ms=11,
                color=SEGMENT["V"] if rss else "white",
                markeredgecolor="black", markeredgewidth=1.2, zorder=5)
    ax.axhline(0, color="black", lw=1.1, zorder=1)

    ax.set_xlim(-1.0, len(order) + 0.2)
    ax.set_ylim(-0.30, 1.35)
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ("left", "right", "top", "bottom"):
        ax.spines[s].set_visible(False)

    n_possible = len(order) - 1
    ax.set_title(f"{L(letter)}  {locus} donor → parent network\n"
                 f"{len(pairs)} of {n_possible} possible pairs",
                 fontsize=10, loc="left", x=title_x, color="black")


# ─── C: sequence evidence ────────────────────────────────────────────────────

def panel_sequences(ax, events, letter, n_distinct, n_calls, title_x=0.0):
    """n_distinct / n_calls are dataset TOTALS; only len(events) rows are drawn."""
    width = max(len(e["cols"]) for e in events)
    # Where parent and donor agree the transcript base carries no information, so
    # it is drawn in the same neutral tone as those two rows. Only the three
    # informative outcomes -- follows the donor, follows the parent, follows
    # neither -- get a colour of their own.
    colour = {"parent": INK, "donor": DONOR_C, "other": OTHER_C,
              "uninf": GREY_DARK, "gap": GREY}
    fs = 6.6 if width <= 60 else 5.6

    for k, ev in enumerate(events):
        y = -k * 1.0
        s_off = ev["start"] - ev["lo"]
        e_off = ev["end"] - ev["lo"]
        ax.add_patch(Rectangle((s_off - .5, y - 0.40), e_off - s_off + 1, 0.80,
                               facecolor=DONOR_C, alpha=0.12, edgecolor=DONOR_C,
                               lw=0.8, zorder=1))
        for j, (_, pb, db, tb, cls) in enumerate(ev["cols"]):
            ax.text(j, y + 0.22, pb, ha="center", va="center", fontsize=fs,
                    family="monospace", color=GREY_DARK, zorder=3)
            ax.text(j, y + 0.01, db, ha="center", va="center", fontsize=fs,
                    family="monospace", color=GREY_DARK, zorder=3)
            ax.text(j, y - 0.20, tb, ha="center", va="center", fontsize=fs,
                    family="monospace", color=colour[cls], zorder=3,
                    fontweight="bold" if cls in ("donor", "other") else "normal")
        # Gene ids are not the point of this panel and only add noise; the row
        # identity is what matters, plus which locus the event came from.
        for dy, lab in ((0.22, "parent"), (0.01, "donor")):
            ax.text(-1.4, y + dy, lab, ha="right", va="center", fontsize=6.8,
                    color=GREY_DARK)
        # Black, not INK: navy is the legend's "transcript follows the parent"
        # colour, so using it for the row label makes the label look like data.
        ax.text(-1.4, y - 0.20, "transcript", ha="right", va="center",
                fontsize=6.8, color="black", fontweight="bold")
        ax.text(-11.5, y, ev["locus"], ha="left", va="center", fontsize=8.5,
                fontweight="bold", color=LOCUS.get(ev["locus"], INK))
        ax.text(width + 1.0, y, f"×{ev['n_tx']}   m={ev['m']}",
                ha="left", va="center", fontsize=6.8, fontweight="bold",
                color=DONOR_C)

    ax.set_xlim(-12, width + 9)
    ax.set_ylim(-1.0 * len(events) + 0.45, 0.75)
    ax.axis("off")
    # Say plainly that this is a subset. The counts are dataset totals and the
    # rows are the best-supported few, so a reader could otherwise take the panel
    # for the whole catalogue.
    ax.set_title(f"{L(letter)}  {len(events)} of {n_distinct} distinct tracts "
                 f"(from {n_calls} transcript-level calls), highest support first\n"
                 "shaded = called tract · ×N = transcripts carrying it · "
                 "m = donor-diagnostic positions supporting it",
                 fontsize=10, loc="left", x=title_x, color="black")


# ─── D: AID spectrum ─────────────────────────────────────────────────────────

def panel_aid(ax, stats, letter):
    cases = [("outside", "hotspot"), ("outside", "coldspot"),
             ("inside", "hotspot"), ("inside", "coldspot")]
    loci = sorted(stats)
    w = 0.38

    def mark(pv):
        return ("***" if pv <= 0.005 else "**" if pv <= 0.01 else
                "*" if pv <= 0.05 else "n.s.")

    all_ps = [float(stats[l][c][f"{m}_p_strict"])
              for l in loci for c, m in cases]
    for li, locus in enumerate(loci):
        d = stats[locus]
        vals, ps = [], []
        for cls, motif in cases:
            vals.append(float(d[cls][f"{motif}_enrichment_strict"]))
            ps.append(float(d[cls][f"{motif}_p_strict"]))
        xs = np.arange(len(cases)) + (li - (len(loci) - 1) / 2) * w
        bars = ax.bar(xs, vals, w * 0.92, color=LOCUS.get(locus, INK),
                      edgecolor="black", lw=.6, label=locus)
        # Coldspots are the same measurement on the same locus, so they read
        # better as a lighter weight of the same colour than as a second texture.
        for b, (cls, motif) in zip(bars, cases):
            if motif == "coldspot":
                b.set_alpha(0.42)
        for x, v, p in zip(xs, vals, ps):
            star = mark(p)
            ax.text(x, v * (1.07 if v >= 1 else 0.93), star, ha="center",
                    va="bottom" if v >= 1 else "top", fontsize=7.5,
                    fontweight="bold")

    ax.axhline(1.0, color="black", lw=1.2)
    ax.set_yscale("log")
    ax.set_yticks([0.25, 0.5, 1, 2, 4])
    ax.set_yticklabels(["0.25×", "0.5×", "1×", "2×", "4×"], fontsize=8)
    ax.set_ylim(0.22, 6.5)
    ax.set_xticks(range(len(cases)))
    ax.set_xticklabels(["hotspot\nOUTSIDE", "coldspot\nOUTSIDE",
                        "hotspot\nINSIDE", "coldspot\nINSIDE"], fontsize=7.5)
    ax.set_ylabel("observed / expected", fontsize=8.5)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    handles = [Patch(facecolor=LOCUS.get(l, INK), edgecolor="black", label=l)
               for l in loci]
    handles += [Patch(facecolor=GREY_DARK, edgecolor="black", label="hotspot"),
                Patch(facecolor=GREY_DARK, edgecolor="black", alpha=0.42,
                      label="coldspot")]
    ax.legend(handles=handles, fontsize=6.6, ncol=2, loc="upper right",
              framealpha=0.93)
    # Spell out only the levels this figure actually uses -- a key describing
    # "*" when nothing is starred once sends the reader hunting for it.
    used, seen = [], set()
    for lv, sym in ((0.005, "***"), (0.01, "**"), (0.05, "*")):
        for pv in all_ps:
            if sym == mark(pv) and sym not in seen:
                seen.add(sym)
                used.append(f"{sym} p≤{lv:g}")
    if any(mark(pv) == "n.s." for pv in all_ps):
        used.append("n.s.")
    ax.set_title(f"{L(letter)}  SHM leaves AID's footprint · "
                 "gene conversion does not\n"
                 "tract-restricted null · " + " · ".join(used),
                 fontsize=9.5, loc="left", x=-0.135, color="black")


# ─── assembly ────────────────────────────────────────────────────────────────

def build_events(paf, vgenes, parent_of, tracts, flank, n, locus):
    sig = [t for t in tracts if t.get("significant") == "True"
           and t.get("primary_donor", "True") == "True"]
    groups = defaultdict(list)
    for t in sig:
        groups[(t["parent"], t["donor"], t["start"], t["end"])].append(t)
    reps = []
    for members in groups.values():
        r = dict(max(members, key=lambda t: int(t["n_support"])))
        r["n_transcripts"] = len(members)
        reps.append(r)
    n_distinct, n_calls = len(reps), len(sig)
    reps.sort(key=lambda t: (-int(t["n_support"]), -t["n_transcripts"]))
    reps = reps[:n]

    want = {t["transcript"] for t in reps}
    proj_of = {}
    for rec in parse_paf(paf):
        if rec.query not in want:
            continue
        p = parent_of.get(rec.query)
        if p is None or rec.target != p:
            continue
        pr = projected_query(rec, vgenes[p])
        if pr is None:
            continue
        prev = proj_of.get(rec.query)
        if prev is None or sum(x is not None for x in pr) > sum(x is not None for x in prev):
            proj_of[rec.query] = pr

    aligner = Align.PairwiseAligner(mode="global", open_gap_score=-10,
                                    extend_gap_score=-0.5, match_score=2,
                                    mismatch_score=-1)
    cache, events = {}, []
    for t in reps:
        tid, par, don = t["transcript"], t["parent"], t["donor"]
        proj = proj_of.get(tid)
        if proj is None:
            continue
        if (par, don) not in cache:
            aln = aligner.align(vgenes[par], vgenes[don])[0]
            dp = [None] * len(vgenes[par])
            for (ps, pe), (ds, de) in zip(aln.aligned[0], aln.aligned[1]):
                for off in range(pe - ps):
                    dp[ps + off] = vgenes[don][ds + off]
            cache[(par, don)] = dp
        dproj = cache[(par, don)]
        s, e = int(t["start"]), int(t["end"])
        lo = max(0, s - flank)
        hi = min(len(vgenes[par]) - 1, e + flank)
        cols = []
        for i in range(lo, hi + 1):
            pb, db, tb = vgenes[par][i], dproj[i], proj[i]
            if tb is None or tb == "-":
                cls = "gap"
            elif db is None or db == pb:
                cls = "uninf"
            elif tb == db:
                cls = "donor"
            elif tb == pb:
                cls = "parent"
            else:
                cls = "other"
            cols.append((i, pb, db or ".", tb or "-", cls))
        events.append({"parent": par, "donor": don, "start": s, "end": e,
                       "m": int(t["n_support"]), "n_tx": t["n_transcripts"],
                       "cols": cols, "lo": lo, "locus": locus})
    return events, n_distinct, n_calls


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rss-annotation", required=True)
    ap.add_argument("--usage-assignments", required=True)
    ap.add_argument("--assignments", required=True)
    ap.add_argument("--donor-pool", required=True, help="IGL donor pool")
    ap.add_argument("--tracts", required=True, help="IGL tracts")
    ap.add_argument("--paf", required=True, help="IGL detailed PAF")
    ap.add_argument("--vgene-fasta", required=True, help="IGL V genes")
    ap.add_argument("--aid-stats", nargs="+", required=True,
                    help="<locus>=<path> aid_spectrum.tsv")
    ap.add_argument("--igh-j", type=int, required=True)
    ap.add_argument("--igl-j", type=int, required=True)
    ap.add_argument("--network-locus", default="IGL")
    ap.add_argument("--event-source", nargs="+", default=[],
                    help="Extra loci for panel E, as LOCUS=paf,vgene_fasta,tracts. "
                         "Panel E should show both loci -- a reader seeing only "
                         "IGL cannot tell whether the pattern generalises.")
    ap.add_argument("--locus-order", nargs="+", default=["IGH", "IGL"],
                    help="order loci appear in panel C; keep it the same as the "
                         "order of the architecture panels")
    ap.add_argument("--n-events", type=int, default=3,
                    help="events drawn PER locus in panel E")
    ap.add_argument("--flank", type=int, default=12)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    rss = read_tsv(args.rss_annotation)
    usage = Counter(r["best_gene"] for r in read_tsv(args.usage_assignments))
    rssmap = {r["gene"]: r for r in rss}
    genes_by_locus = defaultdict(list)
    for r in rss:
        genes_by_locus[r["locus"]].append(
            {"gene": r["gene"], "pos": int(r["pos"]), "strand": r["strand"],
             "rss": r["rss_state"] == "rss_present"})

    pool = {r["rearranged_gene"]: r for r in read_tsv(args.donor_pool)}
    tracts = read_tsv(args.tracts)
    vgenes = read_fasta(args.vgene_fasta)
    parent_of = {r["transcript"]: r["best_gene"]
                 for r in read_tsv(args.assignments)
                 if r["locus"] == args.network_locus and r["best_gene"] in vgenes}
    events, n_distinct, n_calls = build_events(
        args.paf, vgenes, parent_of, tracts, args.flank, args.n_events,
        args.network_locus)
    for spec in args.event_source:
        loc, rest = spec.split("=", 1)
        e_paf, e_fa, e_tr = rest.split(",")
        e_vg = read_fasta(e_fa)
        e_par = {r["transcript"]: r["best_gene"]
                 for r in read_tsv(args.assignments)
                 if r["locus"] == loc and r["best_gene"] in e_vg}
        ev, nd, nc = build_events(e_paf, e_vg, e_par, read_tsv(e_tr),
                                  args.flank, args.n_events, loc)
        events += ev
        n_distinct += nd
        n_calls += nc

    # Panel C lists loci in the same order as panels A and B, so the reader is
    # not asked to switch orientation halfway down the figure.
    order = {loc: i for i, loc in enumerate(args.locus_order)}
    events.sort(key=lambda e: order.get(e["locus"], len(order)))

    aid = {}
    for spec in args.aid_stats:
        loc, path = spec.split("=", 1)
        aid[loc] = kv(path)

    # One vmax across both loci, so equal stem heights mean equal counts.
    vmax = max([usage.get(g["gene"], 0)
                for gs_ in genes_by_locus.values() for g in gs_] + [1])
    ext = {loc: arch_extent(gs_, usage, vmax)
           for loc, gs_ in genes_by_locus.items()}
    # IGL needs almost no negative half (one minus-strand gene), so give it a
    # proportionally shorter row. Scaling the row height by the axis range keeps
    # data-units-per-inch identical in the two panels -- without that, sharing a
    # vmax would still leave the same count drawn at different sizes.
    hA = ext["IGH"][1] - ext["IGH"][0]
    hB = ext["IGL"][1] - ext["IGL"][0]
    rowA = 0.60
    rowB = rowA * hB / hA

    fig = plt.figure(figsize=(15.2, 12.4))
    gs = GridSpec(4, 2, figure=fig,
                  height_ratios=[rowA, rowB, 1.45, 1.30], hspace=0.62,
                  wspace=0.18, top=0.955, bottom=0.075)

    panel_architecture(fig.add_subplot(gs[0, :]), genes_by_locus["IGH"], "IGH",
                       args.igh_j, usage, "A", vmax, ext["IGH"])
    panel_architecture(fig.add_subplot(gs[1, :]), genes_by_locus["IGL"], "IGL",
                       args.igl_j, usage, "B", vmax, ext["IGL"])
    # C before D/E on purpose: those two panels talk about tracts, and about
    # differences INSIDE versus OUTSIDE them, so the reader needs to have seen
    # what a tract actually is first.
    # -0.048 on a full-width panel and -0.105 on a half-width one work out to the
    # same absolute indent, so every panel title starts at the figure margin.
    panel_sequences(fig.add_subplot(gs[2, :]), events, "C", n_distinct, n_calls,
                    title_x=-0.048)
    panel_network(fig.add_subplot(gs[3, 0]), tracts, pool, usage, rssmap,
                  args.network_locus, "D", title_x=-0.105)
    panel_aid(fig.add_subplot(gs[3, 1]), aid, "E")

    key = [
        Line2D([], [], marker="o", ls="none", ms=8, color=SEGMENT["V"],
               markeredgecolor="black", label="has an RSS — can be rearranged"),
        Line2D([], [], marker="o", ls="none", ms=5, color="white",
               markeredgecolor=GREY_DARK, markeredgewidth=0.9,
               label="silent gene — offset shows strand"),
        Line2D([], [], marker="o", ls="none", ms=8, color="white",
               markeredgecolor="black", label="no RSS (panel D)"),
        Patch(facecolor=DONOR_C, label="transcript follows the DONOR"),
        Patch(facecolor=INK, label="transcript follows the parent"),
        Patch(facecolor=OTHER_C, label="matches neither (point mutation)"),
        Patch(facecolor=GREY_DARK, label="parent = donor (uninformative)"),
    ]
    fig.legend(handles=key, loc="lower center", ncol=4, fontsize=8.5,
               frameon=True, framealpha=0.95, bbox_to_anchor=(0.5, 0.005))

    save_figure(fig, args.out)
    print(f"wrote {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
