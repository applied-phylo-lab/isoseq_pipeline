"""
Donor usage drawn along the locus, so the distance question can be read off the
picture instead of off a rate ratio.

Why this exists
---------------
gc_donor_distance.py answers "are closer donors used more often" with a slope
and a confidence interval.  That is the right way to test it and the wrong way
to show it.  These are the same numbers arranged so the array itself is the
x-axis: the parent sits at a marked position and the question becomes whether
the bars pile up near it.

Why raw usage bars would be a lie
---------------------------------
A tract is called from positions where parent and donor differ, so a donor that
resembles its parent is intrinsically harder to detect.  Detection power varies
by orders of magnitude across an array, and it varies along the very axis being
plotted.  A bar chart of raw counts therefore shows the detector's reach as much
as the biology.

Every panel here plots observed usage against the usage that detectability alone
predicts.  Expected counts come from the same model the test conditions on: each
donor's share of the parent's total detection opportunity, weighted by the
parent's own mutation profile taken from outside called tracts.  Bars above the
grey outline are donors used more than their detectability explains; bars below
it are donors used less.  Genes with no detection opportunity at all are drawn
as open marks on the baseline, because "never assessable" and "assessed and
never used" are different statements and must not look the same.

Two layouts
-----------
  one parent   the array as a track, one bar per V gene at its true coordinate
  many parents a donor x parent heatmap, both axes ordered along the locus, with
               each parent's own position marked in its column -- so physical
               distance is vertical distance from that mark

Both are followed by the distance profile: candidate pairs binned into equal
slices of total detection opportunity, ordered by distance from the parent.
Equal-opportunity binning makes the null flat, so any real distance preference
shows up as bars departing from a horizontal line.
"""
import argparse
import csv
import math
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.colors
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from gc_palette import save_figure, GREY, GREY_DARK, INK, NO, CLASS_B, SEGMENT


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def short(name):
    parts = name.split(".")
    return parts[1] if len(parts) > 1 else name


def read_headline(path, unit):
    """
    Pull the physical-distance slope out of the report that produced it.

    The partial slope is the one to quote -- distance at fixed sequence
    divergence -- because the two are correlated in a tandem array.
    """
    want, cov = f"unit={unit}", None
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith("[") and "]" in line:
            cov = line[1:line.index("]")]
        elif cov == "distance_kb_given_divergence" and want in line:
            f = dict(kv.split("=", 1) for kv in line.split() if "=" in kv)
            return (f"{float(f['RR']):.2f}× per 10 kb at fixed divergence "
                    f"(p = {float(f['p_lrt']):.3g})")
    return ""


def expected_counts(rows, total_events):
    """
    Split a parent's observed events across its donors in proportion to
    detection opportunity.  This is the null the test uses, made concrete:
    what the bar chart would look like if donor choice were nothing but
    detectability.
    """
    logs = np.array([float(r["log_expected_hits"])
                     if r["log_expected_hits"] != "-inf" else -np.inf
                     for r in rows])
    if not np.isfinite(logs).any():
        return np.zeros(len(rows))
    w = np.exp(logs - np.nanmax(logs[np.isfinite(logs)]))
    w[~np.isfinite(logs)] = 0.0
    return total_events * w / w.sum() if w.sum() else np.zeros(len(rows))


# ─── layout 1: the array as a track ──────────────────────────────────────────

def draw_track(ax, rows, parent, j_pos):
    """One bar per V gene at its true coordinate; grey outline = expected."""
    rows = sorted(rows, key=lambda r: int(r["donor_pos"]))
    x = np.array([int(r["donor_pos"]) for r in rows]) / 1000.0
    obs = np.array([float(r["n_events"]) for r in rows])
    exp = expected_counts(rows, obs.sum())
    live = np.array([int(r["opportunity"]) > 0 for r in rows])
    ppos = int(rows[0]["parent_pos"]) / 1000.0

    spacing = np.median(np.diff(np.sort(x))) if len(x) > 1 else 1.0
    w = max(spacing * 0.55, 0.15)

    ax.bar(x[live], obs[live], width=w, color=CLASS_B, zorder=3,
           edgecolor="white", linewidth=0.4, label="observed events")
    ax.bar(x[live], exp[live], width=w, facecolor="none", edgecolor=GREY_DARK,
           linewidth=1.1, zorder=4, label="expected from detectability alone")
    if (~live).any():
        ax.plot(x[~live], np.zeros((~live).sum()), marker="v", ms=5,
                color=GREY, ls="none", zorder=5, markeredgecolor=GREY_DARK,
                markeredgewidth=0.6, label="no detection opportunity")

    ax.axvline(ppos, color=INK, lw=1.8, zorder=6)
    ax.annotate(f"parent\n{short(parent)}", (ppos, ax.get_ylim()[1]),
                xytext=(4, -2), textcoords="offset points", fontsize=7.5,
                color=INK, va="top", ha="left", fontweight="bold")
    if j_pos is not None:
        ax.axvline(j_pos / 1000.0, color=SEGMENT["J"], lw=1.8, ls="--", zorder=6)
        ax.annotate("J", (j_pos / 1000.0, ax.get_ylim()[1]), xytext=(3, -2),
                    textcoords="offset points", fontsize=8,
                    color=SEGMENT["J"], va="top", fontweight="bold")
    ax.set_xlabel("position along the contig (kb)")
    ax.set_ylabel("conversion events")
    ax.spines[["top", "right"]].set_visible(False)


# ─── layout 2: donor x parent heatmap ────────────────────────────────────────

def draw_heatmap(fig, ax, by_parent, parents):
    """
    Parents are rows, donors are columns ordered along the locus, and each
    parent's own position is marked in its row.  Physical distance is then
    horizontal distance from that mark, which is the quantity under test.

    Three cell states have to stay distinguishable, and a plain heatmap runs
    them together: used (coloured, one step per event), assessed but never used
    (grey), and never assessable because parent and donor barely differ (white
    with a cross).  The last is not a zero -- it is a missing observation -- and
    colouring it like one would overstate how much of the array was searched.
    """
    donors = sorted({(int(r["donor_pos"]), r["donor"])
                     for rows in by_parent.values() for r in rows})
    xidx = {name: i for i, (_, name) in enumerate(donors)}
    nx, ny = len(donors), len(parents)

    counts = np.zeros((ny, nx))
    assessed = np.zeros((ny, nx), dtype=bool)
    for i, p in enumerate(parents):
        for r in by_parent[p]:
            j = xidx[r["donor"]]
            if int(r["opportunity"]) > 0:
                assessed[i, j] = True
                counts[i, j] = float(r["n_events"])

    vmax = max(1, int(counts.max()))
    # counts here run 1-4, so a continuous ramp wastes its range on values that
    # never occur and leaves "used once" almost invisible against the grey
    steps = ["#EA7580", "#B25D91", "#5A4A8A", "#172869"]
    cmap = matplotlib.colors.ListedColormap(
        [steps[min(i, len(steps) - 1)] for i in range(vmax)])
    ax.set_facecolor("white")
    # assessed-but-unused first, as the ground the coloured cells sit on
    ax.imshow(np.ma.masked_where(~assessed, np.zeros_like(counts)),
              aspect="auto", origin="lower", cmap=matplotlib.colors.ListedColormap([GREY]),
              vmin=0, vmax=1, interpolation="nearest")
    im = ax.imshow(np.ma.masked_where(counts < 1, counts), aspect="auto",
                   origin="lower", cmap=cmap, vmin=0.5, vmax=vmax + 0.5,
                   interpolation="nearest")
    ys, xs = np.where(~assessed)
    ax.plot(xs, ys, marker="x", ms=3.2, color=GREY_DARK, ls="none", zorder=3,
            markeredgewidth=0.7)

    for i, p in enumerate(parents):
        ppos = int(by_parent[p][0]["parent_pos"])
        near = min(range(nx), key=lambda k: abs(donors[k][0] - ppos))
        ax.plot([near], [i], marker="D", ms=7, color=INK, zorder=6,
                markeredgecolor="white", markeredgewidth=1.0)

    ax.set_yticks(range(ny))
    ax.set_yticklabels([short(p) for p in parents], fontsize=8)
    step = max(1, nx // 18)
    ticks = list(range(0, nx, step))
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{donors[k][0] / 1000:.0f}" for k in ticks], fontsize=7.5)
    ax.set_xlabel("donor V gene, ordered along the locus (contig position, kb)")
    ax.set_ylabel("parent")
    ax.set_ylim(-0.6, ny - 0.4)

    cb = fig.colorbar(im, ax=ax, fraction=0.02, pad=0.015,
                      ticks=range(1, vmax + 1))
    cb.set_label("events", fontsize=8)
    cb.ax.tick_params(labelsize=7)
    ax.legend(handles=[
        Line2D([], [], marker="D", ms=7, color=INK, ls="none",
               markeredgecolor="white", label="the parent's own position"),
        Patch(facecolor=GREY, label="assessed, never used"),
        Line2D([], [], marker="x", ms=5, color=GREY_DARK, ls="none",
               label="too similar to assess"),
    ], fontsize=7.5, frameon=False, ncol=3, loc="upper left",
        bbox_to_anchor=(0, -0.22))
    return nx


# ─── the distance profile, which is the answer in one panel ──────────────────

def draw_profile(ax, live, nbins):
    """
    Candidate pairs sliced into equal shares of total detection opportunity,
    ordered by distance from the parent.

    Equal-opportunity slices are the point: under "usage tracks detectability
    and nothing else" every slice expects the same number of events, so the null
    is a flat line and any distance preference is a departure from it.  Binning
    by distance instead would leave a sloped null that has to be explained.
    """
    d = np.array([int(r["genomic_distance_bp"]) for r in live]) / 1000.0
    w = np.array([r["_w"] for r in live])
    obs = np.array([float(r["n_events"]) for r in live])

    order = np.argsort(d)
    d, w, obs = d[order], w[order], obs[order]
    cum = np.cumsum(w) / w.sum()
    edges = [0]
    for k in range(1, nbins):
        edges.append(int(np.searchsorted(cum, k / nbins)))
    edges.append(len(d))
    edges = sorted(set(edges))

    centres, o_bin, labels = [], [], []
    for a, b in zip(edges, edges[1:]):
        if b <= a:
            continue
        centres.append(len(centres))
        o_bin.append(obs[a:b].sum())
        lo, hi = d[a], d[b - 1]
        labels.append(f"{lo:.1f}" if abs(hi - lo) < 0.05
                      else f"{lo:.1f}–{hi:.1f}")
    expected = obs.sum() / len(o_bin)

    ax.bar(centres, o_bin, width=0.72, color=CLASS_B, zorder=3,
           edgecolor="white", linewidth=0.5)
    ax.axhline(expected, color=GREY_DARK, lw=1.6, ls="--", zorder=4,
               label="expected if only detectability mattered")
    ax.legend(fontsize=7.5, frameon=False, loc="upper left")
    ax.set_xticks(centres)
    ax.set_xticklabels(labels, fontsize=7.5)
    ax.set_xlabel("distance from the parent (kb) — equal-opportunity bins")
    ax.set_ylabel("conversion events")
    ax.spines[["top", "right"]].set_visible(False)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--table", required=True,
                    help="{locus}_donor_distance.tsv from gc_donor_distance.py")
    ap.add_argument("--locus", required=True)
    ap.add_argument("--dataset", default="")
    ap.add_argument("--j-pos", type=int)
    ap.add_argument("--bins", type=int, default=5)
    ap.add_argument("--report",
                    help="{locus}_donor_distance.txt. The headline rate ratio is "
                         "read from it rather than typed in, so the number on "
                         "the figure cannot drift from the number in the test")
    ap.add_argument("--unit", default="event")
    ap.add_argument("--out-figure", required=True)
    args = ap.parse_args()

    subtitle = read_headline(args.report, args.unit) if args.report else ""
    rows = [r for r in read_tsv(args.table) if r["locus"] == args.locus]
    by_parent = defaultdict(list)
    for r in rows:
        by_parent[r["parent"]].append(r)
    parents = sorted(by_parent, key=lambda p: int(by_parent[p][0]["parent_pos"]))

    # weight each pair by its detectability, for the equal-opportunity binning
    live = []
    for p in parents:
        pr = [r for r in by_parent[p] if int(r["opportunity"]) > 0]
        if not pr:
            continue
        w = expected_counts(pr, 1.0)
        for r, wi in zip(pr, w):
            r["_w"] = wi
            live.append(r)

    n_events = sum(float(r["n_events"]) for r in rows)
    single = len(parents) == 1

    if single:
        fig, axes = plt.subplots(2, 1, figsize=(11.5, 7.2),
                                 gridspec_kw={"height_ratios": [1.35, 1]})
        draw_track(axes[0], by_parent[parents[0]], parents[0], args.j_pos)
        axes[0].legend(fontsize=7.5, frameon=False, loc="upper left")
        axes[0].set_title(
            f"{args.dataset} {args.locus} — every V gene at its coordinate, "
            f"bar height = times used as a donor",
            fontsize=10, loc="left", color=INK)
    else:
        fig = plt.figure(figsize=(13.0, 8.2))
        gs = fig.add_gridspec(2, 1, height_ratios=[1.15, 1], hspace=0.55)
        ax0 = fig.add_subplot(gs[0])
        draw_heatmap(fig, ax0, by_parent, parents)
        ax0.set_title(
            f"{args.dataset} {args.locus} — who donates to whom, laid out along "
            f"the locus\nclose donors would sit next to the navy diamond",
            fontsize=10, loc="left", color=INK)
        axes = [ax0, fig.add_subplot(gs[1])]

    draw_profile(axes[1], live, args.bins)
    axes[1].set_title(
        f"{n_events:g} events over {len(live)} detectable donor–parent pairs"
        + (f" · {subtitle}" if subtitle else ""),
        fontsize=9.5, loc="left", color=INK)

    fig.tight_layout()
    save_figure(fig, args.out_figure)


if __name__ == "__main__":
    main()
