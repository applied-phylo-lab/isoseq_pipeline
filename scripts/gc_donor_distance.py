"""
Are the donors that get used the ones closest to the parent?

The claim being tested
----------------------
Work on chicken IGL has long held that gene conversion prefers donors that are
close to the rearranged gene -- close in sequence (the most similar pseudogenes
supply most tracts) and close on the chromosome (the proximal end of the array
is used more than the distal end).  Both are testable here, because every donor
that was NOT used is as much a datapoint as every donor that was.

Why the naive correlation is worthless
--------------------------------------
Tracts are called from *informative positions*: places where parent and donor
differ.  A donor 99% identical to its parent leaves almost nothing to detect,
and cannot produce a significant tract however often it is used; a donor 85%
identical offers informative positions everywhere.  Detection power therefore
rises with distance, which is exactly the axis under test, and in the direction
opposite to the published claim.  Plotting usage against identity and reading
off a slope measures the detector, not the biology.

So every test here is conditioned on detection opportunity.  For a parent P and
donor D, opportunity is the number of places a significant tract COULD have been
called: the number of windows of m consecutive informative positions, spaced no
more than max_gap_bp apart, inside the region the transcripts actually cover
(m and max_gap_bp are the same values the detector ran with).  A donor with zero
opportunity is dropped -- it was never in the running.

Three tests, same conditioning
------------------------------
  1. Conditional logit (the primary test).  Each observed event is treated as a
     choice of one donor out of the parent's candidate pool, with the choice
     probability proportional to opportunity x exp(beta * x), where x is
     sequence divergence or log10 genomic distance.  beta = 0 is "usage tracks
     detectability and nothing else"; beta < 0 is the published claim.  The
     parent's own pool is the stratum, so parent-to-parent differences in pool
     composition, array position and transcript depth all cancel.  p from a
     likelihood-ratio test against beta = 0.

  2. Permutation.  The same null, stated as a resampling: draw each event's
     donor from its parent's pool with probability proportional to opportunity,
     and compare the observed mean divergence / mean log10 distance of used
     donors to that null.  Same conditioning, no model, no asymptotics -- which
     matters at these event counts.

  3. Spearman.  Usage count against distance across all candidate pairs, and
     against opportunity, quoted only to show what the uncorrected analysis
     would have said and how much of it is detection bias.

Unit of analysis
----------------
Transcripts are clonally related, so counting them counts clone size, not
conversion events.  The unit is the distinct event -- one (parent, donor, start,
end) -- with transcript-level counts reported alongside and never used as the
primary n.

Ambiguous tracts are excluded by default.  Where several donors explain a tract
equally well, the detector's tie-break is alphabetical, i.e. by contig position,
which would manufacture a genomic-distance signal out of nothing.  --ties keep
includes them (all tied donors, each weighted 1/k) as a sensitivity check.
"""
import argparse
import csv
import math
import sys
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from Bio import Align
from scipy import stats
from scipy.optimize import minimize

from gc_lib import read_fasta, parse_paf, projected_query, parse_gene_names
from gc_palette import save_figure, GREY, GREY_DARK, INK, NO, CLASS_B


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# ─── parent/donor geometry ───────────────────────────────────────────────────

def align_donor_to_parent(parent_seq, donor_seq, aligner):
    """Donor base aligned at each parent position, None where gapped.

    Same routine the detector uses, so 'informative position' means the same
    thing here as it does in the tract calls.
    """
    aln = aligner.align(parent_seq, donor_seq)[0]
    out = [None] * len(parent_seq)
    for (ps, pe), (ds, de) in zip(aln.aligned[0], aln.aligned[1]):
        for off in range(pe - ps):
            out[ps + off] = donor_seq[ds + off]
    return out


def identity(parent_seq, donor_proj):
    """Identity over aligned, ungapped columns -- as in gc_gene_similarity."""
    n = same = 0
    for pb, db in zip(parent_seq, donor_proj):
        if db is None:
            continue
        n += 1
        same += (pb == db)
    return same / n if n else float("nan")


def tract_windows(informative_positions, m, max_gap_bp):
    """
    Every place a significant tract could have been called.

    A significant tract needs m informative positions in a row, consecutive ones
    no further apart than max_gap_bp.  Yield each window satisfying that.
    """
    pos = informative_positions
    for i in range(len(pos) - m + 1):
        w = pos[i:i + m]
        if all(b - a <= max_gap_bp for a, b in zip(w, w[1:])):
            yield w


def opportunity(informative_positions, m, max_gap_bp):
    """Number of windows: how many chances the detector had to fire."""
    return sum(1 for _ in tract_windows(informative_positions, m, max_gap_bp))


def expected_hits(informative_positions, m, max_gap_bp, mut_freq):
    """
    The same count, weighted by where this parent's transcripts actually vary.

    Counting windows treats every informative position as equally likely to be
    hit, and that is not safe for the comparison being made here.  A donor that
    is nearly identical to its parent differs from it only in a few clusters,
    and those clusters are the hypervariable ones -- which is also where the
    transcripts differ from the parent.  So the few windows a similar donor
    offers sit exactly where a hit is most likely, and plain window counting
    would understate its detectability and manufacture the result.

    mut_freq[i] is the fraction of this parent's transcripts differing from it
    at position i, counted OUTSIDE their own called tracts.  Excluding the
    tracts matters: they are the events under test, and a profile that included
    them would say "the transcripts vary here" about variation the conversion
    itself put there, conditioning on the outcome and correcting the effect
    away.  What is left is the SHM profile, the same inside/outside split
    gc_aid_spectrum.py uses.

    Under the detector's own SHM null a difference lands on the donor's
    particular base with chance 1/3, so a window is hit with probability
    prod(mut_freq[i]/3) and the sum over windows is the expected number of
    chance tracts.  Returned as a log, since the products underflow.
    """
    tot = None
    for w in tract_windows(informative_positions, m, max_gap_bp):
        lp = sum(math.log(mut_freq[i] / 3.0) for i in w)
        tot = lp if tot is None else np.logaddexp(tot, lp)
    return float(tot) if tot is not None else float("-inf")


def offset_of(row, kind):
    """Log detectability of this donor, under the chosen model."""
    if kind == "opportunity":
        return math.log(row["opportunity"])
    return float(row["log_expected_hits"])


def unit_weight(row, unit):
    """How much one candidate pair counts, under each unit of analysis."""
    if unit == "pair":
        return 1.0 if row["n_events"] > 0 else 0.0
    if unit == "event":
        return float(row["n_events"])
    return float(row["n_calls"])


# ─── conditional logit ───────────────────────────────────────────────────────

def conditional_logit(choices):
    """
    Fit  P(donor d chosen | parent's pool) = exp(offset_d + b*x_d) / sum(pool)
    by maximum likelihood, one term per observed event.

    The offset is log detectability, so b measures what is left after detection
    bias: b = 0 means donors are used exactly in proportion to how easily a
    tract from them could have been seen.

    choices  one tuple per observation --
             (X over the pool, shape (n_pool, k); offset over the pool; weight;
              covariate row of the donor actually used; its offset)

    With k > 1 each slope is partial: the effect of one covariate holding the
    others fixed.  That is the only honest way to separate sequence similarity
    from physical position in a tandem array, where adjacent genes are recent
    duplicates and the two are correlated by construction.

    Returns (betas, ses, p_lrt per term, loglik, loglik_null).
    """
    k = choices[0][0].shape[1]

    def negll(b, cols=None):
        b = np.atleast_1d(np.asarray(b, dtype=float))
        tot = 0.0
        for X, off, w, chosen_x, chosen_off in choices:
            if cols is None:
                lin = off + X @ b
                pick = chosen_off + chosen_x @ b
            else:
                lin = off + X[:, cols] @ b
                pick = chosen_off + chosen_x[cols] @ b
            mx = lin.max()
            tot += w * (pick - (mx + math.log(np.exp(lin - mx).sum())))
        return -tot

    res = minimize(negll, x0=np.zeros(k), method="BFGS")
    b = np.atleast_1d(res.x).astype(float)
    ll = -float(res.fun)

    # observed information by central difference, not BFGS's running estimate:
    # these standard errors are what the pooled analysis weights by, so they
    # should be the curvature of the actual likelihood at the fit
    h = 1e-4 * np.maximum(np.abs(b), 1.0)
    H = np.empty((k, k))
    for i in range(k):
        for j in range(k):
            bpp, bpm, bmp, bmm = (b.copy() for _ in range(4))
            bpp[i] += h[i]; bpp[j] += h[j]
            bpm[i] += h[i]; bpm[j] -= h[j]
            bmp[i] -= h[i]; bmp[j] += h[j]
            bmm[i] -= h[i]; bmm[j] -= h[j]
            H[i, j] = (negll(bpp) - negll(bpm) - negll(bmp) + negll(bmm)) / (
                4 * h[i] * h[j])
    try:
        cov = np.linalg.inv(H)
        se = np.sqrt(np.clip(np.diag(cov), 0, None))
        se = np.where(np.diag(cov) > 0, se, np.nan)
    except np.linalg.LinAlgError:
        se = np.full(k, np.nan)

    # p per term: likelihood ratio against the model with that term dropped,
    # so a partial slope is tested holding the other covariates where they are
    ps = np.empty(k)
    for i in range(k):
        cols = [j for j in range(k) if j != i]
        if cols:
            r0 = minimize(lambda bb: negll(bb, cols), x0=np.zeros(len(cols)),
                          method="BFGS")
            ll0 = -float(r0.fun)
        else:
            ll0 = -negll(np.zeros(k))
        ps[i] = stats.chi2.sf(max(2.0 * (ll - ll0), 0.0), 1)

    ll_null = -negll(np.zeros(k))
    return b, se, ps, ll, ll_null


def permute_null(pools, observed_mean, rng, n_iter):
    """
    Draw one donor per observation from its parent's pool, with probability
    proportional to detection opportunity, and record the weighted mean of the
    covariate.  Returns (p_two_sided, null_means).

    One draw per observation, not per used pair: the null has to have the same
    number of independent choices in it as the data, or its spread is wrong.
    """
    null = np.empty(n_iter)
    for it in range(n_iter):
        vals, ws = [], []
        for x, prob, w in pools:
            j = rng.choice(len(x), p=prob)
            vals.append(x[j])
            ws.append(w)
        null[it] = np.average(vals, weights=ws)
    centred = np.abs(null - null.mean())
    p = (np.sum(centred >= abs(observed_mean - null.mean())) + 1) / (n_iter + 1)
    return p, null


# ─── pooling across datasets ─────────────────────────────────────────────────

def parse_report(path):
    """Pull (beta, se, p) per covariate and unit back out of a report file."""
    out, meta, cov = {}, {}, None
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith("[") and "]" in line:
            cov = line[1:line.index("]")]
        elif line.startswith("  unit=") and cov:
            f = dict(kv.split("=", 1) for kv in line.split() if "=" in kv)
            out[(cov, f["unit"])] = (float(f["beta"]), float(f["se"]),
                                     int(f["n"]), float(f["p_lrt"]),
                                     float(f["p_perm"]))
        elif "\t" in line and not line.startswith(" "):
            k, v = line.split("\t", 1)
            meta[k] = v
    return meta, out


def too_few(n, se, args):
    """Say why a locus is out of the pool -- the two reasons are different."""
    if not np.isfinite(se) or se <= 0:
        return "   [not pooled: slope not identified]"
    return f"   [not pooled: fewer than {args.min_events} events]"


def combine(args):
    """
    Fixed-effect pool of the per-dataset slopes.

    No single locus here has the events to settle a question like this, so the
    reason to pool is power, and the reason it is legitimate is that beta means
    the same thing everywhere: the log change in usage rate per 1% divergence,
    or per tenfold genomic distance, each already conditioned on that locus's
    own detectability.  Inverse-variance weights, plus a sign test, which
    assumes nothing about the standard errors at all -- worth having, because
    a curvature-based SE from thirteen observations is not something to lean on.
    """
    labels, entries = [], []
    for path in args.reports:
        meta, res = parse_report(path)
        labels.append(f"{meta.get('dataset', '')} {meta.get('locus', '?')}".strip())
        entries.append(res)

    # covariates in the order the first report lists them, de-duplicated
    seen = set()
    covs = [c for c, u in entries[0]
            if u == args.unit and not (c in seen or seen.add(c))]
    missing = [c for c in covs if any((c, args.unit) not in e for e in entries)]
    covs = [c for c in covs if c not in missing]
    if missing:
        print(f"note: {', '.join(missing)} not present in every report, skipped")

    lines, pooled = [], {}
    for cov in covs:
        b = np.array([e[(cov, args.unit)][0] for e in entries])
        se = np.array([e[(cov, args.unit)][1] for e in entries])
        n = np.array([e[(cov, args.unit)][2] for e in entries])
        # A locus with only a couple of events does not identify a slope --
        # its estimate runs off to whatever value separates its two donors --
        # so it is shown but not pooled.
        ok = np.isfinite(se) & (se > 0) & (n >= args.min_events)
        w = np.where(ok, 1.0 / np.where(ok, se, 1.0) ** 2, 0.0)
        if w.sum() == 0:
            sys.exit("no dataset has enough events to pool")
        bpool = float((w * b).sum() / w.sum())
        sepool = float(math.sqrt(1.0 / w.sum()))
        z = bpool / sepool
        p = 2 * stats.norm.sf(abs(z))
        q = float((w * (b - bpool) ** 2).sum())
        p_het = stats.chi2.sf(q, max(int(ok.sum()) - 1, 1))
        # the sign test uses every locus: it needs only the direction, which a
        # two-event locus still reports even though its magnitude is unusable
        n_neg = int((b < 0).sum())
        p_sign = stats.binomtest(n_neg, len(b), 0.5,
                                 alternative="greater").pvalue
        pooled[cov] = dict(b=b, se=se, n=n, ok=ok, bpool=bpool, sepool=sepool,
                           p=p, p_het=p_het, n_neg=n_neg, p_sign=p_sign)
        lines.append(f"[{cov}] unit={args.unit}")
        for lab, bi, sei, ni, oki in zip(labels, b, se, n, ok):
            lines.append(f"  {lab:<16s} n={ni:<4d} beta={bi:+.4f} se={sei:.4f} "
                         f"RR={math.exp(bi):.3f}"
                         f"{'' if oki else too_few(ni, sei, args)}")
        lines.append(f"  POOLED beta={bpool:+.4f} se={sepool:.4f} "
                     f"RR={math.exp(bpool):.3f} "
                     f"CI=[{math.exp(bpool - 1.96 * sepool):.3f},"
                     f"{math.exp(bpool + 1.96 * sepool):.3f}] p={p:.4g}")
        lines.append(f"  heterogeneity Q={q:.2f} p={p_het:.3g}")
        lines.append(f"  sign test {n_neg}/{len(b)} negative p={p_sign:.4g}")
        lines.append("")

    with open(args.out_report, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    forest(labels, pooled, covs, args)


TITLES = {
    "divergence_given_distance_kb": "per 1% divergence, at fixed kb distance",
    "distance_kb_given_divergence": "per 10 kb, at fixed sequence divergence",
    "divergence_given_rank_distance": "per 1% divergence, at fixed array rank",
    "rank_distance_given_divergence": "per intervening gene, at fixed divergence",
    "divergence": "per 1% sequence divergence from the parent",
    "distance_kb": "per 10 kb along the locus",
    "rank_distance": "per intervening V gene in the array",
    "log10_distance": "per tenfold distance along the locus",
    "toward_j": "donor on the J side of the parent",
}


def forest(labels, pooled, covs, args):
    ncol = min(len(covs), 3)
    nrow = math.ceil(len(covs) / ncol)
    fig, axes = plt.subplots(nrow, ncol, squeeze=False, sharey=True,
                             figsize=(5.3 * ncol,
                                      (0.5 * len(labels) + 2.3) * nrow))
    flat = [a for r in axes for a in r]
    for ax in flat[len(covs):]:
        ax.set_visible(False)
    for ax, cov in zip(flat, covs):
        P = pooled[cov]
        y = np.arange(len(labels))[::-1]
        clip = lambda v: np.exp(np.clip(v, -700, 700))
        rr = clip(P["b"])
        lo = clip(P["b"] - 1.96 * P["se"])
        hi = clip(P["b"] + 1.96 * P["se"])
        prr = math.exp(P["bpool"])
        plo = math.exp(P["bpool"] - 1.96 * P["sepool"])
        phi = math.exp(P["bpool"] + 1.96 * P["sepool"])

        # A dataset with two events has an interval many orders of magnitude
        # wide.  Drawn to scale it flattens every other interval into a dot, so
        # the axis is set by the point estimates and the pooled interval, and
        # anything running off the end is capped with an arrow.
        finite = np.concatenate([rr[P["ok"]], [prr, plo, phi]])
        span = max(finite.max() / finite.min(), 4.0)
        left = finite.min() / span ** 0.6
        right = finite.max() * span ** 0.6
        for yi, r, l, h in zip(y, rr, lo, hi):
            ax.hlines(yi, max(l, left), min(h, right), color=GREY_DARK,
                      lw=1.6, zorder=2)
            if l < left:
                ax.plot(left, yi, marker="<", ms=5, color=GREY_DARK, zorder=2)
            if h > right:
                ax.plot(right, yi, marker=">", ms=5, color=GREY_DARK, zorder=2)
        for yi, r, ni, oki in zip(y, rr, P["n"], P["ok"]):
            col = INK if oki else GREY
            ax.scatter([min(max(r, left), right)], [yi], s=18 + 4 * ni,
                       color=col, zorder=3, edgecolor="white", linewidths=0.6)
        if ax is flat[len(covs) - 1] or ax is flat[min(ncol, len(covs)) - 1]:
            for yi, ni, oki in zip(y, P["n"], P["ok"]):
                ax.annotate(f"n={ni}" + ("" if oki else ", not pooled"),
                            (1.0, yi), xycoords=("axes fraction", "data"),
                            xytext=(6, 0), textcoords="offset points",
                            fontsize=7.5, color=GREY_DARK, va="center")
        ax.axvline(1.0, color=GREY_DARK, ls=":", lw=1, zorder=1)
        ax.hlines([-1], max(plo, left), min(phi, right), color=NO, lw=2.4,
                  zorder=4)
        ax.scatter([prr], [-1], marker="D", s=64, color=NO, zorder=5,
                   edgecolor="white", linewidths=0.6)
        ax.set_xscale("log")
        ax.set_xlim(left, right)
        ax.set_yticks(list(y) + [-1])
        ax.set_yticklabels(labels + ["pooled"])
        ax.set_title(f"{TITLES.get(cov, cov)}\n"
                     f"pooled {prr:.3f} [{plo:.3f}, {phi:.3f}], "
                     f"p = {P['p']:.3g}\n"
                     f"sign test {P['n_neg']}/{len(labels)} below 1, "
                     f"p = {P['p_sign']:.3g}",
                     fontsize=8, loc="left", color=INK)
        ax.spines[["top", "right", "left"]].set_visible(False)
        ax.tick_params(labelsize=8.5)
        ax.set_ylim(-1.8, len(labels) - 0.3)
    # the axis means the same thing in every panel, so label it once per column
    for ax in flat[max(0, len(covs) - ncol):len(covs)]:
        ax.set_xlabel("usage rate ratio, conditioned on detectability\n"
                      "<1 = closer donors preferred (the published claim)",
                      fontsize=8.5)
    fig.suptitle("Does donor choice depend on distance from the parent?  "
                 f"One point per locus, unit = {args.unit}, "
                 "fixed-effect pool in rose",
                 fontsize=10.5, color=INK, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 1 - 0.32 / nrow))
    save_figure(fig, args.out_figure)


# ─── main ────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--combine", action="store_true",
                    help="Pool the slopes from --reports written by earlier "
                         "runs, instead of analysing one locus")
    ap.add_argument("--reports", nargs="+", help="Reports to pool under --combine")
    ap.add_argument("--dataset", default="", help="Label used in the pooled figure")
    ap.add_argument("--min-events", type=int, default=3,
                    help="Under --combine, the fewest events a locus needs "
                         "before its slope is pooled (default 3). Below that "
                         "the slope is not identified; the locus is still drawn "
                         "and still counted in the sign test")
    ap.add_argument("--tracts")
    ap.add_argument("--vgene-fasta")
    ap.add_argument("--paf")
    ap.add_argument("--assignments")
    ap.add_argument("--donor-pool")
    ap.add_argument("--locus")
    ap.add_argument("--min-support", type=int,
                    help="m used to call significance in this locus; opportunity "
                         "windows are m informative positions wide")
    ap.add_argument("--max-gap-bp", type=int, default=5)
    ap.add_argument("--pool", choices=("all", "allowed"), default="all",
                    help="Candidate donors per parent: every other V gene, or "
                         "only those the rearrangement left on the chromosome. "
                         "'allowed' is the biologically correct pool but assumes "
                         "the topology, which is unsafe in IGH")
    ap.add_argument("--ties", choices=("drop", "keep"), default="drop",
                    help="Tracts explained equally well by several donors. The "
                         "detector breaks ties alphabetically = by position, so "
                         "keeping them fabricates a distance signal")
    ap.add_argument("--j-pos", type=int,
                    help="J coordinate for this locus. Adds a directional test: "
                         "are donors on the J side of the parent preferred, "
                         "independent of how far away they are")
    ap.add_argument("--min-coverage", type=float, default=0.5,
                    help="A parent position counts as assessable if this "
                         "fraction of the parent's transcripts cover it")
    ap.add_argument("--offset", choices=("opportunity", "expected"),
                    default="expected",
                    help="Detectability model everything is conditioned on. "
                         "'opportunity' counts the windows where a tract could "
                         "have been called; 'expected' (default) weights those "
                         "windows by the parent's own per-position mutation "
                         "frequency, which is the stricter control -- see "
                         "expected_hits()")
    ap.add_argument("--unit", choices=("pair", "event", "call"), default="event",
                    help="What one observation is. 'event' (default) is one "
                         "distinct (parent, donor, start, end) tract; 'pair' "
                         "counts each donor once however often it was used, the "
                         "most conservative reading of clonal data; 'call' "
                         "counts transcripts and is reported only for "
                         "completeness -- it counts clone size. All three are "
                         "in the report whichever is chosen")
    ap.add_argument("--n-permute", type=int, default=20000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out-table", help="Per-pair table; not used by --combine")
    ap.add_argument("--out-report", required=True)
    ap.add_argument("--out-figure", required=True)
    args = ap.parse_args()

    if args.combine:
        if not args.reports:
            ap.error("--combine needs --reports")
        return combine(args)
    missing = [f"--{n.replace('_', '-')}" for n in
               ("tracts", "vgene_fasta", "paf", "assignments", "donor_pool",
                "locus", "min_support", "out_table") if getattr(args, n) is None]
    if missing:
        ap.error("missing " + ", ".join(missing))

    vgenes = read_fasta(args.vgene_fasta)
    info = parse_gene_names(vgenes)

    parent_of = {}
    for r in read_tsv(args.assignments):
        if r["locus"] == args.locus and r["best_gene"] in vgenes:
            parent_of[r["transcript"]] = r["best_gene"]

    allowed_of = {}
    for r in read_tsv(args.donor_pool):
        allowed_of[r["rearranged_gene"]] = (
            set() if r["allowed_donors"] == "NONE"
            else set(r["allowed_donors"].split(",")))

    # transcripts projected onto their parent, exactly as the detector had them
    tx_proj = {}
    for rec in parse_paf(args.paf):
        parent = parent_of.get(rec.query)
        if parent is None or rec.target != parent:
            continue
        proj = projected_query(rec, vgenes[parent])
        if proj is None:
            continue
        prev = tx_proj.get(rec.query)
        if prev is None or sum(b is not None for b in proj) > sum(b is not None for b in prev):
            tx_proj[rec.query] = proj

    by_parent = defaultdict(list)
    tx_of_parent = defaultdict(list)
    for tid, proj in tx_proj.items():
        by_parent[parent_of[tid]].append(proj)
        tx_of_parent[parent_of[tid]].append((tid, proj))

    # ── observed usage ───────────────────────────────────────────────────────
    tracts = [t for t in read_tsv(args.tracts) if t.get("significant") == "True"]
    if args.ties == "drop":
        sig = [t for t in tracts if t.get("primary_donor") == "True"
               and int(t.get("n_candidate_donors", 1) or 1) <= 1]
    else:
        sig = tracts

    # distinct events: one (parent, donor, start, end); tied donors share weight
    ev_weight = defaultdict(float)
    calls = defaultdict(int)
    region_of_event = defaultdict(set)
    for t in sig:
        region_of_event[(t["transcript"], t["start"], t["end"])].add(t["donor"])
    seen = {}
    for t in sig:
        key = (t["parent"], t["donor"], t["start"], t["end"])
        calls[(t["parent"], t["donor"])] += 1
        if key in seen:
            continue
        seen[key] = True
        # a tie is per transcript: several donors explaining the SAME tract
        k = len(region_of_event[(t["transcript"], t["start"], t["end"])])
        ev_weight[(t["parent"], t["donor"])] += 1.0 / k

    # every significant tract of a transcript, masked out of its SHM profile
    tract_spans = defaultdict(set)
    for t in tracts:
        tract_spans[t["transcript"]].add((int(t["start"]), int(t["end"])))

    parents = sorted({p for p, _ in ev_weight})
    if not parents:
        sys.exit(f"{args.locus}: no significant tracts to analyse")

    # ── per (parent, donor) geometry, identity and opportunity ───────────────
    aligner = Align.PairwiseAligner(mode="global", open_gap_score=-10,
                                    extend_gap_score=-0.5, match_score=2,
                                    mismatch_score=-1)
    # rank in the array: physical distance counted in genes rather than in bp.
    # The published claim is usually stated this way -- "the nearest pseudogenes
    # supply most tracts" -- and bp is a poor proxy for it when gene spacing is
    # uneven, which it is in both loci.
    order = sorted(info.values(), key=lambda g: g.pos)
    rank_of = {g.name: i for i, g in enumerate(order)}

    rows = []
    for parent in parents:
        pseq = vgenes[parent]
        projs = by_parent.get(parent, [])
        if not projs:
            continue
        need = max(1, math.ceil(args.min_coverage * len(projs)))
        cov = np.zeros(len(pseq), dtype=int)
        for proj in projs:
            for i, b in enumerate(proj):
                if b is not None and b != "-":
                    cov[i] += 1
        assessable = cov >= need
        # SHM profile: how often this parent's transcripts differ from it at
        # each position, ignoring their own called tracts (see expected_hits).
        # Jeffreys pseudocount so a never-mutated position is not impossible.
        diff = np.zeros(len(pseq), dtype=float)
        shm_cov = np.zeros(len(pseq), dtype=float)
        for tid, proj in tx_of_parent[parent]:
            masked = np.zeros(len(pseq), dtype=bool)
            for a, b in tract_spans.get(tid, ()):
                masked[a:b + 1] = True
            for i, base in enumerate(proj):
                if base is None or base == "-" or masked[i]:
                    continue
                shm_cov[i] += 1
                if base != pseq[i]:
                    diff[i] += 1
        mut_freq = (diff + 0.5) / (shm_cov + 1.0)
        allowed = allowed_of.get(parent, set())
        for donor in vgenes:
            if donor == parent:
                continue
            if args.pool == "allowed" and donor not in allowed:
                continue
            dproj = align_donor_to_parent(pseq, vgenes[donor], aligner)
            inf = [i for i, (pb, db) in enumerate(zip(pseq, dproj))
                   if db is not None and db != pb and assessable[i]]
            opp = opportunity(inf, args.min_support, args.max_gap_bp)
            log_exp = expected_hits(inf, args.min_support, args.max_gap_bp, mut_freq)
            gdist = abs(info[parent].pos - info[donor].pos)
            rank_gap = abs(rank_of[parent] - rank_of[donor])
            # is the donor on the J side of the parent, or the far side?
            if args.j_pos is None:
                toward_j = ""
            elif args.j_pos > info[parent].pos:
                toward_j = int(info[donor].pos > info[parent].pos)
            else:
                toward_j = int(info[donor].pos < info[parent].pos)
            rows.append({
                "locus": args.locus,
                "parent": parent, "donor": donor,
                "parent_pos": info[parent].pos, "donor_pos": info[donor].pos,
                "genomic_distance_bp": gdist,
                "rank_distance": rank_gap,
                "donor_toward_j": toward_j,
                "same_strand": info[parent].strand == info[donor].strand,
                "donor_allowed": donor in allowed,
                "identity": round(identity(pseq, dproj), 5),
                "n_informative": len(inf),
                "opportunity": opp,
                "log_expected_hits": round(log_exp, 4) if opp else "-inf",
                "n_parent_transcripts": len(projs),
                "n_events": round(ev_weight.get((parent, donor), 0.0), 4),
                "n_calls": calls.get((parent, donor), 0),
            })

    with open(args.out_table, "w") as fh:
        cols = list(rows[0].keys())
        fh.write("\t".join(cols) + "\n")
        for r in sorted(rows, key=lambda r: (r["parent"], r["donor_pos"])):
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    # ── build the strata: donors that were actually detectable ───────────────
    live = [r for r in rows if r["opportunity"] > 0]
    dropped_used = sum(1 for r in rows if r["opportunity"] == 0 and r["n_events"] > 0)
    pool_of = defaultdict(list)
    for r in live:
        pool_of[r["parent"]].append(r)

    # Physical distance gets three parameterisations rather than one, because
    # they are not interchangeable and the published claim does not say which it
    # means.  log10 assumes every tenfold matters equally, which is a strange
    # model for an array spanning one decade (blackbird IGL runs 1.9-22 kb).
    # Linear kb is the natural scale for a locus that size.  Rank distance
    # counts intervening genes, which is how "the closest pseudogenes are used
    # most" is normally meant, and is immune to uneven gene spacing.
    covariates = {
        "divergence": ("sequence divergence (1 − identity, %)",
                       lambda r: 100.0 * (1.0 - r["identity"])),
        "distance_kb": ("physical distance along the locus (per 10 kb)",
                        lambda r: r["genomic_distance_bp"] / 10000.0),
        "rank_distance": ("distance in the array (per intervening V gene)",
                          lambda r: float(r["rank_distance"])),
        "log10_distance": ("physical distance along the locus (per tenfold bp)",
                           lambda r: math.log10(max(r["genomic_distance_bp"], 1))),
    }
    if args.j_pos is not None:
        # not a distance: which SIDE of the parent the donor sits on. A locus
        # can show no distance effect and still be polarised.
        covariates["toward_j"] = ("donor lies on the J side of the parent (0/1)",
                                  lambda r: float(r["donor_toward_j"]))

    report, results = [], {}
    n_events = sum(r["n_events"] for r in live)
    n_pairs_used = sum(1 for r in live if r["n_events"] > 0)
    rng = np.random.default_rng(args.seed)

    report.append(f"dataset\t{args.dataset}")
    report.append(f"locus\t{args.locus}")
    report.append(f"pool\t{args.pool}")
    report.append(f"ties\t{args.ties}")
    report.append(f"min_support_m\t{args.min_support}")
    report.append(f"n_parents\t{len(pool_of)}")
    report.append(f"n_candidate_pairs\t{len(rows)}")
    report.append(f"n_pairs_with_opportunity\t{len(live)}")
    report.append(f"n_pairs_used\t{n_pairs_used}")
    report.append(f"n_events\t{n_events:g}")
    report.append(f"n_transcript_calls\t{sum(r['n_calls'] for r in live)}")
    report.append(f"used_pairs_dropped_zero_opportunity\t{dropped_used}")
    report.append(f"primary_unit\t{args.unit}")
    report.append(f"detectability_offset\t{args.offset}")
    blind = [r for r in rows if r["opportunity"] == 0]
    if blind:
        report.append(f"pairs_undetectable_by_construction\t{len(blind)}"
                      f"\tidentity_range\t{min(r['identity'] for r in blind):.4f}"
                      f"-{max(r['identity'] for r in blind):.4f}")

    for name, (label, fn) in covariates.items():
        # events as choices within their parent's pool
        by_unit = {}
        for unit in ("pair", "event", "call"):
            choices, pools = [], []
            for parent, pool in pool_of.items():
                x = np.array([fn(r) for r in pool])
                off = np.array([offset_of(r, args.offset) for r in pool])
                prob = np.exp(off - off.max())
                prob = prob / prob.sum()
                for r in pool:
                    total = unit_weight(r, unit)
                    if total <= 0:
                        continue
                    # split into unit-weight observations so that the model and
                    # the permutation both see as many choices as the data has
                    n_obs = max(1, int(round(total)))
                    w = total / n_obs
                    for _ in range(n_obs):
                        choices.append((x[:, None], off, w,
                                        np.array([fn(r)]),
                                        offset_of(r, args.offset)))
                        pools.append((x, prob, w))
            bs, ses, pss, _, _ = conditional_logit(choices)
            beta, se, p_lrt = float(bs[0]), float(ses[0]), float(pss[0])
            obs = np.average([float(c[3][0]) for c in choices],
                             weights=[c[2] for c in choices])
            p_perm, null = permute_null(pools, obs, rng, args.n_permute)
            by_unit[unit] = dict(beta=beta, se=se, p_lrt=p_lrt, obs=obs,
                                 p_perm=p_perm, null=null, n_obs=len(choices))

        beta = by_unit[args.unit]["beta"]
        se = by_unit[args.unit]["se"]
        p_lrt = by_unit[args.unit]["p_lrt"]
        obs = by_unit[args.unit]["obs"]
        p_perm = by_unit[args.unit]["p_perm"]
        null = by_unit[args.unit]["null"]

        allx = np.array([fn(r) for r in live])
        cnt = np.array([r["n_events"] for r in live])
        rho, p_rho = stats.spearmanr(allx, cnt)
        detect = np.array([offset_of(r, args.offset) for r in live])
        rho_opp, p_opp = stats.spearmanr(allx, detect)
        used = np.array([fn(r) for r in live if r["n_events"] > 0])
        unused = np.array([fn(r) for r in live if r["n_events"] <= 0])
        if len(used) and len(unused):
            u_stat, p_mw = stats.mannwhitneyu(used, unused, alternative="two-sided")
        else:
            p_mw = float("nan")

        results[name] = dict(label=label, beta=beta, se=se, p_lrt=p_lrt,
                             by_unit=by_unit,
                             obs=obs, null=null, p_perm=p_perm,
                             rho=rho, p_rho=p_rho, rho_opp=rho_opp, p_opp=p_opp,
                             used=used, unused=unused, p_mw=p_mw,
                             x=allx, count=cnt,
                             opp=detect,
                             detect_label=("log detection opportunity"
                                           if args.offset == "opportunity"
                                           else "log expected chance tracts"))

        report.append("")
        report.append(f"[{name}] {label}")
        report.append(f"  conditional_logit_beta\t{beta:+.4f}")
        report.append(f"  conditional_logit_se\t{se:.4f}")
        report.append(f"  rate_ratio_per_unit\t{math.exp(beta):.3f}")
        report.append(f"  p_likelihood_ratio\t{p_lrt:.4g}")
        report.append(f"  observed_weighted_mean_used\t{obs:.4f}")
        report.append(f"  null_mean_opportunity_weighted\t{null.mean():.4f}")
        report.append(f"  p_permutation\t{p_perm:.4g}")
        report.append(f"  mean_used\t{used.mean() if len(used) else float('nan'):.4f}")
        report.append(f"  mean_unused\t{unused.mean() if len(unused) else float('nan'):.4f}")
        report.append(f"  p_mannwhitney_used_vs_unused\t{p_mw:.4g}")
        report.append(f"  spearman_usage_vs_x\t{rho:+.3f}\tp\t{p_rho:.4g}")
        report.append(f"  spearman_opportunity_vs_x\t{rho_opp:+.3f}\tp\t{p_opp:.4g}")
        for unit in ("pair", "event", "call"):
            u = by_unit[unit]
            report.append(f"  unit={unit:<5s} n={u['n_obs']:<4d}"
                          f" beta={u['beta']:+.4f} se={u['se']:.4f}"
                          f" RR={math.exp(u['beta']):.3f}"
                          f" p_lrt={u['p_lrt']:.4g} p_perm={u['p_perm']:.4g}")

    # ── partial slopes ───────────────────────────────────────────────────────
    # Physical distance and sequence divergence are not independent in a tandem
    # array: neighbours are recent duplicates, so they are also similar.  A
    # marginal slope on either one carries the other inside it.  These models
    # fit both at once, so each slope is "at fixed <the other>".
    joint = [("divergence", "distance_kb"), ("divergence", "rank_distance")]
    for pair in joint:
        if any(c not in covariates for c in pair):
            continue
        choices = []
        for parent, pool in pool_of.items():
            X = np.column_stack([[covariates[c][1](r) for r in pool]
                                 for c in pair])
            off = np.array([offset_of(r, args.offset) for r in pool])
            for i, r in enumerate(pool):
                total = unit_weight(r, args.unit)
                if total <= 0:
                    continue
                n_obs = max(1, int(round(total)))
                w = total / n_obs
                for _ in range(n_obs):
                    choices.append((X, off, w, X[i], off[i]))
        bs, ses, pss, _, _ = conditional_logit(choices)
        for c, b_, s_, p_ in zip(pair, bs, ses, pss):
            other = "_".join(x for x in pair if x != c)
            report.append("")
            report.append(f"[{c}_given_{other}] {covariates[c][0]}, "
                          f"holding {other} fixed")
            report.append(f"  unit={args.unit:<5s} n={len(choices):<4d}"
                          f" beta={b_:+.4f} se={s_:.4f} RR={math.exp(b_):.3f}"
                          f" p_lrt={p_:.4g} p_perm=nan")

    with open(args.out_report, "w") as fh:
        fh.write("\n".join(report) + "\n")
    print("\n".join(report))

    make_figure(args, results, live, n_events)


# ─── figure ──────────────────────────────────────────────────────────────────

def make_figure(args, results, live, n_events):
    names = list(results)
    fig, axes = plt.subplots(len(names), 3, figsize=(13.5, 3.8 * len(names)),
                             squeeze=False)

    for row, name in enumerate(names):
        R = results[name]
        # (1) the confound itself: opportunity against the covariate
        ax = axes[row][0]
        ax.scatter(R["x"], R["opp"], s=14, color=GREY_DARK, alpha=0.55,
                   linewidths=0, zorder=2)
        u = R["count"] > 0
        ax.scatter(R["x"][u], R["opp"][u], s=34, color=CLASS_B,
                   edgecolor="white", linewidths=0.5, zorder=3)
        ax.set_xlabel(R["label"])
        ax.set_ylabel(R["detect_label"] + "\n(what the test conditions on)")
        ax.set_title(f"the confound: ρ = {R['rho_opp']:+.2f}, "
                     f"p = {R['p_opp']:.3g}", fontsize=9, loc="left", color=INK)

        # (2) raw usage against the covariate -- what the naive analysis sees
        ax = axes[row][1]
        ax.scatter(R["x"], R["count"], s=18, color=GREY_DARK, alpha=0.55,
                   linewidths=0, zorder=2)
        ax.scatter(R["x"][u], R["count"][u], s=38, color=CLASS_B,
                   edgecolor="white", linewidths=0.5, zorder=3)
        ax.set_xlabel(R["label"])
        ax.set_ylabel("events with this donor")
        ax.set_title(f"uncorrected: ρ = {R['rho']:+.2f}, p = {R['p_rho']:.3g}",
                     fontsize=9, loc="left", color=INK)

        # (3) the test: observed against the opportunity-matched null
        ax = axes[row][2]
        ax.hist(R["null"], bins=40, color=GREY, edgecolor=GREY_DARK,
                linewidth=0.4, zorder=2)
        ax.axvline(R["obs"], color=NO, lw=2, zorder=3)
        ax.set_xlabel(f"mean {R['label']} of used donors")
        ax.set_ylabel("permutation replicates")
        rr = math.exp(R["beta"])
        ax.set_title(f"corrected: rate ratio {rr:.2f}/unit, "
                     f"p = {R['p_lrt']:.3g}\npermutation p = {R['p_perm']:.3g}",
                     fontsize=9, loc="left", color=INK)

    for ax in axes.flat:
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(labelsize=8)
        ax.xaxis.label.set_size(8.5)
        ax.yaxis.label.set_size(8.5)

    fig.suptitle(
        f"{args.locus}: is donor choice biased by distance from the parent?  "
        f"{n_events:g} distinct events, {len(live)} detectable donor–parent pairs "
        f"(unit={args.unit}, offset={args.offset}, ties={args.ties})",
        fontsize=10.5, color=INK, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 1 - 0.5 / len(names)))
    save_figure(fig, args.out_figure)


if __name__ == "__main__":
    main()
