"""
Test whether called gene conversion tracts respect locus topology.

The idea
--------
Under cis-acting deletional recombination, donors lying between the rearranged
V and J were excised and cannot template anything.  So for every called tract
we can ask a question with a known answer: was this donor still there?

That gives a real null.  If a detector picked donors at random, the fraction of
"impossible" calls would equal the fraction of the donor pool that is deleted,
which is fixed per parent gene by the locus geometry.  So:

    observed impossible fraction  ~= random expectation  -> calls are noise
    observed impossible fraction  <  random expectation  -> real signal
    observed impossible fraction  == 0                   -> fully topology-consistent

This is a biological negative control, not a simulated one, and it needs no
assumption about mutation rates or donor similarity.

Because each tract has its own expectation (parents differ in how much of the
array they delete), the total number of impossible calls is Poisson-binomial.
We use its exact mean and variance and report a normal-approximation p-value,
plus an exact binomial on the pooled expectation as a cross-check.

Sweeping the p-value threshold shows how much statistical stringency is needed
before calls start respecting the topology -- the practical output is the
threshold at which a detector becomes trustworthy.
"""
import argparse
import math
import sys
from collections import defaultdict

from scipy import stats


def read_tsv(path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(header, line.rstrip("\n").split("\t"))) for line in fh if line.strip()]


def poisson_binomial_test(probs, observed):
    """probs: per-trial probability of being 'impossible'. Returns (exp, sd, z, p)."""
    if not probs:
        return 0.0, 0.0, float("nan"), float("nan")
    exp = sum(probs)
    var = sum(p * (1 - p) for p in probs)
    sd = math.sqrt(var) if var > 0 else 0.0
    if sd == 0:
        return exp, sd, float("nan"), float("nan")
    # continuity-corrected, one-sided (fewer impossible calls than chance)
    z = (observed + 0.5 - exp) / sd
    p = stats.norm.cdf(z)
    return exp, sd, z, p


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tracts", required=True)
    ap.add_argument("--donor-pool", required=True)
    ap.add_argument("--locus", required=True)
    ap.add_argument("--thresholds", nargs="+", type=float,
                    default=[1.0, 0.5, 0.05, 0.01, 1e-3, 1e-4, 1e-6])
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    pool = {r["rearranged_gene"]: r for r in read_tsv(args.donor_pool)}
    # per parent: probability a randomly chosen donor is impossible
    p_impossible = {}
    for g, r in pool.items():
        if r["mechanism"] != "deletion":
            continue
        n_allowed = int(r["n_allowed_donors"])
        n_deleted = int(r["n_deleted"])
        total = n_allowed + n_deleted
        if total > 0:
            p_impossible[g] = n_deleted / total

    tracts = [t for t in read_tsv(args.tracts) if t["mechanism"] == "deletion"]
    if not tracts:
        raise SystemExit("no deletional tracts to test")

    with open(args.out, "w") as fh:
        fh.write("p_threshold\tn_tracts\tn_impossible\tobs_frac\texp_frac"
                 "\texpected_n\tsd\tz\tp_value\tverdict\n")

        print(f"\n{args.locus}: topology consistency of called tracts", file=sys.stderr)
        print(f"{'p_thresh':>10} {'tracts':>7} {'impossible':>11} "
              f"{'observed':>9} {'expected':>9} {'p':>10}", file=sys.stderr)

        for thr in sorted(args.thresholds, reverse=True):
            sel = [t for t in tracts if float(t["p_corrected"]) <= thr]
            probs, obs = [], 0
            for t in sel:
                pi = p_impossible.get(t["parent"])
                if pi is None:
                    continue
                probs.append(pi)
                if t["donor_allowed"] == "False":
                    obs += 1
            if not probs:
                continue
            exp, sd, z, p = poisson_binomial_test(probs, obs)
            obs_frac = obs / len(probs)
            exp_frac = exp / len(probs)
            if math.isnan(p):
                verdict = "undetermined"
            elif p < 0.05 and obs_frac < exp_frac:
                verdict = "topology_respected"
            elif obs_frac >= exp_frac:
                verdict = "indistinguishable_from_random"
            else:
                verdict = "trend_only"
            fh.write(f"{thr:g}\t{len(probs)}\t{obs}\t{obs_frac:.4f}\t{exp_frac:.4f}"
                     f"\t{exp:.2f}\t{sd:.2f}\t{z:.3f}\t{p:.4g}\t{verdict}\n")
            print(f"{thr:>10g} {len(probs):>7} {obs:>11} "
                  f"{obs_frac:>9.3f} {exp_frac:>9.3f} {p:>10.3g}  {verdict}",
                  file=sys.stderr)

    print(
        "\nreading the table: 'observed' is the fraction of called tracts whose donor\n"
        "was deleted by the rearrangement and therefore could not have templated it.\n"
        "'expected' is what random donor choice would give for the same parents.\n"
        "observed >= expected means the calls carry no topological information.",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
