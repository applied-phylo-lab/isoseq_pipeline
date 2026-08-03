"""
Compare BrepConvert's gene conversion calls against the topology-aware detector,
and put both through the same biological negative control.

The control
-----------
Under cis-acting deletional recombination a donor lying between the rearranged
V and J was excised and cannot have templated anything.  So for every call,
from either method, we can ask whether its donor still existed.  Random donor
choice would give an "impossible" rate equal to the deleted fraction of the
array, which the locus geometry fixes per parent gene.  A method carrying real
signal must come in below that.

Handling BrepConvert's ambiguity
--------------------------------
BrepConvert reports each event as one or more `possibility` rows, and each row's
`gene` field may list several donors separated by ';'.  That ambiguity is part
of the result, so an event is scored three ways:

  lenient  - at least one listed donor survived    (most generous reading)
  strict   - every listed donor survived
  expected - what random donor choice would give for the same parent

Comparing lenient against expected is the fair test: if even the most generous
reading cannot beat chance, the donor assignments carry no information.

BrepConvert does not report which functional allele it chose as the parent, so
the parent assignment from gc_call_functional_genes.py is used for both methods.
That keeps the comparison like-for-like.
"""
import argparse
import csv
import math
import sys
from collections import defaultdict

from scipy import stats


def read_tsv(path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(header, line.rstrip("\n").split("\t"))) for line in fh if line.strip()]


def load_map(path, key, val):
    return {r[key]: r[val] for r in read_tsv(path)}


def poisson_binomial(probs, observed):
    if not probs:
        return float("nan"), float("nan"), float("nan")
    exp = sum(probs)
    var = sum(p * (1 - p) for p in probs)
    sd = math.sqrt(var) if var > 0 else 0.0
    if sd == 0:
        return exp, sd, float("nan")
    z = (observed + 0.5 - exp) / sd
    return exp, sd, stats.norm.cdf(z)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--brepconvert-csv", required=True)
    ap.add_argument("--gene-map", required=True)
    ap.add_argument("--transcript-map", required=True)
    ap.add_argument("--assignments", required=True)
    ap.add_argument("--donor-pool", required=True)
    ap.add_argument("--tracts", required=True, help="our detector's tracts TSV")
    ap.add_argument("--locus", required=True)
    ap.add_argument("--p-threshold", type=float, default=0.05)
    ap.add_argument("--out-events", required=True)
    ap.add_argument("--out-summary", required=True)
    args = ap.parse_args()

    alias2gene = load_map(args.gene_map, "alias", "gene")
    alias2tx = load_map(args.transcript_map, "alias", "transcript")
    parent_of = {r["transcript"]: r["best_gene"] for r in read_tsv(args.assignments)
                 if r["locus"] == args.locus}

    pool = {r["rearranged_gene"]: r for r in read_tsv(args.donor_pool)}
    allowed_of, mech_of, p_imp = {}, {}, {}
    for g, r in pool.items():
        allowed_of[g] = (set() if r["allowed_donors"] == "NONE"
                         else set(r["allowed_donors"].split(",")))
        mech_of[g] = r["mechanism"]
        n_a, n_d = int(r["n_allowed_donors"]), int(r["n_deleted"])
        if r["mechanism"] == "deletion" and (n_a + n_d) > 0:
            p_imp[g] = n_d / (n_a + n_d)

    # ── BrepConvert events ────────────────────────────────────────────────────
    # BrepConvert reports the functional allele it picked in the 'allele'
    # column; prefer that over our own assignment so each method is judged on
    # its own parent call. Rows with gene == NA are events where it found a
    # clustered-mismatch stretch but could not name a donor -- they are counted
    # separately and excluded from the topology test, which needs a donor.
    events = defaultdict(list)          # (tx, event) -> list of donor sets
    spans, bc_parent = {}, {}
    n_rows = n_no_donor = 0
    with open(args.brepconvert_csv) as fh:
        rdr = csv.DictReader(fh)
        if rdr.fieldnames and "SeqID" in rdr.fieldnames:
            for row in rdr:
                n_rows += 1
                tx = alias2tx.get(row["SeqID"], row["SeqID"])
                key = (tx, row["event"])
                gene_field = (row.get("gene") or "").strip()
                if gene_field in ("", "NA"):
                    n_no_donor += 1
                    donors = set()
                else:
                    donors = {alias2gene.get(a.strip(), a.strip())
                              for a in gene_field.split(";") if a.strip()}
                events[key].append(donors)
                spans.setdefault(key, (row["start"], row["end"]))
                allele = (row.get("allele") or "").strip()
                if allele and key not in bc_parent:
                    mapped = alias2gene.get(allele)
                    if mapped:
                        bc_parent[key] = mapped

    # ── our detector's significant tracts ─────────────────────────────────────
    ours = [t for t in read_tsv(args.tracts)
            if t["significant"] == "True"]
    our_sites = defaultdict(set)
    for t in ours:
        our_sites[t["transcript"]].add((int(t["start"]), int(t["end"])))

    # ── score BrepConvert events ──────────────────────────────────────────────
    rows = []
    for (tx, ev), possibilities in sorted(events.items()):
        parent = bc_parent.get((tx, ev)) or parent_of.get(tx)
        if parent is None or mech_of.get(parent) != "deletion":
            continue
        allowed = allowed_of.get(parent, set())
        all_donors = set().union(*possibilities) if possibilities else set()
        if not all_donors:
            continue
        n_allowed = sum(1 for d in all_donors if d in allowed)
        start, end = spans[(tx, ev)]
        rows.append({
            "transcript": tx, "event": ev, "parent": parent,
            "start": start, "end": end,
            "n_donors_listed": len(all_donors),
            "n_donors_allowed": n_allowed,
            "lenient_ok": n_allowed > 0,
            "strict_ok": n_allowed == len(all_donors),
            "p_impossible": p_imp.get(parent, float("nan")),
        })

    with open(args.out_events, "w") as fh:
        fh.write("transcript\tevent\tparent\tstart\tend\tn_donors_listed"
                 "\tn_donors_allowed\tlenient_ok\tstrict_ok\n")
        for r in rows:
            fh.write(f"{r['transcript']}\t{r['event']}\t{r['parent']}\t{r['start']}"
                     f"\t{r['end']}\t{r['n_donors_listed']}\t{r['n_donors_allowed']}"
                     f"\t{r['lenient_ok']}\t{r['strict_ok']}\n")

    probs = [r["p_impossible"] for r in rows if not math.isnan(r["p_impossible"])]
    obs_lenient = sum(1 for r in rows if not r["lenient_ok"])
    obs_strict = sum(1 for r in rows if not r["strict_ok"])
    exp, sd, p_len = poisson_binomial(probs, obs_lenient)

    # our detector, same control
    our_del = [t for t in ours if t["mechanism"] == "deletion"]
    our_probs = [p_imp[t["parent"]] for t in our_del if t["parent"] in p_imp]
    our_obs = sum(1 for t in our_del
                  if t["parent"] in p_imp and t["donor_allowed"] == "False")
    our_exp, our_sd, our_p = poisson_binomial(our_probs, our_obs)

    n = len(rows)
    with open(args.out_summary, "w") as fh:
        fh.write("method\tmetric\tvalue\n")
        fh.write(f"brepconvert\trows_total\t{n_rows}\n")
        fh.write(f"brepconvert\trows_without_donor\t{n_no_donor}\n")
        fh.write(f"brepconvert\tevents_total\t{len(events)}\n")
        fh.write(f"brepconvert\tevents_deletional\t{n}\n")
        if n:
            fh.write(f"brepconvert\timpossible_lenient\t{obs_lenient}\n")
            fh.write(f"brepconvert\timpossible_lenient_frac\t{obs_lenient/n:.4f}\n")
            fh.write(f"brepconvert\timpossible_strict_frac\t{obs_strict/n:.4f}\n")
            fh.write(f"brepconvert\texpected_impossible_frac\t{exp/n:.4f}\n")
            fh.write(f"brepconvert\tp_value\t{p_len:.4g}\n")
            fh.write(f"brepconvert\tmean_donors_listed\t"
                     f"{sum(r['n_donors_listed'] for r in rows)/n:.2f}\n")
        m = len(our_probs)
        fh.write(f"topology_detector\ttracts_deletional\t{m}\n")
        if m:
            fh.write(f"topology_detector\timpossible\t{our_obs}\n")
            fh.write(f"topology_detector\timpossible_frac\t{our_obs/m:.4f}\n")
            fh.write(f"topology_detector\texpected_impossible_frac\t{our_exp/m:.4f}\n")
            fh.write(f"topology_detector\tp_value\t{our_p:.4g}\n")
        # concordance
        bc_tx = {r["transcript"] for r in rows}
        our_tx = set(our_sites)
        fh.write(f"concordance\ttranscripts_brepconvert\t{len(bc_tx)}\n")
        fh.write(f"concordance\ttranscripts_detector\t{len(our_tx)}\n")
        fh.write(f"concordance\ttranscripts_both\t{len(bc_tx & our_tx)}\n")
        fh.write(f"concordance\ttranscripts_brepconvert_only\t{len(bc_tx - our_tx)}\n")
        fh.write(f"concordance\ttranscripts_detector_only\t{len(our_tx - bc_tx)}\n")

    print(f"\n{args.locus} — same negative control applied to both methods", file=sys.stderr)
    if n:
        print(f"  BrepConvert       : {n} deletional events, "
              f"{obs_lenient/n:.1%} impossible (lenient) vs {exp/n:.1%} expected, "
              f"p={p_len:.3g}", file=sys.stderr)
        print(f"                      mean donors listed per event: "
              f"{sum(r['n_donors_listed'] for r in rows)/n:.2f}", file=sys.stderr)
    else:
        print("  BrepConvert       : no deletional events to score", file=sys.stderr)
    if our_probs:
        print(f"  topology detector : {len(our_probs)} deletional tracts, "
              f"{our_obs/len(our_probs):.1%} impossible vs {our_exp/len(our_probs):.1%} "
              f"expected, p={our_p:.3g}", file=sys.stderr)


if __name__ == "__main__":
    main()
