"""
Work out, for every V gene that could be rearranged, which other V genes are
still physically present on the chromosome afterwards and could therefore act
as gene conversion donors.

The model
---------
V(D)J recombination joins a V gene to the (D)J cluster in one of two ways,
decided by the V gene's orientation relative to J:

  * V in the SAME orientation as J -> **deletional** recombination.  The DNA
    between the V and J is excised as a circle and lost.  Any V gene lying
    between the rearranged V and J is gone and cannot be a donor.  Only genes
    on the far side of the rearranged V (distal to J) survive.

  * V in the OPPOSITE orientation to J -> **inversional** recombination.  The
    intervening DNA is inverted rather than excised, so it is retained on the
    chromosome and every V gene remains available as a donor.

Only the retained copy matters: gene conversion copies from a chromosomal
template, and the excised circle is lost from the cell.

Locus geometry is taken from the coordinates encoded in the V gene names plus
the J coordinate supplied on the command line, so the same code handles a
locus whose J sits below the V array (IGL here) or above it (IGH here).

This produces a falsifiable prediction, which is the point: for a deletionally
rearranged transcript, donors from the deleted interval are impossible.  Any
that get called are measurable false positives -> see gc_topology_test.py.
"""
import argparse
import sys
from collections import Counter

from gc_lib import read_fasta, parse_gene_names


def build_pools(info, j_pos, j_strand):
    """
    For each V gene, return (mechanism, allowed_donor_set, deleted_set).

    j_pos/j_strand describe the J gene of the same locus.
    """
    genes = sorted(info.values(), key=lambda g: g.pos)
    pools = {}
    for v in genes:
        deletional = (v.strand == j_strand)
        if not deletional:
            # inversion retains everything
            allowed = {g.name for g in genes if g.name != v.name}
            deleted = set()
            mech = "inversion"
        else:
            mech = "deletion"
            if j_pos < v.pos:
                # J lies below the array: the interval (J, V) is excised, so
                # surviving genes are those ABOVE the rearranged V.
                deleted = {g.name for g in genes
                           if j_pos < g.pos < v.pos}
                allowed = {g.name for g in genes if g.pos > v.pos}
            else:
                # J lies above the array: the interval (V, J) is excised, so
                # surviving genes are those BELOW the rearranged V.
                deleted = {g.name for g in genes
                           if v.pos < g.pos < j_pos}
                allowed = {g.name for g in genes if g.pos < v.pos}
        pools[v.name] = (mech, allowed, deleted)
    return pools


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--vgene-fasta", required=True,
                   help="V gene FASTA for ONE locus")
    p.add_argument("--locus", required=True)
    p.add_argument("--j-pos", type=int, required=True,
                   help="Contig coordinate of the J gene for this locus")
    p.add_argument("--j-strand", required=True, choices=["+", "-"])
    p.add_argument("--d-anchor-pos", type=int,
                   help="Contig coordinate of the V-proximal D gene. V recombines "
                        "with the DJ complex, not with J directly, so the deleted "
                        "interval is bounded by D. Defaults to --j-pos when absent "
                        "(identical result whenever no V gene lies between D and J). "
                        "IGL has no D genes, so J is the correct boundary there.")
    p.add_argument("--out", required=True, help="Donor pool TSV")
    p.add_argument("--out-summary", required=True)
    args = p.parse_args()

    vgenes = read_fasta(args.vgene_fasta)
    info = parse_gene_names(vgenes)
    wrong = [g.name for g in info.values() if g.locus != args.locus]
    if wrong:
        raise SystemExit(f"{args.vgene_fasta} contains non-{args.locus} genes: {wrong[:3]}")

    # Mechanism still keys off J: after D->J joining the DJ complex reads in J's
    # orientation, so that is what V is compared against. The D coordinate only
    # moves the boundary of the deleted interval.
    anchor = args.d_anchor_pos if args.d_anchor_pos is not None else args.j_pos
    pools = build_pools(info, anchor, args.j_strand)

    with open(args.out, "w") as fh:
        fh.write("rearranged_gene\tlocus\tpos\tstrand\tmechanism"
                 "\tn_allowed_donors\tn_deleted\tallowed_donors\tdeleted_genes\n")
        for name in sorted(pools, key=lambda n: info[n].pos):
            mech, allowed, deleted = pools[name]
            gi = info[name]
            fh.write(f"{name}\t{gi.locus}\t{gi.pos}\t{gi.strand}\t{mech}"
                     f"\t{len(allowed)}\t{len(deleted)}"
                     f"\t{','.join(sorted(allowed)) if allowed else 'NONE'}"
                     f"\t{','.join(sorted(deleted)) if deleted else 'NONE'}\n")

    mech_counts = Counter(m for m, _, _ in pools.values())
    n = len(pools)
    avg_allowed = sum(len(a) for _, a, _ in pools.values()) / n if n else 0
    with open(args.out_summary, "w") as fh:
        fh.write("metric\tvalue\n")
        fh.write(f"locus\t{args.locus}\n")
        fh.write(f"j_pos\t{args.j_pos}\n")
        fh.write(f"j_strand\t{args.j_strand}\n")
        fh.write(f"n_v_genes\t{n}\n")
        for mech, c in sorted(mech_counts.items()):
            fh.write(f"n_{mech}\t{c}\n")
        fh.write(f"mean_allowed_donors\t{avg_allowed:.1f}\n")
        fh.write(f"mean_donor_pool_fraction\t{avg_allowed / (n - 1):.3f}\n" if n > 1 else "")

    print(
        f"{args.locus}: {n} V genes, J at {args.j_pos}({args.j_strand})\n"
        f"  deletional : {mech_counts.get('deletion', 0)} genes "
        f"(donor pool restricted)\n"
        f"  inversional: {mech_counts.get('inversion', 0)} genes "
        f"(all donors retained)\n"
        f"  mean allowed donors: {avg_allowed:.1f} of {n - 1}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
