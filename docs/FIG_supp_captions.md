# Supplementary figure captions

All gene counts are **distinct genomic loci**, after collapsing V gene
annotations that are the same DNA entered twice (once per strand). The raw
annotations list 162 IGH and 23 IGL V genes for the blackbird and 60 IGH / 51 IGL
for the duck; collapsing gives 103, 22, 34 and 50. No tract call changes in the
blackbird; in the duck IGH, 26 duplicate donors leave the pool and the call count
falls from 104 to 81 without altering the impossible-donor fraction.

---

## `SUPP_IGH_donor_network`

**Figure S_. Donor → parent relationships in blackbird IGH.**

Only the 31 of 103 IGH V loci that appear as a donor or a parent are drawn, at
their contig positions; 12 of the 31 carry an RSS (navy fill; grey, no RSS).
Marker direction gives orientation relative to J. One arc per event, all the same
width: **teal**, the best-matching donor is one the rearrangement left intact;
**grey**, the best-matching donor was deleted but another donor could have
supplied the tract, and the arc is drawn from that surviving donor; **rose**,
below the axis, every candidate donor had been deleted. 28 events over 28
donor → parent pairs, of which **8 have no possible donor**.

*The J position shown (219,000) is not observed.* No IGH J is present in the
bAgePho2 assembly at all: a J consensus built from 110 transcripts — carrying the
canonical C-terminus `WGSGTTVTVASG` — has no match in either haplotype, while the
same search recovers the IGL J at 100% identity over its full length. The IGH V
array, and the constant-region genes, also sit on three different contigs. The
above/below split therefore carries no topological signal for IGH and is shown
only to display the false-positive load; testing both J strands gives a
false-discovery rate of ≈1 either way.

---

## `FIG_main_gene_conversion` (tufted duck)

**Figure S_. The same analysis in the tufted duck.**

Panels as in Figure X, for 540 IG transcripts (370 IGH, 170 IGL) scored against
the *bAytFul2* assembly. **(A, B)** V array architecture: 5 of 34 IGH and 15 of 50
IGL loci carry an RSS; 18 and 11 respectively are the best match for ≥1
transcript. **(C)** Six example tracts across the whole V gene, of 35 distinct
tracts from 109 transcript-level calls; tracts require m ≥ 4 here, recalibrated by
permutation for this dataset. **(D, E)** Donor → parent networks for IGL (7 of 49
possible pairs, 4 of 28 calls impossible) and IGH (18 of 33 pairs, **45/45 events
with no possible donor**). **(F)** Mutation spectrum against the tract-restricted
null. The direction seen in the blackbird reproduces but weakly: outside-tract
differences are hotspot-enriched and coldspot-depleted in IGH (×1.49 and ×0.81,
both p = 0.002) but not in IGL (×1.12, p = 0.20; ×0.88, p = 0.038), and
inside-tract differences show no hotspot signal in either locus (×0.78 and ×0.83,
n.s.).

*Two caveats limit this dataset.* The Iso-Seq RNA derives from up to ten birds
across six tissues while the assembly is one female from that cohort, so
between-individual germline polymorphism adds hotspot-neutral differences to the
outside-tract class. And every distinct IGH tract is private to a single
sequencing run, so panel E is not interpretable as conversion — under somatic
conversion onto a shared germline the donor → parent *pairs* would recur across
tissues.

---

## `SUPP_tract_call_peaks`

**Figure S_. Conversion tracts concentrate in CDR2.**

Tracts covering each position of the parent V gene. **Tract calls** (filled) count
one per transcript-level event, using the best-supported donor where a donor is
ambiguous; **distinct tracts** (line) count each tract position once, so clonally
related transcripts collapse. Teal blocks mark CDR1 and CDR2.

**(A) IGH**, seven parent genes pooled and drawn relative to each gene's own CDR2
start — the genes fall into two length classes whose CDR2 boundaries differ by
~38 bp, so an absolute trace would be bimodal for that reason alone. FR2 is 39 bp
in all seven, so the CDR blocks are exact rather than averaged. Coverage is
confined to −49…+30 and peaks at **+15 with 21 calls from 14 distinct tracts**,
entirely within CDR2. **(B) IGL**, the single parent gene 6336602 in its own
coordinates. The dominant peak fills CDR2 (**90 calls, 6 distinct**), with smaller
peaks in CDR1 (71–83) and in FR3 (185–202) — the one substantial IGL window
outside a CDR.

CDRs lie mid-gene and so do tracts, so the two are confounded. Controlling for
position leaves IGH significant (88.7% of tract positions in a CDR against a
39.4% positional baseline, p ≤ 0.002 under both a position-matched and a
gene-swap null; 2.17× and 2.04× at CDR2 = 21 and 39 bp). IGL cannot be tested —
four windows in one gene, no swap null with a single parent, positional
p = 0.23–0.46 — so panel B is illustrative, not evidential.

Only the framework boundaries are measured, from protein-level motifs; CDR widths
are canonical (CDR2 30 bp, CDR1 24 bp).

---

## `SUPP_donor_distance_pooled`

**Figure S_. Donor choice tracks sequence similarity weakly and physical
position in the array in the direction opposite to the published claim.**

Conditional-logit slope for donor usage as a usage rate ratio, one panel per
covariate, one row per locus; marker area proportional to distinct conversion
events, whiskers 95% CI arrow-capped where they leave the axis, rose diamond the
fixed-effect pool. Below 1 is the direction reported for chicken IGL. Top row,
marginal slopes; bottom row, partial slopes holding the other covariate fixed —
those are the ones to quote, because in a tandem array sequence similarity and
physical proximity are correlated by construction (here negatively, ρ = −0.31 in
both IGL loci: more distant donors are *more* similar).

Every estimate is conditioned on detection opportunity, which is not optional: a
tract is called from positions where parent and donor differ, so a donor similar
to its parent is intrinsically harder to detect and the raw correlation runs
opposite to the biology (ρ = +0.89 between opportunity and divergence in
blackbird IGL). Opportunity is the number of windows of m informative positions
in which a significant tract could have been called, weighted by the parent's
mutation profile taken from differences **outside** called tracts. Under
unweighted window counting the divergence effect in blackbird IGL is RR 0.78
(p = 2×10⁻⁸); weighted it is RR 0.90 (p = 0.11), so most of the apparent
sequence-proximity effect is detection bias.

**Sequence proximity** (per 1% divergence, at fixed physical distance): all four
slopes negative, 4/4 sign test p = 0.063, pooled RR 0.967 [0.932, 1.003],
p = 0.069. Consistent with the published claim; not on its own evidence for it.

**Physical proximity** (per 10 kb, at fixed divergence): pooled RR **1.487
[1.107, 1.999], p = 0.0085**, all four loci in the same direction, no
heterogeneity (Q = 0.14, p = 0.93); per intervening V gene, RR 1.028
[1.006, 1.050], p = 0.013. Donors *further* along the locus are used more. The
effect strengthens rather than weakens under adjustment for sequence similarity,
holds under both detectability models, and is not produced by censoring — the
near-identical pairs excluded for zero detection opportunity are not
systematically the nearby ones in the three pooled loci (p = 0.90, 0.25, none),
though they are in duck IGL, which is not pooled.

**Directionality** (`toward_j`) is not identified in blackbird IGL or duck IGH:
each has one tract-bearing parent at the end of its array, so every donor lies on
the same side and there is no contrast. Where fitted it is flat (RR 1.23,
p = 0.79).

The unit is the distinct conversion event — one (parent, donor, start, end) —
because transcripts within a clone are not independent; counting transcripts
instead gives p = 9×10⁻⁴⁴ in blackbird IGL from 19 independent events. Tracts
explained equally well by several donors are excluded, since the detector breaks
those ties by gene name, which is gene position. Duck IGL contributes two events:
its direction counts in the sign test but its slope is not identified and it is
left out of the pool. The analysis is silent about the very closest donors — pairs
with no detection opportunity at all, up to 98.3% identity — which is the class
the chicken claim concerns.

---

## `SUPP_human_control`

**Figure S_. The detector on human repertoires, which have no gene conversion.**

Ten individuals from Rodriguez et al. (BioProject PRJNA555323), each with a
targeted IGH germline assembly and matched AIRR-seq. Every sample is scored
**against its own germline**: the assembly was annotated with IgDetective,
allelic copies collapsed at 99% identity (60–100 V genes per individual), and
50,000 merged read pairs ≥450 bp put through the identical detector.

Bars give the percentage of transcripts carrying at least one significant tract —
the human median and each bird locus. Open circles are the individual human
samples. Wilson 95% intervals for every bar are in `human_control_rates.tsv`;
they are quoted below but not drawn, since the separation between the groups is
much larger than the intervals.

Humans return a median of **0.051%** (range 0.004–0.608%, n = 9). Every bird
locus sits above every human sample: tufted duck IGL 5.8% (95% CI 2.7–12.0),
red-winged blackbird IGH 7.0% (3.1–15.5), red-winged blackbird IGL 13.0%
(10.4–16.1), tufted duck IGH 24.0% (18.0–31.4). The closest approach is tufted
duck IGL against the worst human sample, a factor of **9.5**; against the human
median it is 113×. The lower bound of the weakest bird locus (2.7%) still
exceeds the highest single human sample by 4.4×.

The unit is the percentage of transcripts, not distinct events per transcript,
because distinct-event counts saturate with sequencing depth and the libraries
differ 100–700-fold in depth (25,000–49,000 human transcripts against 71–524 in
the birds). A per-transcript rate does not have that problem: `n_support` in the
detector counts donor-supporting informative positions *within one transcript*,
so each transcript is scored independently of library size.

The comparison is loaded in the humans' favour throughout. They are run with
unconstrained parent assignment and a fully permissive donor pool — a targeted
assembly contains no J and spans many contigs, so no deletion/inversion model is
possible — and they carry more detection opportunity per transcript (median 276 bp
of covered V against 238–245 bp in the birds). The bird bars likewise use the
**any-parent** run, which matches the human setting; the RSS-restricted main run
is higher, at 20.6% (IGH) and 23.1% (IGL).

**One sample excluded.** A sample is dropped when more than half its calls come
from a single (parent, donor, start, end) combination, since it is then reporting
one event rather than a rate; the rule is applied to birds identically and no
bird locus exceeds 32.5%. Only W-79 trips it, at 94.8%: 1,764 of its 1,861 calls
are the same parent, the same donor and the same 16 bp window, on a gene
annotated non-productive — one germline allele the assembly missed, so every
transcript from that gene carries an identical mismatch block that a paralogue
supplies. Its rate (4.68%) is still below three of the four bird loci. The
excluded sample, both units and all intervals are in `human_control_rates.tsv`.
