# Figure caption — `FIG_main_gene_conversion`

**Figure X. Immunoglobulin diversification in the red-winged blackbird proceeds by
gene conversion.**

Full-length PacBio Kinnex IsoSeq transcripts (678 immunoglobulin transcripts: 109
IGH, 567 IGL) were assigned to V genes annotated in the *same individual's* genome
assembly, so no difference reported here can arise from germline mismatch between
the sequenced animal and the reference.

**(A, B)** Architecture of the IGH and IGL V arrays. Each gene is drawn as a stem
at its contig position; stem height is the number of transcripts for which that
gene is the best match, on a log scale **shared between the two panels**, so equal
heights mean equal counts and the two loci can be compared directly. Stems point
up for genes on the + strand and down for the − strand. A filled navy dot marks a
recombination signal sequence (RSS); genes with no transcripts are drawn as open
circles on the zero line, offset to their own strand. Only 25 of 162 IGH and 2 of
23 IGL V genes carry an RSS and can therefore be rearranged, yet 35 and 15
respectively are the best match for at least one transcript. Because a gene
without an RSS cannot itself be rearranged, those additional genes are either
unannotated RSSs or — as panel C indicates — donors toward which transcripts have
drifted through conversion. Teal diamond, J gene; the IGH axis is broken (`//`)
because its J lies ~140 kb beyond the end of the V array.

**(C)** Sequence-level evidence for individual conversion tracts. **Six of the 39
distinct tracts are shown** — the three best-supported per locus, ranked by the
number of donor-diagnostic positions; the panel is an illustrative subset, not
the full catalogue. For each event three rows are shown: the parent
(the RSS-bearing gene the transcript rearranged from), the donor, and the
transcript. Transcript bases are coloured by what they agree with — **teal**, the
donor and not the parent (the diagnostic positions); **navy**, the parent;
**rose**, neither, i.e. an independent point mutation; **grey**, positions where
parent and donor are identical and which therefore carry no information. Shading
marks the called tract. ×N, number of transcripts carrying that exact tract; m,
number of donor-diagnostic positions supporting it. Across both loci, 39 distinct tracts were recovered from 150
transcript-level calls in 134 transcripts (IGH 19 tracts / 28 calls; IGL 20 / 122);
those totals, not the six rows displayed, are what the panel heading counts. Tracts
require supporting positions contiguous within 5 bp and m ≥ 5 (IGL) or m ≥ 6
(IGH); both thresholds were set empirically (Methods).

**(D)** Donor → parent relationships in IGL. Genes are ordered by contig position;
arc width is proportional to the number of supporting tracts. Filled navy, RSS
present; open, RSS absent. Ten of the 22 possible donor → parent pairs are
observed, all converging on the single functional gene at the J-proximal end of
the array. None of the 122 IGL tracts uses a donor that the rearrangement had
already removed; note, however, that because the functional gene lies at the
J-proximal end no donor is ever deleted in IGL, so this check has no
discriminating power in this locus (it does in IGH, where 8 of 28 tracts have no
topologically possible donor).

**(E)** Mutation spectrum inside versus outside called tracts, both loci. Bars
give the observed fraction of differences falling at AID hotspot (WRCY/RGYW) or
coldspot (SYC/GRS) motifs, divided by the fraction expected under 1,000
permutations in which each transcript's differences are redistributed at random
among positions **of the same class in that same transcript**, each difference counted once per
clone per tract, clones being defined by shared VDJ junction — a tract-internal
difference stays within its own tract window, so the comparison is against the
local base composition rather than that of the whole gene. Coldspot bars are drawn
in a lighter tone of the locus colour. Outside tracts, differences are enriched at
hotspots and depleted at coldspots in both loci (IGH ×2.23 and ×0.51; IGL ×2.24
and ×0.40; all p = 0.002), the expected signature of somatic hypermutation, in
which the mutation occurs at the AID-targeted base itself. Inside tracts the
pattern is absent or reversed (IGL ×0.55 and ×1.47, both p = 0.002; IGH ×0.78,
n.s., and ×1.68, p = 0.014), as predicted if those differences were copied from a
donor and therefore lie wherever the donor happened to differ from the parent,
independently of where AID acted. Both processes are AID-initiated; what differs
is whether the lesion is resolved by error-prone repair at that base (hypermutation)
or by templated repair from a donor (conversion). p values are two-sided empirical
values; p = 0.002 is the resolution floor for 1,000 permutations.
`***` p ≤ 0.005; n.s., not significant.

---

## Short version

**Figure X. Gene conversion drives immunoglobulin diversification in the
red-winged blackbird.** (**A**, **B**) IGH and IGL V arrays; stem height is
transcript count (log, shared scale), direction is strand, navy dots mark an RSS,
open circles on the zero line mark silent genes. Only 25/162 IGH and 2/23 IGL V
genes carry an RSS. Teal diamond, J; the IGH axis is broken. (**C**) Sequence
evidence, six of 39 distinct tracts (three best-supported per locus): parent, donor and transcript, with transcript
bases coloured teal where they follow the donor, navy the parent, rose neither and
grey where parent and donor agree. Shading, called tract; ×N, transcripts carrying
it; m, donor-diagnostic positions. (**D**) IGL donor → parent arcs, width ∝
supporting tracts; 10 of 22 possible pairs, all converging on the single
functional gene. (**E**) Differences outside tracts are hotspot-enriched and
coldspot-depleted in both loci (all p = 0.002), the hypermutation signature;
differences inside tracts are not, consistent with templated repair from a donor.
Enrichments are relative to 1,000 tract-restricted permutations.
