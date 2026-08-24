# Figure caption — `FIG_main_gene_conversion`

**Figure X. Immunoglobulin diversification in the red-winged blackbird proceeds by
gene conversion.**

Full-length PacBio Kinnex IsoSeq transcripts (595 immunoglobulin transcripts: 71
IGH, 524 IGL) were assigned to V genes annotated in the *same individual's* genome
assembly, so no difference reported here can arise from germline mismatch between
the sequenced animal and the reference.

Gene counts are **distinct genomic loci**. The raw annotation lists 162 IGH and 23
IGL V genes, but 59 IGH pairs and 1 IGL pair are the same DNA entered twice, once
on each strand: their intervals overlap and their forward-strand sequence is
identical over the shared block (58 of 59 exactly). Only one member of each pair
carries the V reading frame. Collapsing them gives **103 IGH and 22 IGL** loci and
changes nothing else — no tract call, and no value in panel E, moves.

**(A, B)** Architecture of the IGH and IGL V arrays. Each gene is drawn as a stem
at its contig position; stem height is the number of transcripts for which that
gene is the best match, on a log scale **shared between the two panels**, so equal
heights mean equal counts and the two loci can be compared directly. Stems point
up for genes on the + strand and down for the − strand. A filled navy dot marks a
recombination signal sequence (RSS); genes with no transcripts are drawn as open
circles on the zero line, offset to their own strand. Only **25 of 103 IGH and 2 of
22 IGL** V genes carry an RSS and can therefore be rearranged, yet 34 and 14
respectively are the best match for at least one transcript. Because a gene
without an RSS cannot itself be rearranged, those additional genes are either
unannotated RSSs or — as panel C indicates — donors toward which transcripts have
drifted through conversion. Teal diamond, J gene; the IGH axis is broken (`//`)
because its J lies ~140 kb beyond the end of the V array.

**(C)** Six example conversion tracts shown across the **whole** V gene, three per
locus. Each row is one event: every position of the parent gene is coloured by
what the transcript matches there — **navy**, the parent; **teal**, the donor;
**rose**, neither, i.e. an independent point mutation; **pale grey**, parent and
donor are identical so the position carries no evidence either way (70–86% of a
typical gene); **darker grey**, not covered by the transcript. The teal bar above
each row marks the called tract. The pattern that matters is parent — DONOR —
parent: the transcript follows its parent on both flanks and the donor only inside
the tract, which is what distinguishes a conversion tract from a misassigned
parent. Events are ranked for display by 5′ anchoring × tract purity, not by
support: the best-supported tracts are the longest ones and leave no flank to
show. Gene identities are omitted because these are illustrations, not the result;
the totals in the heading (39 distinct tracts, 165 transcript-level calls) are.
Tracts require supporting positions contiguous within 5 bp and m ≥ 5 (IGL) or
m ≥ 6 (IGH); both thresholds were set empirically (Methods).

**(D)** Donor → parent relationships in IGL. Genes are ordered by contig position;
arc width is proportional to the number of supporting tracts. Filled navy, RSS
present; open, RSS absent. Ten of the 21 possible donor → parent pairs are
observed, all converging on the single functional gene at the J-proximal end of
the array. None of the 123 IGL tracts uses a donor that the rearrangement had
already removed. This check has no discriminating power in IGL, however, because
the functional gene lies at the J-proximal end and so no donor is ever deleted;
nor does it in IGH, where no J has been located and the deletion/inversion
assignment therefore rests on the D cluster's orientation — testing both J strands
gives a false-discovery rate of ≈1 either way, i.e. IGH donor assignment carries
no topological signal. No IGH result depending on that labelling is quoted.

**(E)** Mutation spectrum inside versus outside called tracts, both loci. Bars
give the observed fraction of differences falling at AID hotspot (WRCY/RGYW) or
coldspot (SYC/GRS) motifs, divided by the fraction expected under 1,000
permutations in which each transcript's differences are redistributed at random
among positions **of the same class in that same transcript**, each difference
counted once per clone per tract, clones being defined by shared VDJ junction
(520 IGL transcripts → 396 clones; 68 IGH → 54) — a tract-internal difference
stays within its own tract window, so the comparison is against the local base
composition rather than that of the whole gene. Coldspot bars are drawn in a
lighter tone of the locus colour. Outside tracts, differences are enriched at
hotspots and depleted at coldspots in both loci (IGH ×2.23 and ×0.51; IGL ×2.24
and ×0.40; all p = 0.002), the expected signature of somatic hypermutation, in
which the mutation occurs at the AID-targeted base itself. Inside tracts the
pattern is absent or reversed (IGL ×0.55 and ×1.47, both p = 0.002; IGH ×0.78,
n.s., and ×1.68, p = 0.014), as predicted if those differences were copied from a
donor and therefore lie wherever the donor happened to differ from the parent,
independently of where AID acted. p values are two-sided empirical values;
p = 0.002 is the resolution floor for 1,000 permutations.
`***` p ≤ 0.005; `*` p ≤ 0.05; n.s., not significant.

---

## Short version

**Figure X. Gene conversion drives immunoglobulin diversification in the
red-winged blackbird.** (**A**, **B**) IGH and IGL V arrays; stem height is
transcript count (log, shared scale), direction is strand, navy dots mark an RSS,
open circles on the zero line mark silent genes. Only 25/103 IGH and 2/22 IGL V
genes carry an RSS; counts are distinct loci after collapsing same-DNA duplicate
annotations. Teal diamond, J; the IGH axis is broken. (**C**) Six example tracts
across the whole V gene, three per locus, of 39 distinct tracts from 165
transcript-level calls: navy where the transcript follows its parent, teal the
donor, rose neither, pale grey where parent and donor agree. Teal bar, called
tract. The transcript follows its parent on both flanks and the donor only inside
the tract. (**D**) IGL donor → parent arcs, width ∝ supporting tracts; 10 of 21
possible pairs, all converging on the single functional gene. (**E**) Differences
outside tracts are hotspot-enriched and coldspot-depleted in both loci (all
p = 0.002), the hypermutation signature; differences inside tracts are not,
consistent with templated repair from a donor. Enrichments are relative to 1,000
tract-restricted permutations, counted once per clone per tract.
