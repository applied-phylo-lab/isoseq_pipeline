# Red-winged blackbird IG pipeline — what was done, and how to read the figures

This walks through the analysis end to end, in the order the data actually moves,
and then gives a panel-by-panel reading guide for every figure produced.

Everything below refers to the **matched run** unless stated otherwise:

| | |
|---|---|
| RNA | `m84094_260720_230105_s3.hifi_reads.bcM0001.bam` (PacBio Kinnex, MAS-8) |
| Germline | `bAgePho2_pri` / `bAgePho2_alt` — **the same bird that was sequenced** |
| Stage 1 results | `/local/storage/kav67/redwinged_blackbird_isoseq_matched/results/` |
| Stage 2 results | `…/results/geneconv/` |
| Configs | `config_bAgePho2_matched.yaml`, `config_geneconv_bAgePho2.yaml` |

The earlier run against the VGP assembly (`bAgePho0`, a **different** bird) lives in
`/local/storage/kav67/redwinged_blackbird_isoseq/results/` and is kept as the
contrast case, not as a result in its own right.

Every figure is written three times: `.pdf` (LaTeX), `.png` (slides, quick
viewing) and `.svg` (editing panels in Illustrator/Inkscape without re-running
anything). Multi-page reports get one SVG/PNG per page, `…_p01.svg`, `…_p02.svg`, …

---

## Stage 1 — from raw reads to IG transcripts

### 1. Split the Kinnex array (`skera`)

Kinnex concatenates ~8 cDNA molecules into one long HiFi read to use the
sequencer efficiently. `skera split` cuts them back apart at the MAS-8 adapters
(`data/mas8_primers.fasta`). One raw read in, up to 8 S-reads out.

### 2. Trim primers and orient (`lima`, `isoseq refine`)

`lima --isoseq` finds the 5′ and 3′ cDNA primers, orients every read 5′→3′ and
discards anything without both ends. `isoseq refine --require-polya` then trims
the polyA tail and drops concatemers. What survives is **FLNC** —
*Full-Length Non-Concatemer* — reads that are one complete cDNA molecule each,
from the 5′ cap to the polyA. That "full-length" guarantee is what lets us treat
a read as a real transcript rather than a fragment.

**91,291,873 FLNC reads.**

### 3. Prescreen for IG before clustering

`isoseq cluster2` on all 91.3M FLNC segfaulted after 4h34m. It also isn't what we
want: this is a whole-transcriptome library and we only care about immunoglobulin.
So FLNC are aligned against the V genes first and only the ones that hit are
clustered.

**4,767 IG FLNC reads** (≥70% identity over ≥150 bp).

### 4. Cluster into transcripts (`isoseq cluster2`)

Even with perfect chemistry, ten reads from the same mRNA differ slightly through
sequencing error. Clustering groups reads that are the same molecule and emits one
polished consensus per group. This matters enormously here: a single-read error
looks exactly like a somatic mutation, and averaging over a cluster removes almost
all of it. Cluster size also becomes the expression measure.

**709 transcripts**, of which **676 are IG** (109 IGH, 567 IGL).

### 5. Two independent IG filters, then merge

Transcripts are kept if *either* minimap2 finds a V gene hit *or* immunotools'
`vj_finder` calls them IG. The two agree on 509; minimap2 adds 167, immunotools
adds 2. **678 IG transcripts total.**

> immunotools needs both a V and a J. We could never locate an IGH J in the
> `bAgePho2` assembly — all 162 IGH V genes are on `ptg000161l`, and no J-like
> sequence with a splice donor is there — so IGH relies on the minimap2 arm.

### 6. Align transcripts to V genes

`minimap2 -x splice:hq --cs` against the V gene FASTA. The `cs` tag records the
exact edit string, which is what every downstream difference is read from. Each
transcript is projected onto V gene coordinates so that position 37 means the same
thing in every transcript.

**Identity is always scored with the last 20 bp of the V gene excluded.** V(D)J
junction formation chews back the 3′ end and adds N nucleotides, so that region is
guaranteed to mismatch no matter how functional the gene is. Including it would
make every gene look mutated.

---

## Stage 2 — gene conversion

### 7. Which genes can even be a parent? (RSS)

A V gene can only be rearranged if it has a **recombination signal sequence**.
Without an RSS it can still donate sequence into someone else's rearranged gene,
but it can never itself be the gene the transcript came from. RSS calls are taken
as-is from the curated screen in `clean_birds/gene_list.csv` — present or absent,
no re-calling, no degeneracy scoring.

| | V genes | with RSS |
|---|---|---|
| IGH | 162 | 25 |
| IGL | 23 | 2 |

Few functional genes is the *expected* picture for a gene-conversion species —
it is the situation in chicken and duck, and it is the reason gene conversion
exists: a small number of rearrangeable genes diversified from a large silent
donor array.

### 8. Which genes are actually used? (expression)

RSS is a property of the DNA. Expression is a property of the animal. A gene is
called **expression-confirmed** when ≥2 transcripts match it at ≥98% identity
over ≥200 bp (junction excluded). Genes with an RSS *and* expression are the
parents; genes with expression but no RSS are "rescued" and flagged, because they
are usually a near-identical donor being mistaken for the parent.

**IGH 19 confirmed, IGL 1 confirmed.**

### 9. Restrict the donor pool by locus topology

This is the step that makes the analysis falsifiable, and it is worth
understanding because two figures are built entirely on it.

Gene conversion here is taken to be **cis-acting** — the donor must be on the same
chromosome as the recipient. So which donors were *physically still present* when
conversion happened depends on how the gene rearranged:

- **Deletion** (V in the same orientation as J): everything between that V and
  the J is looped out and destroyed. Any donor in that interval is **gone**.
- **Inversion** (V opposite to J): the intervening DNA is flipped, not removed.
  Every donor survives.

So for each candidate parent we can compute exactly which donors were still
available. A "conversion tract" whose donor was already deleted is **impossible** —
it is a false positive, and we can count them. That count is the error rate,
measured from the data rather than assumed.

In IGL every V gene is on the same strand as J: deletional, and the one functional
gene sits at the J-proximal end, so **nothing gets deleted and every donor is
available**. In IGH the picture is mixed and both mechanisms occur.

### 10. Detect conversion tracts — two independent methods

**BrepConvert** (Fraternali lab, BLAT-based) was run first, since it is the
published tool. It required three fixes to run at all: its IMGT-style gene-name
regex returned `NA` for our names (we now emit IMGT-shaped aliases); Biostrings
≥2.77 moved `pairwiseAlignment`/`pid`/`stringDist` into `pwalign`; and `/` in
transcript IDs broke its temp FASTA writer.

**Our own detector** was written because the paper itself concedes that template
jumping and ordinary mutation can be misread as conversion. It works differently:

1. Take every position where the transcript differs from its parent.
2. Ask whether a donor gene carries exactly that base at that position
   (a *donor-diagnostic* position).
3. Require the supporting positions to be **contiguous within 5 bp** — a real
   conversion tract is a continuous stretch of donor sequence, not a scatter.
4. Keep tracts with at least **m** supporting positions, where m is set per locus
   by calibration (below).

The contiguity rule is what dropped the counts sharply (e.g. 3,646 → 116 in the
VGP IGL run). Without it the "tracts" had a median span of 53 bp carrying only 7
supporting positions — scattered agreement, not tracts.

#### How the threshold is set — and why the p-value was dropped

The original scoring was a likelihood ratio with a Bonferroni correction over the
donor pool, `p_corrected = n_donors × 3⁻ᵐ`. **That has been abandoned.** It
corrects for the wrong thing: `3⁻ᵐ` models the chance of matching a donor's base,
but says nothing about **contiguity**, which is the constraint actually doing the
work. It was therefore wildly conservative — with 162 IGH donors it demanded
m ≥ 8 and returned 10 tracts in the entire locus. Benjamini–Hochberg does not
help either: with `min_informative = 3` every raw p-value is ≤ 0.037 by
construction, so BH accepts 100% of candidates.

The threshold is now set empirically by two controls that constrain from opposite
sides:

**Permutation null (sets the floor).** Each transcript keeps its exact
differences — same number, same substituted bases — but they are scattered to
random positions. Contiguity is destroyed, everything else preserved, so any
tract found is false by construction. Over 20 replicates:

| m ≥ | IGL real | IGL null | IGH real | IGH null |
|---|---|---|---|---|
| 3 | 1776 | 0.80 | 1959 | 5.60 |
| 4 | 556 | **0.00** | 513 | **0.00** |

**AID spectrum (sets the ceiling).** The permutation scatters mutations, so it
cannot detect *clustered* SHM being admitted as tracts. If the threshold is too
loose, inside-tract differences start carrying an AID signature they should not
have:

| m ≥ | IGL tracts | IGL inside hotspot | IGH tracts | IGH inside hotspot |
|---|---|---|---|---|
| 4 | 556 | 0.96 | 513 | 1.29 |
| 5 | 123 | 0.89 | 99 | 1.12 |
| 6 | 93 | 0.92 | 42 | 0.91 |
| 7 | 36 | 0.90 | 15 | 0.72 |

Recomputed after the clonal-copy fix. An earlier version of this sweep counted
clonal copies and put IGL m ≥ 4 at ×1.26, which is what originally ruled it out.
It should not have.

**Chosen: m ≥ 5 for IGL, m ≥ 6 for IGH.** The permutation null already licenses
m ≥ 4 in both; these are the conservative choice, and the outside/inside contrast
is stable across m = 4–7 (outside 2.13–2.39, inside 0.72–1.29), so the conclusion
does not rest on where in that range the line falls. The loci differ because their donor
pools differ (23 vs 162 genes) — the same logic Bonferroni was reaching for,
calibrated empirically instead of analytically. Against the old rule this is a
real relaxation: **IGH 10 → 28 tracts, IGL 93 → 122**, with no loss of
specificity by any available control. The AID contrast is stable across
m = 5…8, so the conclusion does not depend on where in that range the line falls.

**Topology control.** Separately, a tract whose donor was already deleted is a
known false positive. If a fraction `e` of the pool was deleted, false calls hit
a deleted donor with probability `e`, so `FDR = f/e` with no distributional
assumption. On IGH deletional parents this gives FDR ≈ 1.0 at every threshold —
the donor *attribution* there is near-random, even though the tracts themselves
are real. In IGL nothing is ever deleted (`e = 0`) so this control has no power,
which is precisely why the permutation null exists.

#### Competing donors

Several donors can share the diagnostic bases and explain one tract. Those are
alternative explanations of **one** event, not several events. The detector ranks
them and flags the best as `primary_donor`; the arrow figures draw the primary in
colour and the runners-up in grey. Ambiguity is far worse in IGH (67% of tracts
at m ≥ 4 have a competitor) than in IGL (13%), which is the same fact the
topology control reports from the other direction.

### 11. The AID hotspot spectrum test — the strongest evidence

This is the test that actually separates the two processes, and it does not depend
on the donor pool at all.

#### What the motif names mean

They are not acronyms — each letter is an **IUPAC ambiguity code** standing for a
set of allowed bases:

| Code | Bases | Mnemonic |
|---|---|---|
| **W** | A or T | *weak* — 2 hydrogen bonds |
| **S** | G or C | *strong* — 3 hydrogen bonds |
| **R** | A or G | pu*r*ine |
| **Y** | C or T | p*y*rimidine |

AID deaminates **cytosine → uracil** in single-stranded DNA, and it prefers
certain neighbours. The targeted base is the C (or, read on the other strand,
the G):

| Motif | Expansion | Mutated base |
|---|---|---|
| **WRCY** | `[AT] [AG] C [CT]` | the **C**, 3rd position |
| **RGYW** | `[AG] G [CT] [AT]` | the **G**, 2nd position |
| **SYC** | `[GC] [CT] C` | the **C**, 3rd position |
| **GRS** | `G [AG] [GC]` | the **G**, 1st position |

Each pair is one motif seen from two strands: **RGYW is the reverse complement
of WRCY**, and **GRS is the reverse complement of SYC**. Because transcription
exposes both strands as single-stranded DNA, AID can hit either, so both spellings
have to be counted. Concretely, `AACT` is a WRCY hit and `AGTT` is the same site
read the other way.

So a hotspot and a coldspot are both statements about a C:G pair — one where AID
is likely to strike, one where it is not. That is exactly why moving the two in
**opposite** directions is specific to AID, while anything that merely clusters
mutations would move both the same way.

*Variants in the literature:* some papers write the hotspot as WRCH / DGYW
(H = A/C/T, D = A/G/T), which is slightly broader. This pipeline uses the
classic WRCY / RGYW definition throughout.

**Somatic hypermutation** is done by AID, which is not random: it targets WRCY /
RGYW motifs (hotspots) **and avoids** SYC / GRS (coldspots). **Gene conversion**
copies a block of donor sequence — the resulting differences sit wherever the
donor happened to differ, with no relationship to AID motifs.

Both arms matter. Hotspots alone cannot distinguish "AID acted here" from
"mutations are clustered here for any reason at all"; a process that merely
concentrates mutations would move hotspots and coldspots the **same** direction.
Only AID pushes them in **opposite** directions. So the test has four cases, and
the prediction is directional in each:

|  | hotspot (WRCY/RGYW) | coldspot (SYC/GRS) |
|---|---|---|
| **IGL outside tracts** | ×2.24 ↑ (p = 0.002) | ×0.40 ↓ (p = 0.002) |
| **IGL inside tracts** | ×0.55 ↓ (p = 0.002) | ×1.47 ↑ (p = 0.002) |
| **IGH outside tracts** | ×2.23 ↑ (p = 0.002) | ×0.51 ↓ (p = 0.002) |
| **IGH inside tracts** | ×0.78 (n.s.) | ×1.68 ↑ (p = 0.014) |

> **Counting unit: one observation per CLONE per tract**, clones defined by
> shared VDJ junction. See *The clonality test* below — the choice is not a
> convention, it is measured, and getting it wrong in either direction changes
> the result.


These are against the **tract-restricted null**, which is what panel A plots; the
gene-wide value is drawn on each bar as a dashed tick for comparison.

Outside-tract differences carry a textbook AID signature in both loci — enriched
at hotspots **and** depleted at coldspots, both at the permutation floor.
Inside-tract differences do not: in IGL they run the other way, depleted at
hotspots and enriched at coldspots. That is the prediction for a difference
copied from a donor — it sits wherever the donor happened to differ, unrelated to
where AID bound.

The null redistributes each transcript's mutations at random over the positions
actually covered in that same gene, keeping the count fixed (1000 permutations).
This controls for base composition, without which raw motif percentages are
meaningless.

#### The J strand in IGH is an assumption — and it was tested both ways

No J has been located for blackbird IGH. The **side** is not in doubt: the D
cluster (216.4–218.9 kb) sits above the V array (11.8–77.9 kb), and D must lie
between V and J, so J is above V whichever strand the locus reads. Only the
**strand** is inferred, from the D cluster's heptamer orientation (16/24 canonical
`CAC` on `+` against 8/24 on `−`).

That matters because strand is what assigns deletion versus inversion, and
swapping it exactly swaps the labels — 92 deletional / 70 inversional becomes
70 / 92. So the two hypotheses are distinguishable in principle. Running both
through the topology control:

| m ≥ | `j_strand "+"` | `j_strand "−"` |
|---|---|---|
| 4 | FDR 0.82 | FDR 0.91 |
| 5 | FDR 1.11 | FDR 0.88 |
| 6 | FDR 1.22 | FDR 0.91 |
| 8 | FDR 2.78 | FDR 0.92 |

**FDR ≈ 1 under both.** The impossible fraction simply tracks the expected
fraction, meaning IGH donor assignment carries no topological signal for the
control to work with. The test cannot choose a strand, and that failure is itself
the finding — it agrees with every other route by which IGH donor attribution has
turned out to be unreliable.

**Consequence:** no IGH result that depends on the deletion/inversion labelling
should be quoted. The count of tracts with no topologically possible donor is
8/28 under `+` and 13/28 under `−` — an artefact of the assumption, not a
measurement. The AID spectrum, which is the primary evidence, never touches the
donor pool and is unaffected. IGL is unaffected too: its J is directly located and
confirmed (567 junction-only hits at 6,339 kb).

#### The clonality test — what counts as one observation

Every enrichment in the AID test is a fraction over *differences*, so it depends
on deciding what an independent observation is. A conversion tract carried by 40
transcripts is **one** event if those transcripts are clonal descendants of a
single converted B cell, and **forty** events if they are independent
rearrangements that each acquired the same tract. Counting the first case forty
times is pseudoreplication and makes the permutation null far too narrow;
collapsing the second case to one throws away real replication. Both errors are
large, and they point in opposite directions.

This is decidable rather than assumable. The **VDJ junction** — the N-nucleotide
region immediately 3′ of the V gene — is generated once, at rearrangement, and is
inherited by every daughter cell. Clonal relatives share it; independent
rearrangements essentially never do.

**Method.** Take the 45 bp immediately 3′ of the V alignment end for each
transcript (5′ of the start, for minus-strand alignments), then greedily cluster
transcripts within a parent gene at ≥95% junction identity. An inside-tract
difference is counted **once per clone per tract**.

**Result — the two loci differ, and it matters:**

| | junction identity, transcripts sharing a tract | random pairs | transcripts → clones |
|---|---|---|---|
| **IGL** | median 0.267, 1% ≥95% | median 0.267, 0% ≥95% | 520 → 396 |
| **IGH** | median 0.722, **50% ≥95%** | median 0.267, 1% ≥95% | 68 → 54 |

In IGL, transcripts sharing a tract are **no more related than random pairs**.
They are independent rearrangements carrying the same conversion — *recurrent*
conversion, which is what a locus with one functional gene and a preferred set of
donors should produce. The inside class therefore barely shrinks (779 → 712) and
the result stands. In IGH they are substantially clonal, and collapsing changes
the numbers (inside hotspot ×0.61 → ×0.78, no longer approaching significance).

**Why this was worth doing.** An intermediate version of this analysis collapsed
to one observation per *tract* regardless of clonality. That over-corrected
badly: IGL's inside class fell to 92 differences and a real effect
(hotspot ×0.55, p = 0.002) disappeared into a null result (×0.90, p = 0.77).
Getting the counting unit wrong in either direction changes the conclusion.

**A second control, same question.** Recurrent tracts could also be explained if
the "tract" were simply an allele — the parent gene's other-haplotype version
appearing in transcripts from that chromosome. It is not: of 20 distinct IGL
tracts, **0** have their donor-diagnostic bases present in `bAgePho2_alt`, the
same bird's second haplotype. (This is the failure mode that sank the duck
comparison, so it is worth checking explicitly.)

Flags: `--transcripts` enables clone-aware counting; omitting it falls back to
one observation per tract; `--count-clonal-copies` disables collapsing entirely.

#### The two nulls, in one paragraph

> Both nulls ask the same question — do the differences land on AID motifs more or
> less often than chance? They differ in what "chance" means: which positions a
> difference **could** have occupied.
>
> - **Gene-wide:** anywhere in the V gene.
> - **Strict:** only within the same region it was actually found in — a
>   tract-internal difference stays inside its own tract window, an outside
>   difference stays outside.
>
> Only the strict version matches the background composition of the region being
> tested.

In one line: *gene-wide uses the whole gene as the background; strict uses the
background of the specific region each difference came from.*

Note that **both** nulls are computed for **both** classes — strict is not "the
inside one". It reframes the outside class too; it simply changes almost nothing
there, because outside positions are ~95% of the gene (255 of 258 covered
positions in IGL).

This matters when tract windows are compositionally unusual. In IGH they are —
**1.5× hotspot-rich and 2.4× coldspot-poor** relative to the gene average —
because AID targets the sites where conversion initiates. Using the whole gene as
background there compares inside-tract differences against composition they never
experienced, which flipped the coldspot result from ×0.69 to ×1.66. In IGL the
tract windows are compositionally ordinary, so both nulls agree. The two nulls
diverge **only** where the composition differs, which is what makes this an
explanation rather than a rationalisation.

#### Exactly how the p-values are computed

**1. Collect.** Each transcript is projected onto its parent's coordinates. For
every covered position `i` (excluding the last 20 bp), record whether the
transcript differs from the parent, and classify the difference `inside` or
`outside` according to whether `i` falls in one of that transcript's tract
intervals. The motif is always evaluated on the **parent germline sequence** at
`i`, not on the transcript — the question is whether AID would have targeted that
site, which depends on the DNA it acted on.

**2. Observed statistic.** `hotspot_fraction = (differences at hotspot positions)
/ (all differences in that class)`, and likewise for coldspots. A raw fraction is
meaningless on its own, because motif density varies between genes and along a
gene — hence step 3.

**3. Null distribution (1000 replicates).** The mutation *count* is held fixed and
the *positions* are resampled without replacement. Two framings are computed:

| Null | Sampling frame | Controls for |
|---|---|---|
| **gene-wide** (default) | for each gene, draw that gene's mutation count from all covered positions in it | gene-to-gene composition and how many mutations each gene contributed |
| **tract-restricted** (strict) | for each transcript, draw its inside-tract count from that transcript's *within-tract* positions, and its outside count from its outside positions | additionally, **where the tracts sit** |

The strict null exists because the gene-wide one quietly asks the wrong question
of the inside class: those differences are confined to a few short windows by
construction, so if tracts happen to land where hotspots are scarce, the class
would look AID-depleted because of *where the tracts are* rather than what
happened inside them. This matters here — most IGL events cluster at V:137–149.

**4. p-value.** Empirical, with the standard +1 correction (Davison & Hinkley):

```
p_enriched = (#{null >= observed} + 1) / (n_perm + 1)
p_depleted = (#{null <= observed} + 1) / (n_perm + 1)
p_two_sided = min(1, 2 * min(p_enriched, p_depleted))
```

The +1 stops a p-value of exactly zero, which no finite permutation set can
justify. **With 1000 permutations the smallest achievable two-sided p is
0.002** — that is a resolution floor, not a coincidence, and it is why every
strong result in these tables reads `p = 0.002`. A smaller number would require
more permutations, not better data.

**5. Enrichment** = observed / mean(null).

#### Does the tract-restriction change the answer? No.

| | gene-wide | tract-restricted |
|---|---|---|
| IGL outside hotspot | ×2.23 | ×2.24 |
| IGL outside coldspot | ×0.40 | ×0.40 |
| IGL inside hotspot | ×0.48 | ×0.54 |
| IGL inside coldspot | ×1.56 | ×1.47 |
| IGH outside hotspot | ×2.18 | ×2.23 |
| IGH outside coldspot | ×0.51 | ×0.51 |
| IGH inside hotspot | ×0.90 (n.s.) | ×0.61 (p = 0.09) |
| IGH inside coldspot | ×0.69 (n.s.) | **×1.66 (p = 0.004)** |

The outside class is untouched, as expected — it covers most of the gene either
way. The inside class moves, and it moves *toward* the predicted pattern: under
the correct frame IGH inside-tract differences are significantly **enriched** at
coldspots, the same anti-AID mirror IGL shows. The gene-wide null had been
diluting that signal, not manufacturing it.

Both nulls are written to `{locus}_aid_spectrum.tsv` (`*_strict` suffixes).

**Circularity check:** if the tract detector were simply picking non-hotspot
positions by construction, this would be an artifact. It isn't — donor-difference
positions are hotspot-neutral (10.7% vs 10.2% background).

> Note that the tract threshold was chosen partly using this test (step 10), so
> at the chosen m the AID result is not fully independent evidence. What rescues
> it is stability: the contrast holds across m = 5…8, so no particular choice
> manufactures it.

---

### 12. Are closer donors used more often?

A long-standing claim from chicken IGL is that gene conversion prefers donors
close to the rearranged gene: close in sequence (the most similar pseudogenes
supply most tracts) and close on the chromosome. `gc_donor_distance.py` tests
both, per locus, then pools.

**Why the obvious version of this analysis is worthless.** Tracts are called from
*informative positions* — places where parent and donor differ. A donor 99%
identical to its parent leaves almost nothing to detect and can never produce a
significant tract however often it is used; a donor 85% identical offers
informative positions everywhere. Detection power rises with distance, which is
the axis under test, in the direction opposite to the claim. In blackbird IGL the
correlation between plain window-counting opportunity and divergence is
ρ = +0.89. Plotting usage against identity measures the detector.

So the test conditions on detectability. For each (parent, donor) pair it counts
the places a significant tract *could* have been called — windows of m
informative positions spaced no more than `max_gap_bp` apart, inside the region
the transcripts cover, using the same m and gap the detector ran with. Pairs with
zero such windows never entered the competition and are dropped (4 of 21 in
blackbird IGL, at 84.0–92.3% identity; 52 of 155 in IGH, up to 98.3%).

**A second confound, which the window count does not fix.** A donor nearly
identical to its parent differs from it only in a few clusters, and those
clusters are the hypervariable ones — which is exactly where the transcripts also
differ from the parent. The few windows a similar donor offers therefore sit
where a hit is most likely, so plain window counting understates its
detectability and manufactures the result. The default `--offset expected`
weights each window by the parent's own per-position mutation frequency,
computed **outside** the transcripts' own called tracts so the correction does
not condition on the events under test. This matters enormously: in blackbird IGL
the divergence effect is RR 0.78 per 1% (p = 2×10⁻⁸) under plain window counting
and RR 0.90 (p = 0.11) under the weighted version. Most of the apparent
preference for close donors is detection bias.

**Three tests, same conditioning.** A conditional logit over each parent's donor
pool (the primary test — parent is the stratum, so pool composition, array
position and transcript depth all cancel); the same null as a permutation, which
assumes nothing asymptotic at these counts; and plain Spearman, reported only to
show what the uncorrected analysis would have claimed.

**Physical distance is parameterised three ways, because they are not
interchangeable.** `distance_kb` is linear kb, the natural scale for arrays this
size; `rank_distance` counts intervening V genes, which is how "the nearest
pseudogenes supply most tracts" is normally meant and is immune to uneven gene
spacing; `log10_distance` is per tenfold, which is a poor model here — blackbird
IGL spans 1.9–22 kb, barely one decade. With `--j-pos` a `toward_j` term also
asks which *side* of the parent the donor sits on, since a locus can show no
distance effect and still be polarised.

**Sequence and physical distance are confounded with each other**, so both go
into one model as well as separately. In a tandem array neighbours are recent
duplicates and therefore also similar, and here the correlation is *negative* in
both IGL loci (ρ = −0.31, p = 0.002 in duck IGL): more distant donors are more
similar. A marginal slope on either covariate carries the other inside it. The
`*_given_*` rows in the report are the partial slopes, each holding the other
fixed, and those are the ones to quote.

**The unit is the distinct event**, one (parent, donor, start, end). Transcripts
are clonally related, so counting them counts clone size — the report gives all
three units and the transcript-level p-values are, predictably, absurd
(p = 9×10⁻⁴⁴ in IGL from 19 independent events).

**Ambiguous tracts are dropped by default.** Where several donors explain a tract
equally, the detector's tie-break is alphabetical, i.e. by contig position, which
would fabricate a genomic-distance signal. `--ties keep` includes them at weight
1/k as a sensitivity check.

#### The answer

Partial slopes, event unit, each holding the other covariate fixed:

| locus | events | per 1% divergence | per 10 kb along the locus | per intervening V gene |
|---|---|---|---|---|
| blackbird IGL | 19 | RR 0.92 | RR 1.62 | RR 1.05 |
| blackbird IGH | 13 | RR 0.78 | RR 1.52 | RR 1.03 |
| duck IGL | 2 | RR 0.86 (not pooled) | RR 5.51 (not pooled) | RR 1.17 (not pooled) |
| duck IGH | 15 | RR 0.97 | RR 1.36 | RR 1.02 |
| **pooled** | | **0.967 [0.932, 1.003], p = 0.069** | **1.487 [1.107, 1.999], p = 0.0085** | **1.028 [1.006, 1.050], p = 0.013** |

**Sequence proximity: direction consistent, effect small, not significant.** All
four slopes are negative — closer-in-sequence donors used more — a 4/4 sign test
at p = 0.063, pooled just off the line at p = 0.069, no heterogeneity. Over the
observed divergence range this is roughly a threefold preference for the most
similar donors over the most distant. Consistent with the published claim, and
not on its own evidence for it.

**Physical proximity: the effect is significant and runs the *opposite* way to
the published claim.** Donors further along the locus are used *more*, by about
1.5× per 10 kb and 1.03× per intervening gene, all four loci in the same
direction, with no heterogeneity (Q = 0.14, p = 0.93). This is the most
statistically secure result in the analysis, and it is not what chicken IGL
reports.

Three checks say it is not an artefact. It survives mutual adjustment and in fact
*strengthens* (marginal RR 1.41 → partial 1.49), so it is not the sequence effect
in disguise — the two covariates are largely independent here. It holds under
both detectability models, unweighted window counting and the SHM-weighted one,
which is the correction that demolished most of the sequence effect. And it is
not censoring: the near-identical pairs dropped for zero opportunity are *not*
systematically the nearby ones in the three pooled loci (Mann-Whitney p = 0.90,
0.25, and none dropped), though they are in duck IGL (p = 0.0098), which is
excluded from the pool for having two events.

**Directionality is untestable in two of the four loci.** `toward_j` is not
identified in blackbird IGL or duck IGH: each has a single tract-bearing parent
sitting at the end of its array, so every candidate donor lies on the same side
of it and there is no contrast. Where it can be fitted (blackbird IGH) it is flat
(RR 1.23, p = 0.79). This is the same vacuity the topology test hits, reached
from a third direction — the tract-bearing parents sit exactly where positional
tests have no power.

Three caveats bound all of this. Every effect is measured against a
detectability model, and a model that is wrong in either direction moves the
answer — the swing between the two offsets above shows how much. Duck IGL
contributes two events; it is drawn in the forest plot and counted in the sign
test, but its slope is not identified and is excluded from the pool
(`--min-events`). And the analysis is silent about the very closest donors: pairs
with zero detection opportunity — up to 98.3% identity in IGH — never entered the
competition, and those are precisely the class the chicken claim is about. What
is tested is donor choice among donors that *could* have been seen.

---

## What did *not* work — stated plainly, because two figures show it

**The topology test is vacuous in both loci.** It was designed to test the
hypothesis that a functional gene at the J-proximal end is favoured because all
donors stay available. The result: in IGL the parent has **zero** deleted donors,
and in IGH every tract-bearing parent is inversional. There is no contrast to
test against.

That is not a failed analysis — it *is* the hypothesis confirmed, arrived at from
the other direction. The functional genes sit exactly where the test has no power,
which is precisely what the hypothesis predicts. It just cannot be reported as a
significant p-value, and `IGx_report_*.pdf` panel C says so on its face.

**D gene usage is not determinable from this data.** The IGH junction has a median
length of 32 bp; the longest shared stretch with the best-matching D gene averages
7 bp against 6.5 bp expected by chance. Only 4 of 70 transcripts give a confident
call, and those average 8 tied D genes. The figure is included because the
question was asked and the answer is a real result — the data cannot resolve it —
but the D usage numbers should not be quoted as usage.

**BrepConvert does not pass the topology control.** In IGH, 24% of its events use a
donor that was already deleted. An early reading that it *did* pass was a power
artifact: 95% of its events sat on a parent with zero deleted donors, where the
test cannot fail.

---

## Figure guide

### `FIG_main_gene_conversion.pdf` — the summary figure

One page carrying the whole argument. Built by `scripts/gc_main_figure.py`.

| Panel | What it establishes |
|---|---|
| **A** | IGH architecture — 162 V genes, only 25 with an RSS, 35 expressed |
| **B** | IGL architecture — 23 V genes, only **2** with an RSS, 1 functional parent |
| **C** | Raw sequence evidence. The heading counts the **whole dataset** (39 distinct tracts from 150 transcript-level calls; IGH 19/28, IGL 20/122) but only **six rows are drawn** — the three best-supported per locus, ranked by donor-diagnostic positions. It is an illustrative subset, deliberately the strongest cases, tagged IGH/IGL and listed IGH-then-IGL to match panels A and B. Gene ids are omitted: the row identity is the point, not which pseudogene it was. Where parent and donor agree the transcript base carries no information, so it takes the same neutral tone as those rows; only the three informative outcomes (follows donor / follows parent / follows neither) get their own colour |
| **D** | IGL donor→parent network: 10 of 22 possible pairs, **0/122 impossible**. `--network-extra` adds further loci as stacked rows in the same slot (the duck figure uses it for IGH), which shifts the AID panel's letter accordingly |
| **E** | Mutation spectrum, all four cases, both loci, against the tract-restricted null. Hotspots are drawn solid and coldspots as a lighter weight of the same locus colour — they are the same measurement on the same locus, so a second texture would imply a second variable |

The argument runs A→B (almost nothing can rearrange, so diversity cannot come
from combinatorial V use), then **C** (here are the actual events, as sequence),
then D (every part of the array feeds the one gene that can, and never from a
donor recombination had removed), then E (the differences inside those blocks do
not carry AID's targeting footprint, while those outside do — which is what makes
them conversion rather than hypermutation).

**C deliberately precedes D and E.** Both of those panels talk about tracts, and
E in particular contrasts differences *inside* versus *outside* one. Those words
mean nothing until the reader has seen a tract, so the sequence panel comes
first.

**Panels A and B.** Stem *height* is the transcript count on a **log** axis,
normalised by a `vmax` **shared between the two panels**. Sharing it is what
makes A and B comparable at all: the same stem height means the same number of
transcripts in both, so IGH's ceiling of 7 is visibly nowhere near IGL's 233
rather than being rescaled to fill its own panel. The two rows are also sized in
proportion to their axis ranges, so a given count occupies the same number of
millimetres in each — without that, a shared `vmax` would still draw the same
count at different sizes.

Both panels carry the **same tick values** (1, 5, 10, 50, 200, as far as each
reaches), labelled on **both strands**; the lower half previously had gridlines
but no numbers. Since log compresses exactly the difference the panel is being
asked to show, the top genes additionally carry their **count as a printed
label**. Stem *direction* is strand. A navy dot at the stem tip marks an RSS.

IGL gets a shorter row because it has a single minus-strand gene, so a symmetric
axis would spend half the panel on empty space.

The position axis is in **kb below a megabase and Mb above** — IGH reads 10–80 kb,
IGL 6.315–6.337 Mb. Raw base pairs gave either six-digit ticks or a detached
`1e6` offset in the corner, both harder to read than the number itself. The two
panels are different contigs, so there is nothing to compare between their
absolute coordinates and no reason to force one unit on both. The same rule now
applies to the locus maps in `IG_overview.pdf`.

A gene with **no transcripts** is drawn as a small open circle sitting *on* the
zero line, nudged just far enough above or below it to show which strand it is
on. The offset is kept well under `height(1)` (≈0.13 on the shared scale) — at
the earlier 0.115 a silent gene was drawn level with the "1" tick and read as
having one transcript, which is exactly the thing it does not have. A marker
rather than a stem, because marker size is in points and so survives whatever the
log axis does to the low end. In panel A the J gene is 140 kb beyond the V array, so
the axis breaks (`//`) rather than compressing every gene into the left third.

**Panel C.** Every gene gets the same large ringed marker — navy filled if it has
an RSS, white if not — and there are no position labels. Size previously encoded
expression, which made the donor-only genes (the entire subject of the panel) the
hardest things on it to see; the only distinction kept is the one the panel is
about, namely which gene can be rearranged at all.

> **Panel D is not "AID vs not AID".** AID initiates *both* processes — it makes
> the same lesion either way. What differs is the repair:
>
> - **SHM** — the lesion is resolved by error-prone repair *at that base*, so the
>   resulting mutations pile onto AID's target motifs.
> - **Gene conversion** — the lesion is resolved by copying a donor, so the
>   resulting differences sit wherever the donor happened to differ from the
>   parent, which has nothing to do with where AID bound.
>
> So the panel contrasts **SHM against gene conversion**, and what it detects is
> the *footprint of AID targeting* in the outcome, not the presence or absence of
> AID. Labelling it "not AID" would misstate the biology.

**Panel D.** Coldspots are drawn as a lighter weight of the same locus colour
rather than a hatch: they are the same measurement on the same locus, so a second
texture would imply a second variable. The significance key lists only the levels
that actually occur in the figure — here `***` and `n.s.` — since describing a
`*` that appears nowhere sends the reader hunting for it.

There is deliberately **no overall title** — that belongs in the manuscript text —
and every panel title is black, including the two locus panels, so colour is not
doing double duty as both a locus code and a heading style.

### `SUPP_{IGL,IGH}_donor_distance.pdf` — donor choice vs distance, per locus

One row per covariate — sequence divergence, linear kb along the locus, array
rank, log10 bp, and the J-side indicator. **Left**, the confound — detectability against the covariate,
rose for donors that were used. **Middle**, what the uncorrected analysis sees:
events per donor against the covariate. **Right**, the test — the observed mean
covariate of used donors (rose line) against 20,000 draws that pick donors from
the same pool in proportion to detectability alone. Grey points are candidate
donors never used. See §12 for why the left panel is the whole problem.

### `SUPP_donor_distance_pooled.pdf` — the four loci together

Forest plot of the conditional-logit slope as a usage rate ratio, one panel per
covariate. The bottom row holds the partial slopes — each covariate at fixed
value of the other — and those are the ones to quote. Below 1
means closer donors are used more, i.e. the published claim. Marker area scales
with events; whiskers are 95% CI, arrow-capped where they run off the axis. The
rose diamond is the fixed-effect pool. Duck IGL (2 events) is drawn in grey and
excluded from the pool — its slope is not identified — but still counts in the
sign test, which uses direction only.

### `IG_overview.pdf` — start here

| Panel | What it shows |
|---|---|
| **A** | Funnel from raw reads to IG transcripts — where the 91.3M reads went |
| **B** | IGH vs IGL split |
| **C** | Divergence from germline |
| **D**, **E** | Top-20 IGH / IGL V genes by usage |
| **F** | minimap2 vs immunotools agreement |
| **G**, **H** | **Locus maps.** One tick per V gene along the contig. Up/down = strand; **tall = has an RSS**, short = none; colour = how many transcripts use it as parent (viridis_r, dark purple = high) |

G and H are the ones to look at first — they show the whole architecture at once:
where the rearrangeable genes are, which way they point, and which ones the
repertoire actually draws on.

### `SUPP_IGH_donor_network.pdf` — blackbird IGH network, stripped down

The same detector output as `IGH_arrows_detector.pdf`, drawn for a supplement:
gene positions removed from the x axis, a one-line title, and **every arrow the
same width** so the figure says *which* genes exchange sequence rather than how
often. Gene colour still carries the only distinction the panel is about — navy
if the gene has an RSS and can be rearranged, grey if it can only donate. `--only-involved-genes` restricts the track to genes that appear as a donor or a
parent — 31 of 162 in blackbird IGH — which is what makes the remaining markers
legible at full size. The count is stated in the title. Markers are still sized
from the spacing available, so the option and the sizing work together rather
than fighting.

Built with `--hide-gene-labels --uniform-arrows --uniform-gene-size
--short-title --compact`.

> **Why `--compact` flattens the arcs rather than shortening the canvas.**
> `arc3` is evaluated in **display** coordinates: the apex sits `rad × chord / 2`
> *pixels* above the axis, so whether an arc fits depends on the panel's aspect
> ratio and not at all on `ylim`. A full-width chord at `rad 0.42` needs an axes
> height of `0.42 × width`; at 15 × 5.4 in that is 6.3 in against ~4 in
> available, which is why arcs ran into the title and legend. Compact mode
> therefore drops `rad` to 0.36 at a shorter 15 × 5.4 in canvas. The main figure's network panels go further and
> **measure** the drawn panel, clamping `rad` so the longest chord still lands
> inside — a fixed value cannot work there, because the panels differ in width
> and in whether arcs run one way or both. This is also why *widening* a network
> panel makes overflow worse rather than better.

**What each colour means, and which side it sits on.** Exactly **one arc is
drawn per event**, and its colour says how the sequence evidence and the
topology agree:

| colour | meaning | arc runs from | side |
|---|---|---|---|
| **teal** | the best-matching donor is one the rearrangement left intact | that donor | above |
| **grey** | the best-matching donor was deleted, but another donor could have supplied the tract | **that other, surviving donor** | above |
| **rose** | every candidate donor had been deleted | the best-matching (deleted) donor | below |

Teal and grey are both relationships that *could have happened*, which is why
they share the upper half; the difference is only whether the gene the sequence
matched best is the gene the topology permits. Rose is the one case with no
possible explanation at all, and it gets the lower half to itself. A grey arc is
therefore a weaker claim than a teal one — the event is real, but the donor named
is the fallback rather than the best match — and never an impossible one.

Because the arc for a grey event is redrawn *from the surviving donor*, an
illegal donor is never rendered when a legal alternative exists: drawing it would
assert a relationship that could not have occurred. If a grey pair already
appears as a teal call from some other event the grey copy is skipped (the two
arcs would be identical), and the legend then drops the grey entry rather than
promising an arrow the reader cannot find.

**Keeping every arc visible.** `arc3` curvature is applied in *display* space, so
an arc's apex sits `rad × chord / 2` **pixels** off the axis — independent of
`ylim`. Both scripts therefore measure the drawn panel and clamp `rad` so the
longest chord still lands inside, and the main figure clamps each side against
the height that side will actually be given (the axis opens asymmetrically when
only one side carries arcs, so a single long rose arc must not flatten every teal
one).

> **`--transcripts` is not optional.** Without it every member of an expanded
> clone contributes its differences independently, which is pseudoreplication:
> 520 IGL transcripts are only 396 clones, and the inside-tract class is the one
> dominated by expanded clones. Omitting it changes the IGL inside-tract hotspot
> from ×0.49 to ×0.98 and silently disagrees with the published figure. The
> Snakefile rule now passes it; it previously did not, so any `{locus}_aid_spectrum.tsv`
> produced by an older run of the pipeline should be regenerated.

The main-figure networks draw the **J gene** as a teal diamond at whichever end
it sits. Without it the duck IGH panel showed an unexplained fan of red: J lies
below that V array while the only expressed parents sit at the far end, so a
deletional rearrangement removes every donor between. J is what makes the colour
legible.

### `IGx_arrows_*.pdf` — donor → recipient

The figure asked for: genes laid out along the locus, one arrow per donor→recipient
relationship, labelled with the number of supporting transcripts.

- **Arrow direction**: donor → recipient.
- **Arrow thickness and its number**: transcripts supporting that pair. The number
  now sits directly on its own arc (see note below).
- **Above the axis / teal**: the best-matching donor survived the rearrangement.
- **Above the axis / grey**: the best-matching donor was deleted, but another
  donor could have supplied the tract — the arc is drawn **from that surviving
  donor**, so it is still a possible relationship, just not the best-matching one.
- **Below the axis / rose**: **impossible** — *every* candidate donor for that
  tract had been deleted. This is the visible false-positive load.

> **The colour rule is per tract, not per donor.** If any legal donor can explain
> a tract, that tract is never counted as impossible — even when a
> better-supported donor for it was deleted. The event is real; we just cannot
> say for certain which gene supplied it. Rose is reserved for tracts with *no*
> legal explanation at all.
>
> This matters a great deal for BrepConvert, which lists ~15 donors per event.
> Marking it wrong for every impossible donor it happens to mention — even when
> it also named a perfectly good one — put its IGH false-positive rate at
> **38%**. Under the per-tract rule it is **14%** (25 of 176 tracts). The same
> definition is what `gc_compare_methods.py` already called `impossible_lenient`.
- **Marker colour**: navy = has RSS, grey = none.
- **Marker shape**: points *toward* J = deletional, *away* = inversional.
- **Marker size**: large = at least one transcript best-matches it, scored with
  **every** V gene eligible (`unconstrained_assignments.tsv`). Colour says
  whether a gene *can* rearrange; size says whether the data *point at it*.
  Keeping those independent is the whole value of the two channels — a large
  grey marker (expressed, no RSS) is an interesting object, either a missed RSS
  or a donor so heavily copied that transcripts drift onto it.

  > Earlier versions took size from `functional_genes.tsv`, which only counts
  > transcripts for RSS-bearing candidate parents. A donor-only gene therefore
  > had `n_transcripts = 0` by construction and could never be drawn large, so
  > size silently restated colour and disagreed with the overview locus map.
  > Fixed; pass `--usage-assignments` to `gc_plots.py`.
- **Teal diamond**: the J gene, at whichever end it sits.

Suffixes: `_detector` = our detector, `_brepconvert` = BrepConvert. In the VGP run
the variants are lettered `A`–`E` (raw BrepConvert, donor-resolved, our detector
best-match, our detector RSS-parent, BrepConvert RSS-parent) so the effect of each
choice is separable.

> **Arrow labels.** These were previously misplaced. `FancyArrowPatch` builds its
> arc in *display* coordinates, not data coordinates, so no data-space formula
> predicts where an arc peaks; the labels drifted off their curves and most were
> pushed off-canvas. They are now computed from matplotlib's own `Arc3` Bezier,
> pulled back into data coordinates, and where two would collide the later one
> slides **along its own arc** instead of being lifted off it.

### `IGx_report_*.pdf` — the six diagnostic panels

| Panel | What it shows |
|---|---|
| **A** | Significant tracts per donor–recipient pair |
| **B** | Do functional genes sit where donors are plentiful? |
| **C** | Topology control, power-matched — *this is the panel that reports the test is vacuous* |
| **D** | Where tracts fall along the V gene |
| **E** | Are impossible donors as well supported as possible ones? (if yes, the method is not discriminating) |
| **F** | Expression vs donor supply |

### `IGx_matrix_*.pdf`

Donor × recipient heatmap, viridis_r (dark purple = many tracts).

### `IGx_aid_spectrum.pdf` — the key evidence figure

**Only panel A is evidence.** B and C are reported for completeness and do not
support the hypothesis — see below.

| Panel | What it shows |
|---|---|
| **A** | All four cases: {hotspot, coldspot} × {outside, inside}, as fold-change against the **tract-restricted** permutation null, on a log axis so enrichment and depletion are symmetric about 1×. A dashed tick on each bar marks where the gene-wide null would have put it |
| **B** | Transition bias — **confounded, not evidence** |
| **C** | C:G targeting — **no signal against the correct baseline** |
| **D** | Full substitution spectrum (descriptive) |

Panel A is the result, and the shape to look for is the **mirror image**:
outside-tract bars go up on hotspots and down on coldspots; inside-tract bars do
the reverse (IGL) or nothing (IGH).

**Why panel B is not evidence.** The prediction is that outside (SHM) should be
more transition-biased than inside (conversion). The loci disagree:

| | outside | inside | verdict |
|---|---|---|---|
| IGL | 0.77 | 1.28 | **backwards** |
| IGH | 0.89 | 0.80 | weakly consistent |

The confound: conversion copies differences accumulated between **paralogues over
evolutionary time**, and molecular evolution is itself strongly transition-biased.
Inside-tract differences inherit an ancient transition bias that has nothing to do
with AID, so Ti/Tv cannot separate the two processes here.

**Why panel C is not evidence.** An earlier version drew a hardcoded 0.5 baseline,
which was simply wrong — the correct expectation is the C+G content of the covered
positions (0.586 in IGL, 0.641 in IGH), taken from the same permutation as panel
A. Against the proper baseline there is no C:G enrichment anywhere:

| | observed | expected | ratio |
|---|---|---|---|
| IGL outside | 0.495 | 0.586 | 0.85 |
| IGL inside | 0.558 | 0.586 | 0.95 |
| IGH outside | 0.621 | 0.641 | 0.97 |
| IGH inside | 0.546 | 0.644 | 0.85 |

This does not contradict panel A: hotspot motifs are C/G-centred by definition,
but they account for only ~20% of outside-tract differences, so a 2.2×
enrichment there barely shifts the *global* C:G fraction.

**For the paper:** panel A in the main figure; B and C to supplementary with the
confound stated, or dropped.

### `IGH_d_usage.pdf`

| Panel | What it shows |
|---|---|
| **A** | Junction length distribution (median 32 bp) |
| **B** | Is the D match longer than chance? (shuffled-junction null) — **it is not** |
| **C** | Call confidence — 4 of 70 assignable |
| **D** | D gene usage, weighted by how many D genes tie |

All 51 annotated D genes are used regardless of strand, as requested. Read this
figure as "the junction is too short to identify D genes", not as a usage profile.

### `SUPP_germline_matched_vs_not_{IGH,IGL}.pdf`

The same transcripts scored against a different bird's germline (`bAgePho0`) and
against the matched one, so every difference is caused by the reference alone.

| Panel | What it shows |
|---|---|
| **A** | Apparent divergence from germline — a mismatched reference inflates it, and that inflation gets read as SHM or conversion that never happened |
| **B** | Per-transcript identity gain from switching to the matched reference |
| **C** | **Did the transcript get the same parent?** |
| **D** | **How many V genes look real?** |

**Panel C categories** — genes are paired between the two references by best
reciprocal alignment first, so "different" means genuinely different, not just
differently named:

- *same gene (agree)* — the transcript's best parent in reference A is the
  orthologue of its best parent in reference B. The reference didn't matter.
- *DIFFERENT gene (wrong parent)* — the gene **has** an orthologue in the other
  reference, but the transcript was assigned somewhere else entirely. This is the
  damaging case: parent assignment is what every tract, donor and mechanism call
  is defined relative to, so a wrong parent corrupts all of it.
- *no equivalent gene (unanswerable)* — the gene has **no** best-reciprocal
  partner in the other reference at all. That allele simply isn't in the other
  assembly, so the question "same or different?" has no answer. It is kept as its
  own category rather than folded into either, because calling it agreement would
  flatter the mismatched reference and calling it disagreement would overstate
  the damage.

**Panel D categories** (this was previously labelled "genes called used", which
was opaque):

- *look expressed (≥0.98 identity to a transcript)* — V genes you would call
  expressed on the usual near-perfect-identity criterion. With the wrong
  reference this collapses, because that individual's alleles aren't present to
  match against. In IGH it goes **0 → 2**: against the different bird, *not one*
  IGH V gene clears the bar.
- *attract ≥1 transcript* — V genes that are some transcript's best hit, however
  poor. Much less sensitive to the reference, which is the point: the permissive
  count looks stable while the meaningful one collapses.

### `SUPP_reference_choice_across_haplotypes.pdf`

The supplemental figure making the case that matched germline data is necessary.
Five haplotypes, the same transcripts throughout, all scored against `bAgePho2_pri`.

All rows use **exactly the same transcripts** (IGL n = 520, IGH n = 69). That is
enforced with `--restrict-to`, not assumed: the 200 bp coverage floor is applied
per reference, so a borderline transcript can clear it against one germline and
fall under it against another. Before this was fixed the IGL rows ran 521/521/520,
which is harmless arithmetically but contradicts the figure's own claim that only
the reference changed. The comparison is now run on the intersection across every
reference.

The design point is the **control**: `bAgePho2_alt` is the *same bird's other
haplotype*. Without some such baseline a discordance rate means nothing, because
you cannot separate "wrong bird" from "no two assemblies ever agree perfectly".

**What that control is not.** It is *not* a technical noise floor. The bird is
diploid and both haplotypes are transcribed, so a V gene on the alt haplotype is a
genuinely expressible gene with its own alleles. When a transcript is assigned to
an alt gene rather than the pri orthologue, that can be the **correct** answer —
the transcript may really have come from the alt allele. So the 21.5% mixes:

- real allelic origin (the transcript came from the other haplotype's copy);
- genuine content and copy-number differences between the two haplotypes;
- assembly differences (fragmentation, collapsed or duplicated genes).

What it actually measures is the floor imposed by scoring a **diploid animal
against a haploid reference** — the irreducible cost of picking one haplotype.
That cost applies to the matched reference too: `bAgePho2_pri` is also only half
this bird's germline.

**Why it is still the right control.** Every different-bird reference pays that
same haploid-reference cost *and* adds between-individual divergence on top, so
the excess over the same-bird row still isolates the effect of using a different
animal. If anything the baseline is generous — part of the 21.5% is real biology
rather than error, so the true technical floor is lower and the different-bird
penalty larger than the subtraction implies.

| Reference | Locus | Different parent | Median identity gain |
|---|---|---|---|
| bAgePho1_pri | IGH | 75.4% | +1.35% |
| bAgePho1_alt | IGH | 52.2% | +0.88% |
| bAgePho0_pri | IGH | 43.5% | +1.24% |
| bAgePho1_alt | IGL | 31.1% | +0.00% |
| bAgePho0_pri | IGL | 30.4% | +0.00% |
| bAgePho1_pri | IGL | 30.0% | +0.15% |
| **bAgePho2_alt** | IGL | **21.5%** | +0.00% | ← same bird, other haplotype |


The result to quote: in IGL the haploid-reference floor is 21.5%, and the
different-bird references sit at 29.9–31.1% — roughly **9 points of excess
attributable specifically to using a different animal**, on top of a floor that is
itself partly real biology. And the third panel is the warning:
median identity gain is **0.00%** for IGL. No summary statistic would have told
you anything was wrong. The sequences look fine; the *assignments* are wrong.

The identity-gain scatter that used to sit alongside these two panels has been
removed — it restated what the bar panels and the table already say. The numbers
are still in `SUPP_reference_choice_across_haplotypes.tsv`
(`median_identity_gain_pct`), which is where that comparison belongs.

> These five haplotypes went through the **reference-swap comparison only** — the
> transcripts were re-aligned to each haplotype's V genes and re-scored. They were
> **not** put through tract detection, BrepConvert, the AID spectrum test, D usage
> or immunotools. That is why it took minutes rather than the hours the full
> matched run took.

### Stage-1 figures

| File | Contents |
|---|---|
| `summary_plot.pdf` | Per-sample exact matches and top alignments |
| `combined_summary_plot.pdf` | Same, across samples |
| `identity_analysis_plot.pdf` | Identity distributions |
| `pooled_summary_plot.pdf` | Pooled repertoire summary |
| `alignment_position_map.pdf` | Multi-page: per-gene pileup of where transcripts align along each V gene |
| `rank_comparison_plot.pdf` | Multi-page: rank-1 vs rank-2 gene assignments — how separable the best hit is from the runner-up |

The two multi-page reports are where to look when a gene assignment seems
doubtful: `rank_comparison_plot` shows directly whether the top gene actually beat
the second one or whether they are effectively tied.

---

## Second species: tufted duck

The same stage-2 analysis and the same figures were run on a tufted duck IsoSeq
dataset (`/local/storage/kav67/tufted_duck/results/geneconv/`), config
`config_geneconv_tufted_duck.yaml`.

**Stage 1 was not re-run** — it was already complete for all six SRR accessions
and nothing in it changed. `skera` does not apply at all: this is plain
Sequel-era IsoSeq, not a Kinnex array, so there is no concatenated array to
split. Because that dataset clustered directly without the IG prescreen,
`--prefilter-stats` is now optional in `gc_overview.py` and the funnel simply
starts at the clustered transcripts.

Three things had to be decided for this dataset rather than inherited:

* **Pooling.** The six runs are 4–484 IG transcripts each and split by locus
  (SRR13570364 is 92% IGH, SRR13570377 is 90% IGL). Transcript IDs already carry
  their accession, so the PAFs concatenate without collision: 483 IGH + 283 IGL.
* **The strict transcript set, not the merged one.** `merge_stats` reports
  `in_both = 0` for every run — immunotools ran with different sequence IDs, so
  its 46k calls cannot be reconciled with minimap2's 484 and their "union" is
  meaningless.
* **Recalibrated thresholds.** The blackbird's m ≥ 5/6 left only 2 IGH and 4 IGL
  tracts here. Permutation gives FDR ≈ 0 from m ≥ 4 in both duck loci, and the
  AID check agrees, so m = 4.

### Finding J: a mistake worth not repeating

The IGH J was first placed at ~157 kb by blasting the **whole** post-V portion of
each transcript against the contig. That is wrong: the post-V region contains
D + J + the entire constant region, and the constant region is far longer, so it
dominates the hits. Those 157 kb hits averaged **293 bp beginning 52 bp after the
V end** — an exon downstream of the real J.

The correct test is to blast **only the first ~70 bp after the V end**, which
isolates D+J. That gives 206 hits at 188.6–189.5 kb with a mean length of 37 bp,
i.e. J-sized, and translation of that region carries the canonical IGHJ
C-terminus `...GSIDLWGHGTEVTVS...` (cf. human IGHJ4 `WGQGTLVTVSS`) on the minus
strand, with a `CACAGTG` heptamer upstream. **J is at ~188,668 (−).**

The blackbird IGL J was checked the same way and is confirmed: 567 junction-only
hits at 6,339 kb, mean length 39 bp, matching the position already in use.

### What reproduced, and what did not

| | duck | blackbird |
|---|---|---|
| IGL V genes with an RSS | **15/51** | 2/23 |
| IGH V genes with an RSS | 7/60 | 25/162 |
| IGL outside-tract hotspot | ×1.12 (n.s.) | ×2.24 *** |
| IGH outside-tract hotspot | ×1.48 *** | ×2.23 *** |
| inside-tract hotspot | 0.83 / 0.83 | 0.61 / 0.54 |

The **direction** reproduces — outside-tract differences hotspot-enriched,
inside-tract ones not — but far more weakly, and duck IGL does not reach
significance. Two contributors worth stating.

First, the duck IGL architecture is genuinely different: 15 candidate parents,
not one.

Second, the **sample provenance**. Checking Mueller *et al.* 2021
(GigaScience 10:giab081) directly: the assembly DNA and the Iso-Seq RNA come from
the *same cohort* — ten captive-bred ducks (5 female, 5 male, 12 months) sampled
in one dissection, with the assembly DNA taken from "lung tissue of a female
tufted duck" among them. So this is **not** the different-individual case, and an
earlier version of this guide said so wrongly.

It is not a clean matched-individual design either. Iso-Seq was run once per
*tissue* (six runs: brain, ileum, lung, ovary, spleen, testis), drawn from that
ten-bird cohort, and both ovary and testis appear, so at least one male and one
female contributed. The methods state only that "technical replicates were
pooled" and never say whether RNA was pooled across individuals within a tissue.
Either way the transcripts originate from up to ten animals and are scored
against one haploid reference.

That is a *different* confound from the blackbird's, and arguably a worse one for
this particular test: between-individual germline polymorphism is not
AID-targeted, so it adds hotspot-neutral differences to the outside-tract class
and dilutes exactly the enrichment being measured. It is a plausible contributor
to duck IGL sitting at ×1.12 where the blackbird reaches ×2.24. Consistent with
the reference-choice supplement, it is invisible in identity: duck IGL sits at
0.947 against the blackbird's 0.951.

**Duck IGH gives 104/104 tracts impossible.** The J was then searched for
exhaustively, because the alternative reading would contradict the standard
post-rearrangement model: a junction-only blast across the whole 294.8 kb contig
gives one cluster (206 hits at 188.6 kb) and nothing above the V array; a
six-frame scan of the entire contig on both strands for the full IGHJ C-terminus
`W-G-x-G-T-x-[VLIAM]-T-V-[SA]` returns **exactly one hit**, at 188,649 (−); and
blasting that J against the whole genome returns only itself. There is no J on
the plus strand and none above the V cluster.

**Assembly error is ruled out.** `bAytFul3.1` (GCA_976225485.1, 2025, a different
individual) places D and J on the same side of the V array. Two independent
assemblies agree, so the J-distal functional gene is real geometry.

**But the tracts look like germline polymorphism, not conversion.** Mapping runs
to tissues via the paper's Table 2 subread counts gives ovary, lung, ileum, brain,
testis, spleen — and ovary/testis are sex-specific, so at least two birds
contributed. The IGH result is essentially **one tissue**: 103 of 104 tract calls
come from the ileum run, which supplies 145 of the 154 IGH transcripts analysed.
And **all 60 distinct IGH tracts are private to a single run**; none recurs in
another tissue. IGL is the same — 28 distinct tracts, 0 shared.

That is the signature of germline variation. Under somatic conversion onto a
shared germline, individual tracts would be transcript-private but the
donor→parent *pairs* would recur across tissues, because every bird carries the
same donor array. Nothing recurring at all fits transcripts from different birds
each carrying their own alleles, scored against the single assembly female.

So this is **not** evidence for conversion before rearrangement. Ranked:
(1) transcripts carry other birds' germline alleles; (2) donor attribution is
wrong, as the topology control concluded for blackbird IGH; (3) genuine
pre-rearrangement conversion — possible, least supported, not separable here. The two expressed RSS-bearing
parents sit at 229,918 and 230,950, the far end of the array from J, so a
deletional rearrangement there removes every donor. Setting `j_strand: "+"` would
make that number look good by reclassifying them as inversional, but the WGxG
evidence does not support it. The honest reading is that duck IGH donor
attribution is untrustworthy — the same verdict the topology control reached for
blackbird IGH.

## Duplicate V gene annotations (same DNA, both strands)

The V gene annotation lists some loci **twice**: two entries whose genomic
intervals overlap, on opposite strands, whose forward-strand sequence is
identical over the shared block. A V gene's reverse complement still scores
against V profiles, so an annotator scanning both strands emits a hit on each.

This is **not** the tandem duplication that fills a V array. Tandem duplicates
sit at different coordinates and are 90-99% similar; these sit at the same
coordinates and are 100% identical to themselves.

| | annotated entries | overlapping pairs | distinct loci |
|---|---|---|---|
| bAgePho2 IGH | 162 | **59** (58/59 identical over the overlap) | **103** |
| bAgePho2 IGL | 23 | 1 | 22 |
| duck IGH | 60 | 26 | 34 |
| duck IGL | 51 | 1 | 50 |

**Which member is real.** Only one member can carry the V reading frame -- the
reverse complement of a V exon is not a V exon. `scripts/gc_dedup_annotations.py`
translates each entry in its own declared orientation and scores the V domain
hallmarks (FR2 `W[VILM]RQ`, `[LIV]EW[VILMA]`, FR3 `Y[YFH]C`). On bAgePho2 IGH the
test is unanimous: in all 59 pairs the loser scores **zero** motifs and the winner
1-3. It agrees with two independent annotations -- the winner carries the RSS in
18 of 18 informative pairs and the transcripts in 23 of 24.

**What collapsing changes: denominators, and nothing else.** Rerunning stage 2 on
the deduplicated FASTA (with the PAFs filtered to the kept targets, which is
equivalent to having aligned against it) reproduces the tract calls **exactly** --
IGH 42 calls / 28 pairs, IGL 123 / 10 -- and every AID value is byte-identical.
What moves is "25 of 162 IGH V genes carry an RSS" becoming 25 of 103, and the
`n possible pairs` denominators in the network panels.

> Nor could a tract call have changed: no donor -> parent pair in either locus is
> an overlap pair, and both members point at the same DNA, so as a conversion
> donor they are one option rather than two.

Config `config_geneconv_bAgePho2_dedup.yaml`, output
`results_dedup/`, report `{locus}_dedup_report.tsv`, figure
`SUPP_IGH_overlapping_annotations`.


## Colour conventions

LaCroix palette throughout, matching the bird paper.

| | |
|---|---|
| IGH | `#87b4dc` |
| IGL | `#638E6E` |
| V / D / J | `#172869` / `#088BBE` / `#1BB6AF` |
| Yes / possible | `#088BBE` |
| No / impossible | `#EA7580` |
| Neutral, background | `#D9D9D9` |
| Heatmaps | `viridis_r` — dark purple = high, yellow = low |

Subfigure labels are `A` `B` with no parentheses. Green-vs-red pairings are
avoided everywhere.

---

## Reproducing

```bash
snakemake -s Snakefile --configfile config_bAgePho2_matched.yaml --cores 16
```

```bash
snakemake -s Snakefile.geneconv --configfile config_geneconv_bAgePho2.yaml --cores 16
```

Stage 1 calls `python` from `PATH` rather than a configured interpreter, so the
conda environment has to be active:

```bash
export PATH=/home/kav67/miniconda3/envs/snakemake/bin:$PATH
```

All figure output goes through `save_figure()` (single figures) or
`MultiPageFigures` (paged reports) in `scripts/gc_palette.py`, so PDF, PNG and SVG
are always written together and a figure can never exist in one format but not
another.
