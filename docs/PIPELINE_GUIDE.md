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
4. Score by likelihood ratio with a Bonferroni correction over the donor pool:
   `p_corrected = n_donors × 3⁻ᵐ`, where **m is the number of donor-diagnostic
   positions supporting that tract**. Each such position is one of 3 possible
   non-parent bases, so m positions agreeing by chance costs 3⁻ᵐ; multiplying by
   the number of donors tested corrects for having looked at all of them.

The contiguity rule is what dropped the counts sharply (e.g. 3,646 → 116 in the
VGP IGL run). Without it the "tracts" had a median span of 53 bp carrying only 7
supporting positions — scattered agreement, not tracts.

### 11. The AID hotspot spectrum test — the strongest evidence

This is the test that actually separates the two processes, and it does not depend
on the donor pool at all.

**Somatic hypermutation** is done by AID, which is not random: it targets WRCY /
RGYW motifs (hotspots) and avoids SYC / GRS (coldspots). **Gene conversion** copies
a block of donor sequence — the resulting differences sit wherever the donor
happened to differ, with no relationship to AID motifs.

So: differences *outside* tracts should be hotspot-enriched; differences *inside*
tracts should not be.

| | outside tracts | inside tracts |
|---|---|---|
| **IGL** | **2.20×** enriched (p = 0.001) | 0.40× (n.s.) |
| **IGH** | **2.12×** enriched (p = 0.001) | 0.83× (n.s.) |

Both loci, in the predicted direction, against a permutation null that shuffles
positions within each gene. **Circularity check:** if the tract detector were
simply picking non-hotspot positions by construction, this would be an artifact.
It isn't — donor-difference positions are hotspot-neutral (10.7% vs 10.2%
background).

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

### `IGx_arrows_*.pdf` — donor → recipient

The figure asked for: genes laid out along the locus, one arrow per donor→recipient
relationship, labelled with the number of supporting transcripts.

- **Arrow direction**: donor → recipient.
- **Arrow thickness and its number**: transcripts supporting that pair. The number
  now sits directly on its own arc (see note below).
- **Above the axis / teal**: topologically possible — the donor still existed.
- **Below the axis / rose**: **impossible** — the donor had already been deleted.
  This is the visible false-positive load.
- **Marker colour**: navy = has RSS, grey = none.
- **Marker shape**: points *toward* J = deletional, *away* = inversional.
- **Marker size**: large = at least one transcript best-matches it.
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

| Panel | What it shows |
|---|---|
| **A** | AID hotspot targeting, inside vs outside tracts, against the permutation null |
| **B** | Transition bias |
| **C** | C:G targeting |
| **D** | Full substitution spectrum |

Panel A is the result: outside-tract differences are hotspot-enriched (AID/SHM),
inside-tract differences are not (copied from a donor).

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
| bAgePho1_pri | IGH | 74.3% | +1.32% |
| bAgePho1_alt | IGH | 52.2% | +0.88% |
| bAgePho0_pri | IGH | 42.2% | +1.24% |
| bAgePho1_alt | IGL | 31.1% | +0.00% |
| bAgePho0_pri | IGL | 30.3% | +0.00% |
| bAgePho1_pri | IGL | 29.9% | +0.15% |
| **bAgePho2_alt** | IGL | **21.5%** | +0.00% | ← same bird, other haplotype |

The result to quote: in IGL the haploid-reference floor is 21.5%, and the
different-bird references sit at 29.9–31.1% — roughly **9 points of excess
attributable specifically to using a different animal**, on top of a floor that is
itself partly real biology. And the third panel is the warning:
median identity gain is **0.00%** for IGL. No summary statistic would have told
you anything was wrong. The sequences look fine; the *assignments* are wrong.

The third panel encodes locus with **colour only** (IGH blue, IGL green); marker
shape says whether the reference is the same bird (diamond) or a different one
(circle). An earlier version used both colour *and* shape for locus and then
overrode the colour on the control point, which put a blue square on the plot that
matched nothing in the legend.

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
