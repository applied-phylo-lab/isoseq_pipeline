# IsoSeq IG Transcript Pipeline

A Snakemake pipeline for identifying and analysing avian immunoglobulin (IG) transcripts from PacBio IsoSeq long-read data.  The pipeline handles two input modes (raw BAM or pre-processed FASTQ), filters transcripts against germline V gene databases, detects functionally rearranged loci via exact-match analysis, and quantifies gene conversion signatures across the top-N V gene alignments per transcript.

---

## Requirements

| Tool | Version | Conda env |
|---|---|---|
| Snakemake | ≥ 7.0 | `alignment_env` |
| minimap2 | ≥ 2.26 | `alignment_env` |
| pandas / numpy / matplotlib | — | `alignment_env` |
| lima | any | `isoseq_env` (raw BAM mode only) |
| isoseq3 | any | `isoseq_env` (raw BAM mode only) |
| samtools | any | `isoseq_env` (raw BAM mode only) |

Create environments from the provided YAML files:

```bash
conda env create -f envs/alignment.yaml
conda env create -f envs/isoseq.yaml
```

---

## Quick start

```bash
conda activate alignment_env

# edit config.yaml to match your data, then:
snakemake --cores 16 -n        # dry run
snakemake --cores 16           # run
```

---

## Configuration (`config.yaml`)

| Key | Description |
|---|---|
| `results_dir` | Root output directory for all results and logs |
| `samples` | List of sample IDs |
| `input_mode` | `"raw_bam"` or `"preprocessed_fastq"` (see below) |
| `kinnex` | `true` for Kinnex / MAS-Seq libraries — inserts a `skera split` step before lima to break each concatenated HiFi read into its cDNA segments. Default `false` (plain Iso-Seq) |
| `mas_adapters_fasta` | MAS array adapter FASTA for `skera split`, required when `kinnex: true` (e.g. `data/mas8_primers.fasta`) |
| `hifi_bam` | Path pattern for raw HiFi BAM files (`{sample}` placeholder) |
| `primers_fasta` | Primers FASTA used by lima and isoseq refine |
| `preprocessed_fastq` | Path pattern for pre-processed FASTQ files (`{sample}` placeholder) |
| `vgene_loci` | Dict mapping locus name → FASTA of germline V genes for that locus |
| `require_polya` | Whether to pass `--require-polya` to isoseq refine |
| `use_cluster2` | Use `cluster2` (recommended) vs legacy `cluster` |
| `include_singletons` | Keep singleton reads from clustering |
| `min_alignment_identity` | Minimum identity (0–1) to call a transcript IG-related |
| `min_alignment_coverage` | Minimum fraction of the V gene that must be covered |
| `top_n_alignments` | Number of secondary alignments to retain per transcript |
| `run_immunotools` | Run immunotools `diversity_analyzer` alongside the minimap2 filter (default `true`). Set `false` for species with no immunotools germline database — the merge step then falls back to the minimap2 filter alone |
| `immunotools_org` | Organism name passed to `diversity_analyzer --org` (default `tufted_duck`). Build a new one with `scripts/make_immunotools_org.sh` |
| `immunotools_path` | Path to `diversity_analyzer.py` |
| `minimap2_preset` | minimap2 preset: `splice:hq` (IsoSeq/cDNA, default) or `map-ont` (ONT) — do **not** use `map-hifi`, it is designed for genomic mapping and its k=19 / min-score=80 settings will produce no alignments against short V gene references |

---

## Input modes

### `raw_bam`

Start from raw PacBio HiFi BAM files.  The pipeline runs the full IsoSeq preprocessing stack automatically:

```
hifi_reads.bam
    │
    ▼ lima  (primer removal + demultiplexing)
fl.bam
    │
    ▼ isoseq refine  (poly-A trimming, concatemer removal)
flnc.bam
    │
    ▼ isoseq cluster2  (isoform clustering)
clustered.bam
    │
    ▼ samtools fasta
clustered.fasta  ──► downstream analysis
```

### `preprocessed_fastq`

Start from FASTQ files that have already been through IsoSeq preprocessing (e.g. exported from SMRT Link or a previous run).  The pipeline converts them to FASTA and skips directly to alignment:

```
{sample}.fastq
    │
    ▼ awk (FASTQ → FASTA)
clustered.fasta  ──► downstream analysis
```

---

## Pipeline steps

### Step 5 — Screening alignment (`minimap2_align_filter`)

All clustered transcripts are aligned to each V gene locus FASTA independently using minimap2.  One PAF file is produced per sample × locus.  The `--cs` flag adds a compressed difference string to every alignment, enabling exact-match detection later.

```
minimap2 -cx {preset} --cs {locus_vgenes.fasta} {clustered.fasta} > {locus}.paf
```

### Step 6 — IG transcript filtering (`filter_ig_transcripts`)

`scripts/filter_ig_transcripts.py` reads all screening PAF files and keeps only transcripts that meet **both** thresholds against at least one V gene in any locus:

- identity ≥ `min_alignment_identity`
- V gene coverage ≥ `min_alignment_coverage`

Outputs:
- `ig_transcripts.fasta` — filtered transcript sequences
- `ig_transcript_ids.txt` — transcript IDs that passed
- `filter_stats.tsv` — total vs retained counts, breakdown by best locus

### Step 7 — Detailed alignment (`minimap2_align_detailed`)

The filtered IG transcripts are re-aligned to each locus with secondary alignments enabled (`-N {top_n} -p 0`), retaining up to `top_n_alignments` hits per transcript.  This is the alignment used for all downstream analysis.

### Step 8 — Exact match analysis (`analyze_exact_matches`)

`scripts/analyze_exact_matches.py` identifies cases where a V gene is found **verbatim** inside a transcript — no mismatches, no indels, full gene coverage.  The criterion is:

```
alignment_block_length == V_gene_length  AND  n_matching_bases == V_gene_length
```

A gene meeting this criterion is a **functional VDJ recombination candidate**: the rearranged locus contains that exact germline sequence.

Outputs:
- `exact_match_summary.tsv` — per gene: transcript count, fraction of IG transcripts, functional flag
- `exact_match_per_transcript.tsv` — per transcript: which genes were found exactly

### Step 9 — Top-N alignment table (`analyze_top_alignments`)

`scripts/analyze_top_alignments.py` collects the top-N V gene alignments for every IG transcript and computes per-alignment metrics:

| Column | Meaning |
|---|---|
| `rank` | 1 = best hit by identity |
| `identity` | Fraction of aligned bases that match |
| `v_coverage` | Fraction of the V gene covered by this alignment |
| `query_frac` | Fraction of the transcript covered |
| `delta_identity_from_rank1` | Identity drop from the best hit |
| `is_exact` | True if this alignment is a perfect verbatim match |

**Interpreting gene conversion:** birds use gene conversion to diversify IG sequences — donor pseudogenes overwrite segments of the functional V gene.  A transcript derived from gene conversion will show a small `delta_identity_from_rank1` between rank-1 and rank-2 hits because multiple V genes explain the sequence nearly equally well.  A large delta indicates a single dominant donor.

### Step 10 — Per-sample summary plots (`plot_summary`)

`scripts/plot_summary.py` produces a five-panel PDF (`summary_plot.pdf`):

1. **Exact match counts per gene** — bar chart coloured by locus; any bar > 0 is a functional candidate
2. **Identity distribution rank-1 vs rank-2** — overlapping histograms
3. **Delta-identity histogram** — distribution of rank1 − rank2 identity; low values suggest gene conversion
4. **V gene coverage distribution** — fraction of each V gene covered by the best alignment
5. **Per-transcript identity decay** — identity across ranks 1–N for a random sample of transcripts (median overlaid in red)

### Step 11 — Combined cross-sample plot (`plot_combined`)

`scripts/plot_combined.py` pools results from all samples and produces `combined/combined_summary_plot.pdf` with five panels:

1. **IGH V gene usage heatmap** — genes (rows) × samples (columns), colour = transcript count; top 20 genes shown. Immediately shows which V genes dominate across experiments.
2. **IGL V gene usage heatmap** — same layout for IGL.
3. **Locus assignment stacked bar** — per sample, how many IG transcripts were assigned IGH-only, IGL-only, both loci, or other. Shows the overall composition of each library.
4. **Co-occurring transcripts bar** — count of transcripts with hits to **both** IGH and IGL V genes within the top-N alignments. A transcript showing strong hits to both loci is a particularly reliable antibody transcript candidate, since it is highly unlikely for a non-IG read to align well to genes from two separate loci simultaneously.
5. **IGH × IGL combination heatmap** — for co-occurring transcripts, which IGH gene is paired with which IGL gene (pooled across samples). Identical combinations appearing across multiple samples suggest a dominant functional VDJ + VJ pairing in the population.

---

## Output structure

```
{results_dir}/
├── logs/
│   └── {sample}/
│       ├── fastq_to_fasta.log
│       ├── minimap2_{locus}.log
│       ├── filter_ig.log
│       ├── minimap2_{locus}_detailed.log
│       ├── exact_matches.log
│       ├── top_alignments.log
│       └── plot.log
└── {sample}/
    ├── isoseq/
    │   └── clustered.fasta        # input to alignment steps
    ├── immunotools/
    │   └── cleaned_sequences.fasta          # diversity_analyzer output
    ├── alignment/
    │   ├── {locus}.paf                      # screening alignments (minimap2)
    │   ├── {locus}_detailed.paf             # multi-hit alignments (combined set)
    │   ├── ig_transcripts.fasta             # minimap2-filtered transcripts
    │   ├── ig_transcript_ids.txt
    │   ├── filter_stats.tsv
    │   ├── combined_ig_transcripts.fasta    # union of minimap2 + immunotools
    │   └── merge_stats.tsv                  # overlap counts between the two methods
    └── analysis/
        ├── exact_match_summary.tsv
        ├── exact_match_per_transcript.tsv
        ├── top{N}_alignments.tsv
        └── summary_plot.pdf
└── combined/
    └── combined_summary_plot.pdf  # cross-sample combined analysis
```

---

## V gene databases

Place one FASTA per locus under `data/vgenes/` and register them in `config.yaml`:

```yaml
vgene_loci:
  IGH: "data/vgenes/tufted_duck_IGH.fasta"
  IGL: "data/vgenes/tufted_duck_IGL.fasta"
```

Each sequence header should be a unique gene identifier (e.g. `>IGHV1-1`).  Any number of loci is supported.

---

# Stage 2 — Gene conversion analysis

```bash
snakemake -s Snakefile.geneconv --configfile config_geneconv_bAgePho2.yaml --cores 16
```

Stage 1 turns reads into IG transcripts and says which V gene each one came
from. Stage 2 asks a different question: **where did the rest of the sequence
come from?** In birds the answer is usually gene conversion — a rearranged V
gene is repeatedly overwritten in patches by nearby pseudogenes.

The hard part is that gene conversion and somatic hypermutation (SHM) produce
the same raw observation: a transcript that differs from its germline gene. The
whole of stage 2 is about separating those two, and about not fooling ourselves
while doing it. Read the statistics section below before quoting any number.

---

## Step-by-step

### 1. RSS annotation (`gc_rss_annotation.py`)

A V gene can only be *rearranged* if it carries a recombination signal sequence.
Without one it can still **donate** sequence, but it can never be the gene a
transcript came from. RSS calls are taken as given from the curated screen
(`clean_birds/gene_list.csv`) — binary, present or absent, no re-calling.

|  | V genes | with RSS |
|---|---|---|
| IGH | 162 | 25 |
| IGL | 23 | 2 |

Few functional genes is expected in a gene-conversion species; it is the
situation in chicken and duck, and it is the reason conversion exists.

### 2. Which genes are actually used (`gc_call_functional_genes.py`)

RSS is a property of the DNA; expression is a property of the animal. A gene is
**expression-confirmed** when ≥2 transcripts match it at ≥98% identity over
≥200 bp, scoring with the last 20 bp excluded (V(D)J junction formation destroys
germline identity there regardless of function).

Two assignment tables are produced and they are *not* interchangeable:

| File | Eligible parents | Used for |
|---|---|---|
| `transcript_assignments.tsv` | RSS-bearing genes only | tract detection, AID test |
| `unconstrained_assignments.tsv` | every V gene | expression, RSS validation |

Using the constrained table to test the RSS annotation would be circular.

### 3. Donor pool by locus topology (`gc_donor_pool.py`)

Conversion here is taken to be **cis-acting**, so which donors physically
existed depends on how the gene rearranged:

* **Deletion** (V co-oriented with J) — everything between V and J is looped out.
  Donors in that interval are **gone**.
* **Inversion** (V opposite to J) — the DNA is flipped, not removed. All donors
  survive.

This produces the analysis's most valuable asset: a tract whose donor was
already deleted is a **known false positive**, so the error rate can be measured
rather than assumed.

### 4. Tract detection (`gc_detect_tracts.py`)

For each transcript, against its assigned parent:

1. find every position where transcript ≠ parent;
2. ask whether a donor carries exactly that base there (*donor-diagnostic*);
3. require supporting positions to be **contiguous within 5 bp** — a conversion
   tract is a continuous stretch of donor sequence, not scattered similarity;
4. score it (see the statistics section).

The contiguity rule does most of the work. Without it a "tract" degenerates into
"this window resembles donor D overall", which is a statement about gene
similarity rather than about conversion.

### 5. BrepConvert (`run_brepconvert.R`)

The published tool, run for comparison. Three fixes were needed: its IMGT-style
name regex returned `NA` for our gene names (we emit IMGT-shaped aliases);
Biostrings ≥2.77 moved `pairwiseAlignment`/`pid`/`stringDist` into `pwalign`;
and `/` in transcript IDs broke its temp FASTA writer.

### 6–8. AID spectrum, D usage, figures

Covered in the statistics section and in `docs/PIPELINE_GUIDE.md`.

---

## The statistical tests, in detail

Six independent tests are run. They are listed here with what each one assumes,
what it can prove, and — importantly — what it cannot.

### Test 1 — Is a tract more than chance? (analytic, **superseded**)

*Statistic* `m` = number of donor-diagnostic positions supporting the tract.
*Null* differences are independent point mutations; given that a position
differs, the chance it happens to carry the donor's specific base is ~1/3, so

```
p_raw = 3^-m          p_corrected = n_donors × 3^-m      (Bonferroni)
```

**Why it was abandoned.** It corrects for the wrong thing. `3^-m` models the
chance of *matching a donor base*, but says nothing about **contiguity**, which
is the constraint actually doing the work. The result is wildly conservative:
with 162 IGH donors it demands m ≥ 8 and returned 10 tracts in the whole locus.
Benjamini–Hochberg is no help either — with `min_informative = 3` every raw
p-value is ≤ 0.037 by construction, so BH accepts 100% of candidates.

### Test 2 — Permutation null for tract calling (**this sets the floor**)

Each transcript keeps its exact differences from its parent — same number, same
substituted bases — but they are scattered to random positions along the gene.
Contiguity is destroyed; everything else is preserved. Any tract found in
permuted data is false by construction, so

```
FDR(m) = (null tracts per replicate at ≥ m) / (real tracts at ≥ m)
```

20 replicates:

| m ≥ | IGL real | IGL null | IGH real | IGH null |
|---|---|---|---|---|
| 3 | 1776 | 0.80 | 1959 | 5.60 |
| 4 | 556 | **0.00** | 513 | **0.00** |
| 5 | 123 | 0.00 | 99 | 0.00 |

So contiguity alone is very hard to achieve by chance, and m ≥ 4 already has an
FDR indistinguishable from zero **against scattered mutation**.

### Test 3 — Topology control (a *direct* false-positive measurement)

If a fraction `e` of the donor pool was deleted for a given parent, false calls
land on a deleted donor with probability `e` while true calls never do. With
observed impossible fraction `f`:

```
FDR = f / e
```

No distributional assumption at all. On **IGH deletional parents** this gives
FDR ≈ 1.0 at every threshold — the donor *attribution* there is essentially
random. In **IGL** the single parent sits at the J-proximal end so nothing is
ever deleted (`e = 0`) and this test has **no power**; that is why test 2 exists.

### Test 4 — AID hotspot / coldspot spectrum (**the primary evidence**)

AID initiates both processes but only SHM leaves its fingerprint on the result:
SHM mutates *at* the AID-targeted base, conversion copies a donor and the
differences land wherever the donor happened to differ.

*Null* each transcript's mutations are redistributed at random over the
positions actually covered in that same gene, keeping the count fixed
(1000 permutations). This controls for base composition, which otherwise makes
raw motif percentages meaningless.

Prediction, and result at the calibrated thresholds:

| | hotspot (WRCY/RGYW) | coldspot (SYC/GRS) |
|---|---|---|
| **IGL outside tracts** | ×2.23 (p = 0.001) ↑ | ×0.40 (p = 0.001) ↓ |
| **IGL inside tracts** | ×0.48 (p = 0.001) ↓ | ×1.56 (p = 0.001) ↑ |
| **IGH outside tracts** | ×2.18 (p = 0.002) ↑ | ×0.51 (p = 0.002) ↓ |
| **IGH inside tracts** | ×0.91 (n.s.) | ×0.69 (p = 0.084) |

Outside-tract differences carry a textbook AID signature — enriched at hotspots
**and** depleted at coldspots. Inside-tract differences do not; in IGL they run
the other way entirely. The coldspot arm is what makes this decisive: a process
that merely clustered mutations would move hotspots and coldspots the *same*
direction, and only AID pushes them in opposite directions.

*Circularity check:* donor-difference positions are hotspot-neutral (10.7% vs
10.2% background), so the detector is not simply selecting non-hotspot positions.

### Test 5 — Transition bias and C:G targeting (**not evidence — see below**)

Reported for completeness. Neither supports the hypothesis; see "What the
spectrum panels do and do not show".

### Test 6 — Topology / evolutionary hypothesis (**vacuous, honestly reported**)

Mann–Whitney U comparing donor availability for functional vs non-functional
genes. It returns nothing usable, because in IGL the parent has **zero** deleted
donors and in IGH every tract-bearing parent is inversional. There is no
contrast to test.

That is not a failed analysis — it is the hypothesis confirmed from the other
side. Functional genes sit exactly where the test has no power, which is what
the hypothesis predicts. It cannot be reported as a p-value, and
`IGx_report_*.pdf` panel C says so on its face.

### Test 7 — D gene usage (negative result)

Longest shared stretch between the junction and each D gene, against a
shuffled-junction null. Median junction is 32 bp; observed match 7 bp vs 6.5 bp
by chance; only 4/70 transcripts give a confident call and those average 8 tied
D genes. **D usage is not determinable from this data.** The figure is included
because the answer is a real result, not because the numbers are usable.

### Test 8 — Germline reference comparison

Same transcripts scored against each reference; genes paired by **best reciprocal
alignment** so "different gene" means genuinely different, not differently named.
The control is `bAgePho2_alt`, the same bird's other haplotype — see
`docs/PIPELINE_GUIDE.md` for why that is a haploid-reference floor rather than a
technical noise floor.

---

## How the tract threshold was chosen

Tests 2 and 3 constrain from below and above, but neither is sufficient alone.
The permutation null (test 2) scatters mutations, so it cannot detect
**clustered SHM** being admitted as tracts. The AID spectrum supplies that
orthogonal check: if the threshold is too permissive, inside-tract differences
start carrying an AID signature they should not have.

Sweeping the support cutoff:

| m ≥ | IGL tracts | IGL inside hotspot | IGH tracts | IGH inside hotspot |
|---|---|---|---|---|
| 4 | 556 | **1.26** ✗ | 513 | **1.45** ✗ |
| 5 | 123 | 0.48 ✓ | 99 | 1.07 ✗ |
| 6 | 93 | 0.40 ✓ | 42 | 0.88 ✓ |
| 7 | 36 | 0.64 ✓ | 15 | 0.83 ✓ |
| 8 | 19 | 0.00 ✓ | 10 | 0.84 ✓ |

**Chosen: m ≥ 5 for IGL, m ≥ 6 for IGH** — the point at which inside-tract
differences stop looking AID-driven. The loci differ because their donor pools
differ (23 vs 162 genes), which is the same logic Bonferroni was reaching for,
calibrated empirically instead of analytically.

Relative to the old Bonferroni rule this is a real relaxation: IGH goes from
**10 → 42** tracts and IGL from **93 → 123**, with no loss of specificity by any
control available.

*On circularity:* the AID contrast is stable across m = 5…8 (outside stays
2.10–2.23 throughout, inside stays below 1), so the qualitative result does not
depend on where in that range the line is drawn.

---

## What the spectrum panels do and do not show

`IGx_aid_spectrum.pdf` has four panels. **Only panel A is evidence.**

**Panel A — hotspot/coldspot.** The real test. Composition-matched permutation
null, two directions, both loci, all four cases. See test 4.

**Panel B — transition bias.** Does *not* support the hypothesis, and the two
loci disagree:

| | outside | inside | prediction |
|---|---|---|---|
| IGL | 0.77 | 1.28 | outside should be higher — **it is lower** |
| IGH | 0.89 | 0.80 | outside higher — weakly consistent |

The reason is a confound. Conversion copies differences that accumulated between
**paralogues over evolutionary time**, and molecular evolution is itself strongly
transition-biased. Inside-tract differences therefore inherit an ancient
transition bias that has nothing to do with AID. Ti/Tv cannot separate the two
processes here.

**Panel C — C:G targeting.** Also not supporting. The old version drew a
hardcoded 0.5 baseline, which was simply wrong — the correct expectation is the
C+G content of the covered positions (0.59 in IGL, 0.64 in IGH). Against the
proper baseline:

| | observed | expected | ratio |
|---|---|---|---|
| IGL outside | 0.495 | 0.586 | **0.85** |
| IGL inside | 0.558 | 0.586 | 0.95 |
| IGH outside | 0.621 | 0.641 | **0.97** |
| IGH inside | 0.546 | 0.644 | 0.85 |

No C:G enrichment anywhere — mild depletion if anything. This is not
contradictory with panel A: hotspot motifs are C/G-centred by definition, but
they account for only ~20% of outside-tract differences, so a 2.2× enrichment
there barely moves the *global* C:G fraction.

**Recommendation:** keep panel A in the main figure; move B and C to
supplementary with the confound stated, or drop them. Presenting them as
supporting evidence would not survive review.

---

## Gene conversion outputs

```
{results_dir}/geneconv/
├── rss_annotation.tsv                 # RSS present/absent per V gene
├── functional_genes.tsv               # RSS ∩ expression
├── transcript_assignments.tsv         # parent per transcript (RSS-constrained)
├── unconstrained_assignments.tsv      # best gene per transcript (all genes)
├── {locus}_donor_pool.tsv             # surviving donors per rearranged gene
├── {locus}_tracts.tsv                 # tracts + q-values + donor ambiguity
├── {locus}_tracts_brepconvert.tsv
├── {locus}_aid_spectrum.{pdf,tsv}     # the primary evidence figure
├── {locus}_arrows_{detector,brepconvert}.pdf
├── {locus}_report_{detector,brepconvert}.pdf
├── {locus}_tract_alignments.pdf       # per-event sequence evidence
├── RSS_validation.{pdf,tsv}           # RSS calls vs expression
├── IG_overview.pdf                    # start here
├── SUPP_germline_matched_vs_not_{IGH,IGL}.pdf
└── SUPP_reference_choice_across_haplotypes.pdf
```

Every figure is written as `.pdf`, `.png` and `.svg`. Multi-page reports get one
SVG/PNG per page (`…_p01.svg`).
