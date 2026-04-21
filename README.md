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

### Step 10 — Summary plots (`plot_summary`)

`scripts/plot_summary.py` produces a five-panel PDF (`summary_plot.pdf`):

1. **Exact match counts per gene** — bar chart coloured by locus; any bar > 0 is a functional candidate
2. **Identity distribution rank-1 vs rank-2** — overlapping histograms
3. **Delta-identity histogram** — distribution of rank1 − rank2 identity; low values suggest gene conversion
4. **V gene coverage distribution** — fraction of each V gene covered by the best alignment
5. **Per-transcript identity decay** — identity across ranks 1–N for a random sample of transcripts (median overlaid in red)

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
    ├── alignment/
    │   ├── {locus}.paf            # screening alignments
    │   ├── {locus}_detailed.paf   # multi-hit alignments
    │   ├── ig_transcripts.fasta   # filtered IG transcripts
    │   ├── ig_transcript_ids.txt
    │   └── filter_stats.tsv
    └── analysis/
        ├── exact_match_summary.tsv
        ├── exact_match_per_transcript.tsv
        ├── top{N}_alignments.tsv
        └── summary_plot.pdf
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
