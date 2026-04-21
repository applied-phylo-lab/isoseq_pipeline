from pathlib import Path

configfile: "config.yaml"

SAMPLES   = config["samples"]
LOCI      = list(config["vgene_loci"].keys())
INPUT_MODE = config["input_mode"]
RESULTS   = config["results_dir"]
TOP_N     = config["top_n_alignments"]

if INPUT_MODE not in ("raw_bam", "preprocessed_fastq"):
    raise ValueError(f"input_mode must be 'raw_bam' or 'preprocessed_fastq', got: {INPUT_MODE!r}")


# ─── Helpers ──────────────────────────────────────────────────────────────────

def preprocessed_fastq(wildcards):
    return config["preprocessed_fastq"].format(sample=wildcards.sample)


# ─── Target ───────────────────────────────────────────────────────────────────

rule all:
    input:
        expand(f"{RESULTS}/{{sample}}/alignment/ig_transcripts.fasta", sample=SAMPLES),
        expand(f"{RESULTS}/{{sample}}/analysis/exact_match_summary.tsv", sample=SAMPLES),
        expand(f"{RESULTS}/{{sample}}/analysis/top{TOP_N}_alignments.tsv", sample=SAMPLES),
        expand(f"{RESULTS}/{{sample}}/analysis/summary_plot.pdf", sample=SAMPLES),


# ══════════════════════════════════════════════════════════════════════════════
#  BRANCH A — raw BAM input: full IsoSeq preprocessing
# ══════════════════════════════════════════════════════════════════════════════

rule lima:
    input:
        bam=config["hifi_bam"],
        primers=config["primers_fasta"],
    output:
        bam=f"{RESULTS}/{{sample}}/isoseq/fl.bam",
    log:
        f"{RESULTS}/logs/{{sample}}/lima.log",
    shell:
        "lima {input.bam} {input.primers} {output.bam} "
        "--isoseq --peek-guess 2>{log}"


rule isoseq_refine:
    input:
        bam=f"{RESULTS}/{{sample}}/isoseq/fl.bam",
        primers=config["primers_fasta"],
    output:
        bam=f"{RESULTS}/{{sample}}/isoseq/flnc.bam",
    params:
        polya="--require-polya" if config["require_polya"] else "",
    log:
        f"{RESULTS}/logs/{{sample}}/refine.log",
    shell:
        "isoseq refine {input.bam} {input.primers} {output.bam} "
        "{params.polya} 2>{log}"


rule make_fofn:
    input:
        f"{RESULTS}/{{sample}}/isoseq/flnc.bam",
    output:
        f"{RESULTS}/{{sample}}/isoseq/flnc.fofn",
    shell:
        "realpath {input} > {output}"


rule isoseq_cluster:
    input:
        fofn=f"{RESULTS}/{{sample}}/isoseq/flnc.fofn",
    output:
        bam=f"{RESULTS}/{{sample}}/isoseq/clustered.bam",
    params:
        singletons="--singletons" if config["include_singletons"] else "",
        cmd="cluster2" if config["use_cluster2"] else "cluster",
        extra="" if config["use_cluster2"] else "--verbose --use-qvs",
    log:
        f"{RESULTS}/logs/{{sample}}/cluster.log",
    shell:
        "isoseq {params.cmd} {input.fofn} {output.bam} "
        "{params.singletons} {params.extra} 2>{log}"


rule extract_clustered_fasta:
    input:
        f"{RESULTS}/{{sample}}/isoseq/clustered.bam",
    output:
        f"{RESULTS}/{{sample}}/isoseq/clustered.fasta",
    log:
        f"{RESULTS}/logs/{{sample}}/extract_fasta.log",
    shell:
        "samtools fasta {input} > {output} 2>{log}"


# ══════════════════════════════════════════════════════════════════════════════
#  BRANCH B — preprocessed FASTQ input: convert to FASTA only
# ══════════════════════════════════════════════════════════════════════════════

rule fastq_to_fasta:
    input:
        preprocessed_fastq,
    output:
        f"{RESULTS}/{{sample}}/isoseq/clustered.fasta",
    log:
        f"{RESULTS}/logs/{{sample}}/fastq_to_fasta.log",
    shell:
        "awk 'NR%4==1{{sub(\"^@\",\">\"); print}} NR%4==2{{print}}' "
        "{input} > {output} 2>{log}"


# ══════════════════════════════════════════════════════════════════════════════
#  Snakemake rule priority: only activate the branch that matches input_mode
# ══════════════════════════════════════════════════════════════════════════════

if INPUT_MODE == "raw_bam":
    ruleorder: extract_clustered_fasta > fastq_to_fasta
else:
    ruleorder: fastq_to_fasta > extract_clustered_fasta


# ══════════════════════════════════════════════════════════════════════════════
#  Shared downstream: alignment, filtering, analysis
# ══════════════════════════════════════════════════════════════════════════════

# ─── Step 5: Align all transcripts to V gene loci (screening pass) ────────────

rule minimap2_align_filter:
    input:
        transcripts=f"{RESULTS}/{{sample}}/isoseq/clustered.fasta",
        vgenes=lambda wc: config["vgene_loci"][wc.locus],
    output:
        paf=f"{RESULTS}/{{sample}}/alignment/{{locus}}.paf",
    params:
        preset=config["minimap2_preset"],
    log:
        f"{RESULTS}/logs/{{sample}}/minimap2_{{locus}}.log",
    threads: 8
    shell:
        "minimap2 -cx {params.preset} --cs -t {threads} "
        "{input.vgenes} {input.transcripts} "
        "> {output.paf} 2>{log}"


# ─── Step 6: Filter — keep only IG transcripts ────────────────────────────────

rule filter_ig_transcripts:
    input:
        transcripts=f"{RESULTS}/{{sample}}/isoseq/clustered.fasta",
        pafs=expand(f"{RESULTS}/{{sample}}/alignment/{{locus}}.paf",
                    locus=LOCI, allow_missing=True),
    output:
        fasta=f"{RESULTS}/{{sample}}/alignment/ig_transcripts.fasta",
        ids=f"{RESULTS}/{{sample}}/alignment/ig_transcript_ids.txt",
        stats=f"{RESULTS}/{{sample}}/alignment/filter_stats.tsv",
    params:
        min_identity=config["min_alignment_identity"],
        min_coverage=config["min_alignment_coverage"],
    log:
        f"{RESULTS}/logs/{{sample}}/filter_ig.log",
    shell:
        "python scripts/filter_ig_transcripts.py "
        "--transcripts {input.transcripts} "
        "--pafs {input.pafs} "
        "--min-identity {params.min_identity} "
        "--min-coverage {params.min_coverage} "
        "--out-fasta {output.fasta} "
        "--out-ids {output.ids} "
        "--out-stats {output.stats} "
        "2>{log}"


# ─── Step 7: Detailed alignment with secondary hits (IG transcripts only) ─────

rule minimap2_align_detailed:
    input:
        transcripts=f"{RESULTS}/{{sample}}/alignment/ig_transcripts.fasta",
        vgenes=lambda wc: config["vgene_loci"][wc.locus],
    output:
        paf=f"{RESULTS}/{{sample}}/alignment/{{locus}}_detailed.paf",
    params:
        preset=config["minimap2_preset"],
        top_n=TOP_N,
    log:
        f"{RESULTS}/logs/{{sample}}/minimap2_{{locus}}_detailed.log",
    threads: 8
    shell:
        "minimap2 -cx {params.preset} --cs "
        "-N {params.top_n} -p 0 "
        "-t {threads} "
        "{input.vgenes} {input.transcripts} "
        "> {output.paf} 2>{log}"


# ─── Step 8: Exact match analysis ─────────────────────────────────────────────

rule analyze_exact_matches:
    input:
        pafs=expand(f"{RESULTS}/{{sample}}/alignment/{{locus}}_detailed.paf",
                    locus=LOCI, allow_missing=True),
        vgene_fastas=list(config["vgene_loci"].values()),
    output:
        summary=f"{RESULTS}/{{sample}}/analysis/exact_match_summary.tsv",
        per_transcript=f"{RESULTS}/{{sample}}/analysis/exact_match_per_transcript.tsv",
    log:
        f"{RESULTS}/logs/{{sample}}/exact_matches.log",
    shell:
        "python scripts/analyze_exact_matches.py "
        "--pafs {input.pafs} "
        "--vgene-fastas {input.vgene_fastas} "
        "--out-summary {output.summary} "
        "--out-per-transcript {output.per_transcript} "
        "2>{log}"


# ─── Step 9: Top-N alignment table (gene conversion quantification) ────────────

rule analyze_top_alignments:
    input:
        pafs=expand(f"{RESULTS}/{{sample}}/alignment/{{locus}}_detailed.paf",
                    locus=LOCI, allow_missing=True),
    output:
        table=f"{RESULTS}/{{sample}}/analysis/top{TOP_N}_alignments.tsv",
    params:
        top_n=TOP_N,
    log:
        f"{RESULTS}/logs/{{sample}}/top_alignments.log",
    shell:
        "python scripts/analyze_top_alignments.py "
        "--pafs {input.pafs} "
        "--top-n {params.top_n} "
        "--output {output.table} "
        "2>{log}"


# ─── Step 10: Summary plots ────────────────────────────────────────────────────

rule plot_summary:
    input:
        exact=f"{RESULTS}/{{sample}}/analysis/exact_match_summary.tsv",
        top_aln=f"{RESULTS}/{{sample}}/analysis/top{TOP_N}_alignments.tsv",
    output:
        f"{RESULTS}/{{sample}}/analysis/summary_plot.pdf",
    log:
        f"{RESULTS}/logs/{{sample}}/plot.log",
    shell:
        "python scripts/plot_summary.py "
        "--exact-matches {input.exact} "
        "--top-alignments {input.top_aln} "
        "--output {output} "
        "2>{log}"
