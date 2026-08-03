from pathlib import Path

configfile: "config.yaml"

SAMPLES   = config["samples"]
LOCI      = list(config["vgene_loci"].keys())
INPUT_MODE = config["input_mode"]
RESULTS   = config["results_dir"]
TOP_N     = config["top_n_alignments"]

# Kinnex / MAS-Seq libraries concatenate several cDNA molecules into one HiFi
# read, so they need a `skera split` pass before lima.  Plain Iso-Seq libraries
# do not.
KINNEX = config.get("kinnex", False)

# Screen FLNC reads against the V gene databases and cluster only the IG hits.
# See the rule block below for why this exists.
PREFILTER_IG = config.get("prefilter_ig_before_cluster", False)

# immunotools/diversity_analyzer needs a germline database (V *and* J genes) for
# the organism.  Species where only V genes are available can switch it off and
# run on the minimap2 filter alone.
RUN_IMMUNOTOOLS = config.get("run_immunotools", True)
IMMUNOTOOLS_ORG = config.get("immunotools_org", "tufted_duck")

if INPUT_MODE not in ("raw_bam", "preprocessed_fastq"):
    raise ValueError(f"input_mode must be 'raw_bam' or 'preprocessed_fastq', got: {INPUT_MODE!r}")


# ─── Helpers ──────────────────────────────────────────────────────────────────

def preprocessed_fastq(wildcards):
    return config["preprocessed_fastq"].format(sample=wildcards.sample)


def hifi_bam(wildcards):
    return config["hifi_bam"].format(sample=wildcards.sample)


def lima_input_bam(wildcards):
    """Kinnex reads go through skera first; plain Iso-Seq reads go straight in."""
    if KINNEX:
        return f"{RESULTS}/{wildcards.sample}/isoseq/segmented.bam"
    return hifi_bam(wildcards)


# ─── Target ───────────────────────────────────────────────────────────────────

rule all:
    input:
        expand(f"{RESULTS}/{{sample}}/alignment/ig_transcripts.fasta", sample=SAMPLES),
        expand(f"{RESULTS}/{{sample}}/alignment/combined_ig_transcripts.fasta", sample=SAMPLES),
        expand(f"{RESULTS}/{{sample}}/alignment/merge_stats.tsv", sample=SAMPLES),
        expand(f"{RESULTS}/{{sample}}/analysis/exact_match_summary.tsv", sample=SAMPLES),
        expand(f"{RESULTS}/{{sample}}/analysis/top{TOP_N}_alignments.tsv", sample=SAMPLES),
        expand(f"{RESULTS}/{{sample}}/analysis/summary_plot.pdf", sample=SAMPLES),
        f"{RESULTS}/combined/combined_summary_plot.pdf",
        f"{RESULTS}/combined/pooled_summary_plot.pdf",
        f"{RESULTS}/combined/IGH_gene_usage.tsv",
        f"{RESULTS}/combined/IGL_gene_usage.tsv",
        f"{RESULTS}/combined/identity_analysis_plot.pdf",
        f"{RESULTS}/combined/gene_identity_stats.tsv",
        f"{RESULTS}/combined/alignment_position_map.pdf",
        f"{RESULTS}/combined/rank_comparison_plot.pdf",


# ══════════════════════════════════════════════════════════════════════════════
#  BRANCH A — raw BAM input: full IsoSeq preprocessing
# ══════════════════════════════════════════════════════════════════════════════

# Kinnex / MAS-Seq only: split each concatenated HiFi read into its cDNA segments.
rule skera_split:
    input:
        bam=hifi_bam,
        adapters=lambda wc: config["mas_adapters_fasta"],
    output:
        bam=f"{RESULTS}/{{sample}}/isoseq/segmented.bam",
    log:
        f"{RESULTS}/logs/{{sample}}/skera.log",
    threads: 32
    shell:
        "skera split {input.bam} {input.adapters} {output.bam} "
        "-j {threads} 2>{log}"


rule lima:
    input:
        bam=lima_input_bam,
        primers=config["primers_fasta"],
    output:
        bam=f"{RESULTS}/{{sample}}/isoseq/fl.bam",
    params:
        # lima --isoseq always splits its output by primer pair and names the
        # file <prefix>.<5p>--<3p>.bam, so ask for a prefix and rename after.
        prefix=f"{RESULTS}/{{sample}}/isoseq/lima_out",
    log:
        f"{RESULTS}/logs/{{sample}}/lima.log",
    threads: 32
    shell:
        "lima {input.bam} {input.primers} {params.prefix}.bam "
        "--isoseq --peek-guess -j {threads} 2>{log} && "
        "produced=$(ls {params.prefix}.*--*.bam 2>/dev/null | head -1) && "
        "test -n \"$produced\" && "
        "mv \"$produced\" {output.bam} && "
        "mv \"$produced\".pbi {output.bam}.pbi"


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
    threads: 32
    shell:
        "isoseq refine {input.bam} {input.primers} {output.bam} "
        "{params.polya} -j {threads} 2>{log}"


# ─── Optional: screen FLNC reads for IG before clustering ────────────────────
# isoseq cluster2 segfaults during its Consensus stage on very large Kinnex FLNC
# sets (reproduced at 91.3 M reads, isoseq 4.3.0 — the newest release).  Because
# only IG transcripts are ever analysed downstream, screening the FLNC reads
# against the V gene databases first and clustering only the hits sidesteps the
# crash and cuts the clustering input by several orders of magnitude.  Non-IG
# reads would never cluster together with IG reads, so the IG transcripts this
# produces are equivalent to clustering everything.
#
# Caveat: with this on, clustered.fasta already contains only IG reads, so the
# retention_rate in filter_stats.tsv is ~1.0 and no longer measures IG content
# of the library — read it off prefilter_stats.tsv instead.

rule flnc_fasta:
    input:
        f"{RESULTS}/{{sample}}/isoseq/flnc.bam",
    output:
        temp(f"{RESULTS}/{{sample}}/isoseq/flnc.fasta"),
    log:
        f"{RESULTS}/logs/{{sample}}/flnc_fasta.log",
    shell:
        "samtools fasta {input} > {output} 2>{log}"


rule prefilter_align_flnc:
    input:
        reads=f"{RESULTS}/{{sample}}/isoseq/flnc.fasta",
        vgenes=lambda wc: config["vgene_loci"][wc.locus],
    output:
        paf=f"{RESULTS}/{{sample}}/isoseq/prefilter_{{locus}}.paf",
    params:
        preset=config["minimap2_preset"],
    log:
        f"{RESULTS}/logs/{{sample}}/prefilter_{{locus}}.log",
    threads: 64
    shell:
        "minimap2 -cx {params.preset} --cs -t {threads} "
        "{input.vgenes} {input.reads} > {output.paf} 2>{log}"


rule prefilter_ig_flnc:
    input:
        bam=f"{RESULTS}/{{sample}}/isoseq/flnc.bam",
        pafs=expand(f"{RESULTS}/{{sample}}/isoseq/prefilter_{{locus}}.paf",
                    locus=LOCI, allow_missing=True),
    output:
        bam=f"{RESULTS}/{{sample}}/isoseq/ig_flnc.bam",
        ids=f"{RESULTS}/{{sample}}/isoseq/ig_flnc_ids.txt",
        stats=f"{RESULTS}/{{sample}}/isoseq/prefilter_stats.tsv",
    params:
        # Same thresholds as filter_ig_transcripts.py: identity = matches/block,
        # coverage = bases of the V gene spanned.
        min_identity=config["min_alignment_identity"],
        min_coverage_bp=config["min_alignment_coverage_bp"],
    log:
        f"{RESULTS}/logs/{{sample}}/prefilter_ig.log",
    shell:
        "awk -F'\\t' '$11>0 && $10/$11>={params.min_identity} && "
        "($9-$8)>={params.min_coverage_bp} {{print $1}}' {input.pafs} "
        "| sort -u > {output.ids} 2>{log} && "
        "samtools view -N {output.ids} -b -o {output.bam} {input.bam} 2>>{log} && "
        "{{ printf 'metric\\tvalue\\n'; "
        "printf 'flnc_reads\\t%s\\n' \"$(samtools view -c {input.bam})\"; "
        "printf 'ig_flnc_reads\\t%s\\n' \"$(wc -l < {output.ids})\"; "
        "printf 'min_identity_used\\t{params.min_identity}\\n'; "
        "printf 'min_coverage_bp_used\\t{params.min_coverage_bp}\\n'; "
        "}} > {output.stats}"


def cluster_input_bam(wildcards):
    """Cluster the IG-screened FLNC subset when prefiltering is on."""
    stem = "ig_flnc" if PREFILTER_IG else "flnc"
    return f"{RESULTS}/{wildcards.sample}/isoseq/{stem}.bam"


rule make_fofn:
    input:
        cluster_input_bam,
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
    threads: 64
    shell:
        "isoseq {params.cmd} {input.fofn} {output.bam} "
        "{params.singletons} {params.extra} -j {threads} 2>{log}"


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
        script="scripts/filter_ig_transcripts.py",
        transcripts=f"{RESULTS}/{{sample}}/isoseq/clustered.fasta",
        pafs=expand(f"{RESULTS}/{{sample}}/alignment/{{locus}}.paf",
                    locus=LOCI, allow_missing=True),
    output:
        fasta=f"{RESULTS}/{{sample}}/alignment/ig_transcripts.fasta",
        ids=f"{RESULTS}/{{sample}}/alignment/ig_transcript_ids.txt",
        stats=f"{RESULTS}/{{sample}}/alignment/filter_stats.tsv",
    params:
        min_identity=config["min_alignment_identity"],
        min_coverage_bp=config["min_alignment_coverage_bp"],
    log:
        f"{RESULTS}/logs/{{sample}}/filter_ig.log",
    shell:
        "python scripts/filter_ig_transcripts.py "
        "--transcripts {input.transcripts} "
        "--pafs {input.pafs} "
        "--min-identity {params.min_identity} "
        "--min-coverage-bp {params.min_coverage_bp} "
        "--out-fasta {output.fasta} "
        "--out-ids {output.ids} "
        "--out-stats {output.stats} "
        "2>{log}"


# ─── Step 7a: Diversity Analyzer (immunotools) ───────────────────────────────
# Run IGH and IGL separately because tufted_duck has no IGK germline sequences.
# vj_finder requires equal V and J database counts, which breaks with empty IGK files.

rule run_diversity_analyzer_igh:
    input:
        transcripts=f"{RESULTS}/{{sample}}/isoseq/clustered.fasta",
    output:
        fasta=f"{RESULTS}/{{sample}}/immunotools/igh/cleaned_sequences.fasta",
    params:
        outdir=f"{RESULTS}/{{sample}}/immunotools/igh",
        immunotools=config["immunotools_path"],
        org=IMMUNOTOOLS_ORG,
    log:
        f"{RESULTS}/logs/{{sample}}/diversity_analyzer_igh.log",
    threads: 16
    shell:
        "python {params.immunotools} "
        "-i {input.transcripts} "
        "-o {params.outdir} "
        "--org {params.org} "
        "-l IGH "
        "-t {threads} "
        "--skip-plots "
        "2>{log}"


rule run_diversity_analyzer_igl:
    input:
        transcripts=f"{RESULTS}/{{sample}}/isoseq/clustered.fasta",
    output:
        fasta=f"{RESULTS}/{{sample}}/immunotools/igl/cleaned_sequences.fasta",
    params:
        outdir=f"{RESULTS}/{{sample}}/immunotools/igl",
        immunotools=config["immunotools_path"],
        org=IMMUNOTOOLS_ORG,
    log:
        f"{RESULTS}/logs/{{sample}}/diversity_analyzer_igl.log",
    threads: 16
    shell:
        "python {params.immunotools} "
        "-i {input.transcripts} "
        "-o {params.outdir} "
        "--org {params.org} "
        "-l IGL "
        "-t {threads} "
        "--skip-plots "
        "2>{log}"


rule combine_diversity_analyzer:
    input:
        igh=f"{RESULTS}/{{sample}}/immunotools/igh/cleaned_sequences.fasta",
        igl=f"{RESULTS}/{{sample}}/immunotools/igl/cleaned_sequences.fasta",
    output:
        fasta=f"{RESULTS}/{{sample}}/immunotools/cleaned_sequences.fasta",
    shell:
        "cat {input.igh} {input.igl} > {output.fasta}"


# run_immunotools: false — stand in with an empty FASTA so the merge step (and
# everything downstream of it) keeps the same shape and the run reduces to the
# minimap2 filter alone.
rule skip_diversity_analyzer:
    output:
        fasta=f"{RESULTS}/{{sample}}/immunotools/cleaned_sequences.fasta",
    shell:
        "touch {output.fasta}"


if RUN_IMMUNOTOOLS:
    ruleorder: combine_diversity_analyzer > skip_diversity_analyzer
else:
    ruleorder: skip_diversity_analyzer > combine_diversity_analyzer


# ─── Step 7b: Merge minimap2 + immunotools filtered sets ─────────────────────

rule merge_ig_transcripts:
    input:
        script="scripts/merge_ig_fastas.py",
        minimap_fasta=f"{RESULTS}/{{sample}}/alignment/ig_transcripts.fasta",
        immunotools_fasta=f"{RESULTS}/{{sample}}/immunotools/cleaned_sequences.fasta",
    output:
        fasta=f"{RESULTS}/{{sample}}/alignment/combined_ig_transcripts.fasta",
        stats=f"{RESULTS}/{{sample}}/alignment/merge_stats.tsv",
    log:
        f"{RESULTS}/logs/{{sample}}/merge_ig.log",
    shell:
        "python scripts/merge_ig_fastas.py "
        "--minimap-fasta {input.minimap_fasta} "
        "--immunotools-fasta {input.immunotools_fasta} "
        "--out-fasta {output.fasta} "
        "--out-stats {output.stats} "
        "2>{log}"


# ─── Step 7c: Detailed alignment with secondary hits (combined transcripts) ───

rule minimap2_align_detailed:
    input:
        transcripts=f"{RESULTS}/{{sample}}/alignment/combined_ig_transcripts.fasta",
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
        script="scripts/analyze_exact_matches.py",
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
        script="scripts/analyze_top_alignments.py",
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


# ─── Step 10: Per-sample summary plots ────────────────────────────────────────

rule plot_summary:
    input:
        script="scripts/plot_summary.py",
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


# ─── Step 11: Combined cross-sample plot ──────────────────────────────────────

rule plot_combined:
    input:
        script="scripts/plot_combined.py",
        top_alns=expand(f"{RESULTS}/{{sample}}/analysis/top{TOP_N}_alignments.tsv",
                        sample=SAMPLES),
        exact_per_transcript=expand(
            f"{RESULTS}/{{sample}}/analysis/exact_match_per_transcript.tsv",
            sample=SAMPLES),
        exact_summaries=expand(
            f"{RESULTS}/{{sample}}/analysis/exact_match_summary.tsv",
            sample=SAMPLES),
    output:
        plot=f"{RESULTS}/combined/combined_summary_plot.pdf",
        pooled_summary=f"{RESULTS}/combined/pooled_summary_plot.pdf",
        igh_usage=f"{RESULTS}/combined/IGH_gene_usage.tsv",
        igl_usage=f"{RESULTS}/combined/IGL_gene_usage.tsv",
        identity_plot=f"{RESULTS}/combined/identity_analysis_plot.pdf",
        identity_tsv=f"{RESULTS}/combined/gene_identity_stats.tsv",
    params:
        samples=" ".join(SAMPLES),
    log:
        f"{RESULTS}/logs/combined/combined_plot.log",
    shell:
        "python scripts/plot_combined.py "
        "--top-alns {input.top_alns} "
        "--exact-per-transcript {input.exact_per_transcript} "
        "--exact-summaries {input.exact_summaries} "
        "--samples {params.samples} "
        "--output {output.plot} "
        "--out-pooled-summary {output.pooled_summary} "
        "--out-igh-usage {output.igh_usage} "
        "--out-igl-usage {output.igl_usage} "
        "--out-identity-plot {output.identity_plot} "
        "--out-identity-tsv {output.identity_tsv} "
        "2>{log}"


# ─── Step 12: Alignment position maps (match vs mismatch per V gene) ──────────

rule plot_alignment_positions:
    input:
        script="scripts/plot_alignment_positions.py",
        pafs=expand(f"{RESULTS}/{{sample}}/alignment/{{locus}}_detailed.paf",
                    sample=SAMPLES, locus=LOCI),
    output:
        f"{RESULTS}/combined/alignment_position_map.pdf",
    log:
        f"{RESULTS}/logs/combined/alignment_positions.log",
    shell:
        "python scripts/plot_alignment_positions.py "
        "--pafs {input.pafs} "
        "--output {output} "
        "2>{log}"


# ─── Step 13: Rank-1 vs rank-2 alignment comparison ──────────────────────────

rule plot_rank_comparison:
    input:
        script="scripts/plot_rank_comparison.py",
        top_alns=expand(f"{RESULTS}/{{sample}}/analysis/top{TOP_N}_alignments.tsv",
                        sample=SAMPLES),
        pafs=expand(f"{RESULTS}/{{sample}}/alignment/{{locus}}_detailed.paf",
                    sample=SAMPLES, locus=LOCI),
        vgene_fastas=list(config["vgene_loci"].values()),
    output:
        f"{RESULTS}/combined/rank_comparison_plot.pdf",
    params:
        samples=" ".join(SAMPLES),
    log:
        f"{RESULTS}/logs/combined/rank_comparison.log",
    shell:
        "python scripts/plot_rank_comparison.py "
        "--top-alns {input.top_alns} "
        "--pafs {input.pafs} "
        "--vgene-fastas {input.vgene_fastas} "
        "--samples {params.samples} "
        "--output {output} "
        "2>{log}"
