#!/usr/bin/env bash
# Run the conversion detector on one human sample, against its OWN germline.
#
# The negative control for the bird result: humans have no gene conversion, so
# whatever rate this reports is the detector's false-positive floor measured on
# real repertoire data rather than on permuted data. Every step mirrors the bird
# pipeline so the two numbers are comparable.
#
# Deliberately conservative choices, all of which make a spurious call EASIER:
#   * unconstrained parent assignment. No RSS annotation exists for these
#     assemblies, so any V gene may be the parent -- the same setting as the
#     bird "any-parent" run, not the RSS-restricted main run.
#   * a permissive donor pool: every gene is available to every parent. There is
#     no J in the targeted-capture assembly and the locus lands on many contigs,
#     so no deletion/inversion model is possible.
#   * merged Illumina reads (mean expected error ~0.5) rather than HiFi.
#
# usage: gc_human_control.sh SAMPLE
set -euo pipefail
S="${1:?usage: gc_human_control.sh SAMPLE}"
BASE=/local/storage/kav67/human_control
REPO=/home/kav67/isoseq_pipeline
PY=/home/kav67/miniconda3/envs/snakemake/bin/python
VSEARCH=/home/kav67/miniconda3/envs/snakemake/bin/vsearch
MM2=/home/kav67/miniconda3/envs/alignment_env/bin/minimap2
MAN=$REPO/data/human_control/rodriguez_10_paired.tsv
OUT=$BASE/geneconv/$S
mkdir -p "$OUT" "$REPO/data/vgenes/human"

IGD=$BASE/igdetective/$S/combined_genes_IGH.txt
[ -s "$IGD" ] || { echo "$S: no IgDetective output"; exit 1; }

# ── germline: IgDetective V genes, allelic copies collapsed ──────────────────
# The assembly is diploid and fragmented, so the same gene appears once per
# haplotype on different contigs. Left uncollapsed those near-identical pairs
# sit in the donor pool as separate entries and can manufacture tracts -- the
# same near-twin problem that made blackbird IGH unanswerable.
ALL=$OUT/${S}_IGHV_all.fasta
VG=$REPO/data/vgenes/human/${S}_IGHV.fasta
$PY - "$IGD" "$ALL" <<'PYEOF'
import csv, sys
rows=[r for r in csv.DictReader(open(sys.argv[1]),delimiter="\t") if r["GeneType"]=="V"]
tag=sys.argv[2].split("/")[-1].split("_")[0].replace("-","")
with open(sys.argv[2],"w") as fh:
    for r in rows:
        fh.write(f">{tag}_IGH.{r['Pos']}.{r['Contig']}.V.{r['Productive']}.{r['Strand']}\n"
                 f"{r['Sequence'].upper()}\n")
print(f"  {len(rows)} IGHV genes annotated", file=sys.stderr)
PYEOF
$VSEARCH --cluster_fast "$ALL" --id 0.99 --strand both --centroids "$VG" \
         --qmask none --quiet 2>/dev/null
echo "  $(grep -c '>' "$VG") genes after collapsing allelic pairs at 99%"

# ── repertoire: merge pairs, keep V-spanning reads ──────────────────────────
R1=$(awk -F'\t' -v s="$S" '$1==s && $3=="airr_seq"{print $4}' "$MAN")
MERGED=$OUT/${S}_merged.fasta
: > "$MERGED"
for run in $R1; do
  $VSEARCH --fastq_mergepairs "$BASE/raw/${run}_1.fastq" \
           --reverse "$BASE/raw/${run}_2.fastq" \
           --fastaout - --fastq_allowmergestagger \
           --fastq_maxdiffs 15 --fastq_minovlen 20 --quiet 2>/dev/null >> "$MERGED"
done
READS=$OUT/${S}_reads.fasta
$PY - "$MERGED" "$READS" <<'PYEOF'
import random, sys
sys.path.insert(0,"/home/kav67/isoseq_pipeline/scripts")
from gc_lib import read_fasta
random.seed(0)
seqs=read_fasta(sys.argv[1])
keep=[v for v in seqs.values() if len(v)>=450]      # spans the V region
sub=random.sample(keep, min(50000,len(keep)))
with open(sys.argv[2],"w") as fh:
    for i,s in enumerate(sub): fh.write(f">read{i}\n{s}\n")
print(f"  {len(seqs)} merged, {len(keep)} >=450bp, {len(sub)} used", file=sys.stderr)
PYEOF

# ── align, assign, detect ───────────────────────────────────────────────────
PAF=$OUT/${S}_IGH_detailed.paf
$MM2 -cx splice:hq --cs -t 8 "$VG" "$READS" > "$PAF" 2>/dev/null
export PYTHONPATH=$REPO/scripts
$PY $REPO/scripts/gc_call_functional_genes.py \
  --pafs "$PAF" --vgene-fastas "$VG" --junction-margin 20 \
  --min-identity 0.90 --min-covered-bp 200 --min-transcripts 1 \
  --out-genes "$OUT/functional_genes.tsv" \
  --out-assignments "$OUT/transcript_assignments.tsv" 2>/dev/null

$PY - "$VG" "$OUT/IGH_donor_pool.tsv" <<'PYEOF'
import sys
sys.path.insert(0,"/home/kav67/isoseq_pipeline/scripts")
from gc_lib import read_fasta
names=sorted(read_fasta(sys.argv[1]))
with open(sys.argv[2],"w") as fh:
    fh.write("rearranged_gene\tlocus\tpos\tstrand\tmechanism\tn_allowed_donors\t"
             "n_deleted\tallowed_donors\tdeleted_genes\n")
    for g in names:
        o=[d for d in names if d!=g]
        fh.write(f"{g}\tIGH\t{g.split('.')[1]}\t{g.split('.')[-1]}\tunknown\t"
                 f"{len(o)}\t0\t{','.join(o)}\tNONE\n")
PYEOF

$PY $REPO/scripts/gc_detect_tracts.py \
  --paf "$PAF" --vgene-fasta "$VG" \
  --assignments "$OUT/transcript_assignments.tsv" \
  --donor-pool "$OUT/IGH_donor_pool.tsv" --locus IGH \
  --min-informative 3 --max-gap-bp 5 --p-threshold 0.05 --correction support \
  --min-support-significant 6 \
  --out-tracts "$OUT/IGH_tracts.tsv" \
  --out-summary "$OUT/IGH_tract_summary.tsv" 2>&1 | sed 's/^/  /'
echo "$S COMPLETE"
