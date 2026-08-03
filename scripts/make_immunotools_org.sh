#!/bin/bash
# Build an immunotools organism database (germline + CDR annotation) from a set
# of V gene FASTAs, and register it so `diversity_analyzer.py --org <ORG>` works.
#
#   usage: make_immunotools_org.sh <ORG> <IGHV.fasta> <IGLV.fasta>
#   e.g.   make_immunotools_org.sh VGP_redwinged_blackbird \
#              data/vgenes/VGP_redwinged_blackbird_IGH.fasta \
#              data/vgenes/VGP_redwinged_blackbird_IGL.fasta
#
# What it produces
#   immunotools/data/germline/<ORG>/IG/{IGHV,IGLV,IGHJ,IGLJ,IGKV,IGKJ}.fa
#   immunotools/data/annotation/<ORG>_v_imgt.txt
#   immunotools/data/annotation/<ORG>_j_imgt.txt
#
# J genes: the bird annotations we work from contain V genes only, so the
# Gallus IMGT J genes are borrowed from the chicken_imgt database — the same
# stand-in the existing tufted_duck database uses.  Note that vj_finder's
# min_j_segment_length is 30 and songbird IGL transcripts align only ~28 bp to
# the Gallus IGLJ, so IGL recall through immunotools is poor until real J genes
# are available.  Replace IGLJ.fa/IGHJ.fa here as soon as they are.
#
# V CDR annotation: igblastn transfers IMGT FR/CDR boundaries from the closest
# germline hit, so it needs a reference organism that already has an
# .ndm.imgt.  igblast ships none for birds, so this script builds one on the
# fly from immunotools' own chicken_imgt data — its *_v_imgt.txt is already in
# exactly the .ndm.imgt column layout.  Genes whose FR/CDR set comes back
# incomplete are dropped by generate_v_cdr_labeling_by_igblast.py; that is
# expected (tufted_duck kept 76/121).
set -euo pipefail

ORG=${1:?usage: make_immunotools_org.sh <ORG> <IGHV.fasta> <IGLV.fasta> [IGHJ.fasta] [IGLJ.fasta]}
IGHV_IN=${2:?missing IGHV fasta}
IGLV_IN=${3:?missing IGLV fasta}
# Optional species-specific J genes.  Anything not supplied falls back to the
# Gallus IMGT J gene from chicken_imgt (see the note at the top of this file).
IGHJ_IN=${4:-}
IGLJ_IN=${5:-}

IMMUNOTOOLS=${IMMUNOTOOLS:-/home/kav67/immunotools}
IGBLAST_BIN=${IGBLAST_BIN:-/home/kav67/miniconda3/envs/ig-assembly-eval/bin}
MAKEBLASTDB=${MAKEBLASTDB:-/programs/bin/blast+/makeblastdb}

GERM=$IMMUNOTOOLS/data/germline/$ORG/IG
ANNO=$IMMUNOTOOLS/data/annotation
REF=$IMMUNOTOOLS/data/germline/chicken_imgt/IG
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

# ── 1. germline ──────────────────────────────────────────────────────────────
mkdir -p "$GERM"
cp "$IGHV_IN" "$GERM/IGHV.fa"
cp "$IGLV_IN" "$GERM/IGLV.fa"
cp "${IGHJ_IN:-$REF/IGHJ.fa}" "$GERM/IGHJ.fa"
cp "${IGLJ_IN:-$REF/IGLJ.fa}" "$GERM/IGLJ.fa"
: > "$GERM/IGKV.fa"
: > "$GERM/IGKJ.fa"
for f in "$GERM/IGHJ.fa" "$GERM/IGLJ.fa"; do
    [ -n "$(tail -c1 "$f")" ] && echo >> "$f"
done
echo "germline written to $GERM"

# ── 2. igblast reference built from chicken_imgt ─────────────────────────────
mkdir -p "$WORK/internal_data/chicken" "$WORK/optional_file" "$WORK/db"
cat "$REF/IGHV.fa" "$REF/IGLV.fa" > "$WORK/db/chicken_V.fa"
cat "$GERM/IGHJ.fa" "$GERM/IGLJ.fa" > "$WORK/db/chicken_J.fa"
{
    echo "#gene/allele name, FWR1 start, FWR1 stop, CDR1 start, CDR1 stop, FWR2 start, FWR2 stop, CDR2 start, CDR2 stop, FWR3 start, FWR3 stop, chain type, coding frame start. "
    echo "#FWR/CDR positions are 1-based while the coding frame start positions are 0-based "
    cat "$ANNO/chicken_imgt_v_imgt.txt"
} > "$WORK/internal_data/chicken/chicken.ndm.imgt"
printf '#gene/allele name, first coding frame start position, chain type, CDR3 stop, extra bps beyond J coding end.\n#All positions are 0-based\n\nIMGT000014|IGHJ*01|Gallus\t2\tJH\t22\t0\nIMGT000008|IGLJ*01|Gallus\t1\tJL\t6\t0\n' \
    > "$WORK/optional_file/chicken_gl.aux"
"$MAKEBLASTDB" -in "$WORK/db/chicken_V.fa" -dbtype nucl -parse_seqids -out "$WORK/db/chicken_V" > /dev/null
"$MAKEBLASTDB" -in "$WORK/db/chicken_J.fa" -dbtype nucl -parse_seqids -out "$WORK/db/chicken_J" > /dev/null
"$MAKEBLASTDB" -in "$WORK/db/chicken_V.fa" -dbtype nucl -parse_seqids -out "$WORK/internal_data/chicken/chicken_V" > /dev/null

# ── 3. V CDR labelling via igblast ───────────────────────────────────────────
export IGDATA=$WORK
for locus in IGHV IGLV; do
    "$IGBLAST_BIN/igblastn" \
        -germline_db_V "$WORK/db/chicken_V" \
        -germline_db_J "$WORK/db/chicken_J" \
        -germline_db_D "$WORK/db/chicken_J" \
        -auxiliary_data "$WORK/optional_file/chicken_gl.aux" \
        -organism chicken -domain_system imgt -ig_seqtype Ig \
        -query "$GERM/$locus.fa" -outfmt 7 -num_threads 8 \
        -out "$WORK/$locus.igblast.txt"
    python "$IMMUNOTOOLS/py/cdr_labeling_utils/generate_v_cdr_labeling_by_igblast.py" \
        "$WORK/$locus.igblast.txt" "$WORK/$locus.labels" > /dev/null
done
cat "$WORK/IGHV.labels" "$WORK/IGLV.labels" > "$ANNO/${ORG}_v_imgt.txt"
echo "V annotation: $(wc -l < "$ANNO/${ORG}_v_imgt.txt") of $(grep -hc '>' "$GERM/IGHV.fa" "$GERM/IGLV.fa" | paste -sd+ | bc) V genes labelled"

# ── 4. J CDR3 labelling ──────────────────────────────────────────────────────
python "$IMMUNOTOOLS/py/cdr_labeling_utils/create_j_gene_cdr3_labeling.py" \
    "$GERM/IGHJ.fa" "$WORK/jh.txt" IGH
# The script's IGL motif list (TTAGG/TTCGG/TTCAT) does not match the Gallus
# IGLJ, whose conserved FGXG motif is TTTGGGGC at position 7 — the position the
# existing tufted_duck_j_imgt.txt already records for this same sequence.
python "$IMMUNOTOOLS/py/cdr_labeling_utils/create_j_gene_cdr3_labeling.py" \
    "$GERM/IGLJ.fa" "$WORK/jl.txt" IGL || true
if [ ! -s "$WORK/jl.txt" ]; then
    grep 'LJ$' "$ANNO/tufted_duck_j_imgt.txt" > "$WORK/jl.txt"
fi
cat "$WORK/jh.txt" "$WORK/jl.txt" > "$ANNO/${ORG}_j_imgt.txt"
echo "J annotation: $(wc -l < "$ANNO/${ORG}_j_imgt.txt") J genes labelled"

# ── 5. register the organism ─────────────────────────────────────────────────
if grep -q "'$ORG'" "$IMMUNOTOOLS/diversity_analyzer.py"; then
    echo "$ORG already registered in diversity_analyzer.py"
else
    echo "NOTE: add \"'$ORG' : '$ORG',\" to organism_dict in $IMMUNOTOOLS/diversity_analyzer.py"
fi
echo "done: --org $ORG"
