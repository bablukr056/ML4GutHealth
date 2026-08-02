#!/bin/bash
#
# Metagenomic read preprocessing and taxonomic classification pipeline
#
# Steps:
#   1. Quality trimming and adapter removal (fastp)
#   2. Host (human) read removal (Bowtie2)
#   3. Taxonomic classification (Kraken2)
#   4. Abundance re-estimation (Bracken)
#
# Usage:
#   ./bioinfor_pipeline.sh <BIOPROJECT> <SAMPLE_FOLDER> <SAMPLE>
#
# Arguments:
#   BIOPROJECT     BioProject accession (used for input/output directory naming)
#   SAMPLE_FOLDER  Subfolder under the BioProject containing the raw FASTQ files
#   SAMPLE         Sample identifier (Run accession); expects <SAMPLE>_pass_1.fastq.gz
#                  and, for paired-end data, <SAMPLE>_pass_2.fastq.gz
#
# Execution environment:
#   This script was executed as an HTCondor job on the HPC cluster
#   (submitted via `condor_submit`). Per-sample arguments (BIOPROJECT,
#   SAMPLE_FOLDER, SAMPLE) were supplied through a job submission text
#   file listing one set of arguments per sample/job, enabling batch
#   submission of the pipeline across all samples as independent Condor
#   jobs.
# CPU and memory resources were allocated per job in the Condor submission file.
# 4 CPU cores and 128 GB of memory were requested for each job (or each samples), which is sufficient for the computationally intensive steps of the pipeline (Bowtie2, Kraken2, Bracken).
#
set -uo pipefail

# ------------------------------------------------------------------
# Environment
# ------------------------------------------------------------------
source /lustrehome/babluuniba2022/miniconda3/bin/activate
conda activate my_pipeline # Activate the conda environment containing fastp, bowtie2, kraken2, and bracken

# ------------------------------------------------------------------
# Arguments and paths
# ------------------------------------------------------------------
if [ $# -lt 3 ]; then
  echo "Usage: $0 <BIOPROJECT> <SAMPLE_FOLDER> <SAMPLE>"
  exit 1
fi

BIOPROJECT="$1"
SAMPLE_FOLDER="$2"
SAMPLE="$3"

INPUT_FOLDER="/lustre/SEQUENCING/kumar_phd_prj/${BIOPROJECT}/${SAMPLE_FOLDER}"
OUTPUT_FOLDER="/lustre/home/babluuniba2022/DEC"
KRAKEN2_DB="/lustrehome/babluuniba2022/database/kraken2_mgnify_genome"
BOWTIE2_DB="/lustrehome/babluuniba2022/database/bowtie2_GRCh38/GRCh38_noalt_as/GRCh38_noalt_as"
ERROR_LOG="error_log.txt"
THREADS=128

if [ ! -d "$INPUT_FOLDER" ]; then
  echo "[$(date '+%F %T')] ERROR: input folder does not exist: $INPUT_FOLDER" >> "$ERROR_LOG"
  exit 1
fi

mkdir -p \
  "$OUTPUT_FOLDER/fastp/$BIOPROJECT" \
  "$OUTPUT_FOLDER/bowtie2/$BIOPROJECT" \
  "$OUTPUT_FOLDER/kraken2/$BIOPROJECT" \
  "$OUTPUT_FOLDER/bracken/$BIOPROJECT"

R1="$INPUT_FOLDER/${SAMPLE}_pass_1.fastq.gz"
R2="$INPUT_FOLDER/${SAMPLE}_pass_2.fastq.gz"

IS_PAIRED=true
if [ ! -f "$R2" ]; then
  IS_PAIRED=false
fi

log_error_and_exit() {
  echo "[$(date '+%F %T')] ERROR: $1 (sample: $SAMPLE, bioproject: $BIOPROJECT)" >> "$ERROR_LOG"
  exit 1
}

# ------------------------------------------------------------------
# Step 1: Quality trimming and adapter removal (fastp)
# ------------------------------------------------------------------
echo "[Step 1] Running fastp for sample: $SAMPLE"

FASTP_R1="$OUTPUT_FOLDER/fastp/$BIOPROJECT/${SAMPLE}_trimmed_R1.fastq.gz"
FASTP_R2="$OUTPUT_FOLDER/fastp/$BIOPROJECT/${SAMPLE}_trimmed_R2.fastq.gz"
FASTP_JSON="$OUTPUT_FOLDER/fastp/$BIOPROJECT/${SAMPLE}_fastp.json"
FASTP_HTML="$OUTPUT_FOLDER/fastp/$BIOPROJECT/${SAMPLE}_fastp.html"

if $IS_PAIRED; then
  fastp \
    -i "$R1" -I "$R2" \
    -o "$FASTP_R1" -O "$FASTP_R2" \
    --detect_adapter_for_pe --dont_overwrite \
    -q 20 -u 20 -l 45 \
    --thread "$THREADS" \
    -j "$FASTP_JSON" -h "$FASTP_HTML"
else
  fastp \
    -i "$R1" \
    -o "$FASTP_R1" \
    --dont_overwrite \
    -q 20 -u 20 -l 45 \
    --thread "$THREADS" \
    -j "$FASTP_JSON" -h "$FASTP_HTML"
fi

[ $? -eq 0 ] || log_error_and_exit "fastp failed"

# ------------------------------------------------------------------
# Step 2: Host read removal (Bowtie2 against GRCh38)
# ------------------------------------------------------------------
echo "[Step 2] Running Bowtie2 host filtering for sample: $SAMPLE"

NONHUMAN_PREFIX="$OUTPUT_FOLDER/bowtie2/$BIOPROJECT/${SAMPLE}_nonhuman.fq.gz"

if $IS_PAIRED; then
  bowtie2 \
    -x "$BOWTIE2_DB" \
    -1 "$FASTP_R1" -2 "$FASTP_R2" \
    -p "$THREADS" \
    --un-conc-gz "$NONHUMAN_PREFIX" \
    -S /dev/null
  NONHUMAN_R1="${NONHUMAN_PREFIX}.1"
  NONHUMAN_R2="${NONHUMAN_PREFIX}.2"
else
  bowtie2 \
    -x "$BOWTIE2_DB" \
    -U "$FASTP_R1" \
    -p "$THREADS" \
    --un-gz "$NONHUMAN_PREFIX" \
    -S /dev/null
  NONHUMAN_SE="$NONHUMAN_PREFIX"
fi

[ $? -eq 0 ] || log_error_and_exit "Bowtie2 failed"

# ------------------------------------------------------------------
# Step 3: Taxonomic classification (Kraken2)
# ------------------------------------------------------------------
echo "[Step 3] Running Kraken2 for sample: $SAMPLE"

KRAKEN_OUT="$OUTPUT_FOLDER/kraken2/$BIOPROJECT/${SAMPLE}.kraken2.out"
KRAKEN_REPORT="$OUTPUT_FOLDER/kraken2/$BIOPROJECT/${SAMPLE}.k2report"
KRAKEN_CLASSIFIED="$OUTPUT_FOLDER/kraken2/$BIOPROJECT/${SAMPLE}_classified.fastq"
KRAKEN_UNCLASSIFIED="$OUTPUT_FOLDER/kraken2/$BIOPROJECT/${SAMPLE}_unclassified.fastq"

if $IS_PAIRED; then
  kraken2 \
    --db "$KRAKEN2_DB" \
    --output "$KRAKEN_OUT" \
    --report "$KRAKEN_REPORT" \
    --classified-out "$KRAKEN_CLASSIFIED" \
    --unclassified-out "$KRAKEN_UNCLASSIFIED" \
    --threads "$THREADS" \
    --paired "$NONHUMAN_R1" "$NONHUMAN_R2"
else
  kraken2 \
    --db "$KRAKEN2_DB" \
    --output "$KRAKEN_OUT" \
    --report "$KRAKEN_REPORT" \
    --classified-out "$KRAKEN_CLASSIFIED" \
    --unclassified-out "$KRAKEN_UNCLASSIFIED" \
    --threads "$THREADS" \
    "$NONHUMAN_SE"
fi

[ $? -eq 0 ] || log_error_and_exit "Kraken2 failed"

# ------------------------------------------------------------------
# Step 4: Species-level abundance re-estimation (Bracken)
# ------------------------------------------------------------------
echo "[Step 4] Running Bracken for sample: $SAMPLE"

BRACKEN_OUT="$OUTPUT_FOLDER/bracken/$BIOPROJECT/${SAMPLE}.bracken.out"
READ_LEN=100
TAXONOMIC_LEVEL="S"
READ_THRESHOLD=10

bracken \
  -d "$KRAKEN2_DB" \
  -i "$KRAKEN_REPORT" \
  -o "$BRACKEN_OUT" \
  -r "$READ_LEN" \
  -l "$TAXONOMIC_LEVEL" \
  -t "$READ_THRESHOLD"

[ $? -eq 0 ] || log_error_and_exit "Bracken failed"

echo "[Complete] Pipeline finished successfully for sample: $SAMPLE"

# ------------------------------------------------------------------
