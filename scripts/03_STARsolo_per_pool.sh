#!/bin/bash

# ==============================================================================
# 03 -- STARsolo per pool (HPC)
#
# Runs one STARsolo alignment per pool. Replaces the legacy chain of:
#     02_demultiplex (BRBseqTools)  -> 03_trimming_and_QC  ->
#     04_rRNA_filtering             -> 05_STAR_alignment   ->
#     06_featureCounts              -> 06b_UMI_dedup (would-be) ->
#     06c_featureCounts_UMI (would-be)
#
# One STARsolo invocation handles: sample barcode demux (from R1 first 14 nt) +
# UMI extraction (R1 15-28) + adapter clipping (CellRanger4) + alignment +
# UMI-collapsed gene counting. Dual-dedup mode ("1MM_Directional NoDedup")
# emits BOTH the deduplicated and raw-read matrices in one alignment run so
# we can compare them without re-aligning.
#
# Matches Alithea's July 2026 data-analysis-manual §1.4 command exactly except
# for the dual-dedup flag.
#
# Usage:
#     bash scripts/03_STARsolo_per_pool.sh <POOL_N>
# where POOL_N is 1..4.
#
# Inputs (HPC):
#     $baseDir/trimmed/pool_{N}_R1.fastq.gz              trimmed R1 (from 02b)
#     $baseDir/trimmed/pool_{N}_R2.fastq.gz              trimmed R2 (from 02b)
#     $repoDir/data/starsolo/barcode_whitelist.txt       from 02_prepare_barcodes.R
#     $baseDir/Zea_mays/genomeIndex/                     existing STAR index
#     $baseDir/Zea_mays/Zea_mays.gtf                     annotation
#
# Note: trimming is a separate stage (02b_trim_pools.sh) so we can control
# what happens with the 150 PE overshoot. STARsolo's --clipAdapterType
# CellRanger4 is tuned for the ~90 nt R2 length Alithea's manual assumes,
# and would soft-clip inefficiently on 150 nt R2 with unclipped adapter.
#
# Outputs (HPC, per pool):
#     $baseDir/starsolo/pool_N/Aligned.sortedByCoord.out.bam
#     $baseDir/starsolo/pool_N/Log.final.out
#     $baseDir/starsolo/pool_N/Solo.out/Gene/raw/
#         features.tsv                       gene IDs
#         barcodes.tsv                       barcodes that matched the whitelist
#         umiDedup-1MM_Directional.mtx       UMI-collapsed counts (canonical)
#         umiDedup-NoDedup.mtx               raw read counts (comparison)
#     $baseDir/starsolo/pool_N/Solo.out/Barcodes.stats   barcode QC
#
# Cost: ~2-6 h per pool depending on pool size. Pool 1 is the largest at 40 GB R1.
# ==============================================================================

set -euo pipefail

pool="${1:-}"
if [[ ! "$pool" =~ ^[1-4]$ ]]; then
  echo "usage: $0 <POOL_N>   (POOL_N = 1, 2, 3, or 4)" >&2
  exit 1
fi

# --- paths ----------------------------------------------------------------
baseDir="/rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez/hannah"
repoDir="$baseDir/BZeaBRBseq"
trimDir="${baseDir}/trimmed"

R1="${trimDir}/pool_${pool}_R1.fastq.gz"
R2="${trimDir}/pool_${pool}_R2.fastq.gz"
whitelist="${repoDir}/data/starsolo/barcode_whitelist.txt"
starIndex="/rsstu/users/r/rrellan/sara/ref/STAR_index"    # shared lab STAR index (same one the legacy 05_STAR_alignment.sh used)

outDir="${baseDir}/starsolo/pool_${pool}"
mkdir -p "$outDir"

# --- environment ----------------------------------------------------------
module load conda
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /usr/local/usrapps/maize/hdpil/hdpil

# --- validate inputs ------------------------------------------------------
for f in "$R1" "$R2" "$whitelist"; do
  [ -f "$f" ] || { echo "ERROR: missing input: $f" >&2; exit 1; }
done
[ -d "$starIndex" ] || { echo "ERROR: STAR index dir missing: $starIndex" >&2; exit 1; }
command -v STAR >/dev/null || { echo "ERROR: STAR not on PATH" >&2; exit 1; }
STAR --version

echo "=== STARsolo pool $pool ==="
echo "R1:        $R1  ($(du -h "$R1" | cut -f1))"
echo "R2:        $R2  ($(du -h "$R2" | cut -f1))"
echo "Whitelist: $whitelist  ($(wc -l < "$whitelist") barcodes)"
echo "Index:     $starIndex"
echo "Output:    $outDir/"
echo ""

# --- STARsolo -------------------------------------------------------------
# Note: readFilesIn takes GENOMIC read (R2) first, then BARCODE read (R1).
cd "$outDir"
STAR \
  --runMode alignReads \
  --runThreadN 16 \
  --genomeDir "$starIndex" \
  --readFilesIn "$R2" "$R1" \
  --readFilesCommand zcat \
  --soloType CB_UMI_Simple \
  --soloCBwhitelist "$whitelist" \
  --soloCBstart 1 --soloCBlen 14 \
  --soloUMIstart 15 --soloUMIlen 14 \
  --soloBarcodeReadLength 0 \
  --soloCBmatchWLtype 1MM \
  --soloUMIdedup 1MM_Directional NoDedup \
  --soloStrand Forward \
  --soloCellFilter None \
  --clipAdapterType CellRanger4 \
  --outFilterMultimapNmax 1 \
  --outFilterIntronMotifs None \
  --outSAMmapqUnique 60 \
  --outSAMunmapped Within \
  --outSAMmultNmax 1 \
  --outSAMtype BAM SortedByCoordinate \
  --outBAMsortingThreadN 1 \
  --outBAMsortingBinsN 4 \
  --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
  --limitBAMsortRAM 40000000000 \
  --outFileNamePrefix "${outDir}/"

echo ""
echo "=== pool $pool done ==="
echo "BAM:            ${outDir}/Aligned.sortedByCoord.out.bam"
echo "UMI matrix:     ${outDir}/Solo.out/Gene/raw/umiDedup-1MM_Directional.mtx"
echo "Raw matrix:     ${outDir}/Solo.out/Gene/raw/umiDedup-NoDedup.mtx"
echo "Barcodes:       ${outDir}/Solo.out/Gene/raw/barcodes.tsv"
echo "Features:       ${outDir}/Solo.out/Gene/raw/features.tsv"
echo "Alignment log:  ${outDir}/Log.final.out"
