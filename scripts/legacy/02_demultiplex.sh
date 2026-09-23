#!/bin/bash

# ==============================================================================
# 02 -- BRB-seq demultiplexing
# Split each multiplexed pool's R2 fastq (transcriptome reads) into per-sample
# fastqs using the barcode information from R1.
#
# Tool: BRBseqTools 1.6.1 (https://github.com/DeplanckeLab/BRB-seqTools)
#
# Inputs (at parent-level dir):
#   BZeaBRB[1-4]_S[1-4]_L004_R1_001.fastq.gz   (barcode + UMI reads)
#   BZeaBRB[1-4]_S[1-4]_L004_R2_001.fastq.gz   (transcriptome reads)
# Reference files (in hannah/):
#   BRBseqTools-1.6.1.jar
#   barcodes.txt   (tab-separated: sample_id \t barcode_sequence)
# Output:
#   hannah/raw_reads/{POOL}/    per-sample R2 fastqs, ready for step 03
# ==============================================================================

set -e

# activate conda env (Java + BRBseqTools need to be on PATH)
module load conda
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /usr/local/usrapps/maize/hdpil/hdpil

# ---- paths ------------------------------------------------------------------
baseDir="/rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez/hannah"
poolDir="/rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez"
outDIR="${baseDir}/raw_reads/"
barcodes="${baseDir}/barcodes.txt"
brbTool="${baseDir}/BRBseqTools-1.6.1.jar"

# ---- validate inputs --------------------------------------------------------
for required in "$brbTool" "$barcodes"; do
    if [ ! -f "$required" ]; then
        echo "ERROR: missing required file: $required"
        echo ""
        echo "  BRBseqTools jar: download from"
        echo "    https://github.com/DeplanckeLab/BRB-seqTools/releases"
        echo "  barcodes.txt: tab-separated, sample_id <tab> barcode_sequence"
        exit 1
    fi
done

mkdir -p "$outDIR"

# ---- pools ------------------------------------------------------------------
pools=("BZeaBRB1_S1_L004" "BZeaBRB2_S2_L004" "BZeaBRB3_S3_L004" "BZeaBRB4_S4_L004")

for pool in "${pools[@]}"; do
    r1="${poolDir}/${pool}_R1_001.fastq.gz"
    r2="${poolDir}/${pool}_R2_001.fastq.gz"

    if [ ! -f "$r1" ] || [ ! -f "$r2" ]; then
        echo "Skipping ${pool}: raw pool fastqs not found at ${poolDir}"
        continue
    fi

    poolOut="${outDIR}${pool}"
    if [ -d "$poolOut" ] && [ -n "$(ls -A "$poolOut" 2>/dev/null)" ]; then
        echo "Demultiplex output already exists for ${pool}, skipping..."
        continue
    fi
    mkdir -p "$poolOut"

    echo "=== Demultiplexing ${pool} ==="
    java -jar "$brbTool" Demultiplex \
        -r1 "$r1" \
        -r2 "$r2" \
        -c "$barcodes" \
        -p BU??????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????? \
        -UMI 14 \
        -o "$poolOut"
done

echo "All demultiplexing complete. Per-sample fastqs written under: ${outDIR}"
