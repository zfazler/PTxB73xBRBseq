#!/bin/bash

# ==============================================================================
# REC -- One-off recovery of 5 samples missing from the pool 1 demultiplex.
#
# Background: PN1_SID5, SID6, SID7, SID8, SID9 (plate positions A05-A09) were
# successfully demultiplexed by BRBseqTools 1.6.1 (3-5M reads each per the
# original stats.txt) but their fastqs never made it into hannah/raw_reads/.
# The raw pool R1/R2 files still exist, so we can re-demultiplex just these
# five barcodes, drop the outputs into raw_reads/ under the correct sample
# names, and let the standard pipeline (03 -> 04 -> 05 -> 06) pick them up.
#
# What this script does:
#   1. Writes a minimal barcodes_recovery.txt with only the 5 missing entries,
#      using the target sample IDs as the "Name" column so BRBseqTools names
#      its output fastqs correctly directly (no rename step needed).
#   2. Runs BRBseqTools Demultiplex on pool 1's R1/R2 using that file.
#   3. Moves the 5 resulting fastqs into hannah/raw_reads/ next to the other
#      379 samples.
#
# Downstream: after this succeeds, re-run q_03, q_04, q_05, q_06 in order.
# Scripts 03-05 skip samples whose outputs already exist, so only the 5 new
# samples get processed. Script 06 regenerates counts across all samples.
# ==============================================================================

set -e

# Java runtime for BRBseqTools jar.
# NCSU has java as a system module — no conda env needed for this step.
module load java/17

# ---- paths ------------------------------------------------------------------
baseDir="/rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez/hannah"
poolDir="/rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez"
rawReads="${baseDir}/raw_reads"
brbTool="${baseDir}/BRBseqTools-1.6.1.jar"

# per-pool R1/R2 (only pool 1 needed for this recovery)
r1="${poolDir}/BZeaBRB1_S1_L004_R1_001.fastq.gz"
r2="${poolDir}/BZeaBRB1_S1_L004_R2_001.fastq.gz"

# temp workspace for this one-off recovery
recDir="${baseDir}/REC_recovery"
recBarcodes="${recDir}/barcodes_recovery.txt"

# ---- validate inputs --------------------------------------------------------
for required in "$brbTool" "$r1" "$r2"; do
    if [ ! -f "$required" ]; then
        echo "ERROR: missing required file: $required"
        exit 1
    fi
done

mkdir -p "$recDir"

# ---- write the minimal barcodes file (5 rows, tab-separated) ---------------
# Barcode sequences come from data/barcodes.txt entries
# for positions A05-A09.
cat > "$recBarcodes" <<'EOF'
Name	B1
PN1_SID5	TTCAATCTCCTTAG
PN1_SID6	CTCGGTTCGAATGC
PN1_SID7	TACACTATAGCTAG
PN1_SID8	CAAGTATAAGGAAC
PN1_SID9	CTGATATGCAGCGA
EOF

echo "=== Recovery barcodes file (5 entries) ==="
cat "$recBarcodes"
echo ""

# ---- refuse to re-run if the recovery outputs already exist -----------------
already=0
for sid in PN1_SID5 PN1_SID6 PN1_SID7 PN1_SID8 PN1_SID9; do
    if [ -f "${rawReads}/${sid}.fastq.gz" ]; then
        already=$((already + 1))
    fi
done
if [ "$already" -eq 5 ]; then
    echo "All 5 recovery fastqs already present in raw_reads/. Nothing to do."
    exit 0
fi

# ---- run BRBseqTools on pool 1 with the 5-barcode file ---------------------
echo "=== Running BRBseqTools Demultiplex on pool 1 with 5 target barcodes ==="
java -jar "$brbTool" Demultiplex \
    -r1 "$r1" \
    -r2 "$r2" \
    -c "$recBarcodes" \
    -p BU??????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????? \
    -UMI 14 \
    -o "$recDir"

# ---- move the 5 recovered fastqs into raw_reads/ ---------------------------
echo ""
echo "=== Moving recovered fastqs to ${rawReads}/ ==="
for sid in PN1_SID5 PN1_SID6 PN1_SID7 PN1_SID8 PN1_SID9; do
    src="${recDir}/${sid}.fastq.gz"
    dst="${rawReads}/${sid}.fastq.gz"
    if [ -f "$src" ]; then
        mv "$src" "$dst"
        echo "  moved: ${sid}.fastq.gz ($(du -h "$dst" | cut -f1))"
    else
        echo "  MISSING after demux: $src"
    fi
done

echo ""
echo "=== Done ==="
echo "Recovery workspace kept at: ${recDir}"
echo "(contains stats.txt and any undetermined output; safe to delete once verified)"
echo ""
echo "Next steps: re-run q_03 -> q_04 -> q_05 -> q_06."
echo "Scripts 03-05 will only process the 5 new samples (existing outputs are skipped)."
