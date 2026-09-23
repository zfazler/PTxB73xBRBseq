#!/bin/bash

# activate conda env (minimap2 + samtools + pigz need to be on PATH)
module load conda
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /usr/local/usrapps/maize/hdpil/hdpil

# Set paths
#MINIMAP2_PATH="/programs/minimap2-2.27/minimap2"
REFERENCE_PATH="/rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez/hannah/Arabidopsis_Zea_rRNA.fasta"
TRIMMED_DIR="/rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez/hannah/03_trimmed"
OUTPUT_DIR="/rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez/hannah/clean_reads"
LOG_FILE="/rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez/hannah/rRNA_survival_statistics.log"

# Create output and log directories if they don't exist
mkdir -p "$OUTPUT_DIR"
if [ ! -f "$LOG_FILE" ]; then
    echo "Sample ID, Total Reads, Aligned Reads, Unaligned Reads, % Aligned, % Unaligned" > "$LOG_FILE"
fi

# Loop through all trimmed FASTQ files
for fq_file in "$TRIMMED_DIR"/*.fq.gz; do
    # Extract sample ID
    sample_id=$(basename "$fq_file" "_trimmed.fq.gz")

    # Define output clean FASTQ filename
    clean_fq_file="${OUTPUT_DIR}/${sample_id}_clean.fastq.gz"

    # Check if the clean FASTQ file already exists
    if [ -f "$clean_fq_file" ]; then
        echo "Clean reads for ${sample_id} already exist, skipping..."
        continue
    fi

    # Define intermediate SAM and BAM filenames
    sam_file="${sample_id}_rRNA.sam"
    bam_file="${sample_id}_rRNA.bam"

    # Align reads with minimap2 and convert SAM to BAM
    echo "Aligning reads for ${sample_id}..."
    minimap2 -ax sr -t 60 "$REFERENCE_PATH" "$fq_file" > "$sam_file"
    samtools view -bS "$sam_file" > "$bam_file"

    # Count total, aligned, and unaligned reads
    total_reads=$(samtools view -c "$bam_file")
    aligned_reads=$(samtools view -c -F 4 "$bam_file")
    unaligned_reads=$(samtools view -c -f 4 "$bam_file")

    # Calculate percentages
    percentage_aligned=$(echo "scale=2; $aligned_reads / $total_reads * 100" | bc)
    percentage_unaligned=$(echo "scale=2; $unaligned_reads / $total_reads * 100" | bc)

    # Output clean reads to FASTQ
    echo "Extracting unaligned (clean) reads for ${sample_id}..."
    samtools view -b -f 4 "$bam_file" | samtools fastq - | pigz > "$clean_fq_file"

    # Delete the intermediate SAM and BAM files to clean up
    rm "$sam_file"
    rm "$bam_file"

    # Log the statistics
    echo "${sample_id}, $total_reads, $aligned_reads, $unaligned_reads, $percentage_aligned%, $percentage_unaligned%" >> "$LOG_FILE"

    # Status message
    echo "Finished processing $sample_id. Clean reads saved to $clean_fq_file"
done

echo "All files processed. Statistics saved to $LOG_FILE."
