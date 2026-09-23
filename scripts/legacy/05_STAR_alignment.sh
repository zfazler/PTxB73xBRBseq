#!/bin/bash

# activate conda environment
module load conda
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /usr/local/usrapps/maize/hdpil/hdpil

# Check if a species name is provided

if [ -z "$1" ]; then
    echo "Usage: $0 <Species Name>"
    exit 1
fi

species=$1

# Set directories based on the species name
baseDir="/rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez/hannah"
metadata="${baseDir}/metadata.txt"
inDIR="${baseDir}/clean_reads/"
outDIR="${baseDir}/alignments/"
genomeDir="/rsstu/users/r/rrellan/sara/ref/STAR_index/"
sampleIDs="${baseDir}/sample_ids.txt"

# Ensure output directories exist
mkdir -p "$outDIR"

species="Zea_mays"

# Extract relevant sample IDs based on species from metadata
awk -v sp="$species" -F'\t' '$2 == sp {print $1}' "$metadata" > "$sampleIDs"

# Perform alignments using the extracted sample IDs
while read sample_id; do
    fastqFile="${inDIR}${sample_id}_clean.fastq.gz"

    if [ -f "$fastqFile" ]; then
        if [ ! -f "${outDIR}${sample_id}_Aligned.sortedByCoord.out.bam" ]; then
            echo "Running alignment for $sample_id..."
            STAR \
                --outFileNamePrefix "${outDIR}${sample_id}_" \
                --outSAMtype BAM SortedByCoordinate \
                --outSAMstrandField intronMotif \
                --genomeDir "${genomeDir}" \
                --runThreadN 12 \
                --readFilesCommand zcat \
                --readFilesIn "$fastqFile" \
                --twopassMode Basic
        else
            echo "Alignment already completed for $sample_id, skipping..."
        fi
    else
        echo "FASTQ file not found for $sample_id"
    fi
done < "$sampleIDs"

echo "All alignments for $species completed."

# Call the appropriate R script for counting reads based on species
#if [ "$species" = "Zea_mays" ]; then
#    Rscript /rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez/hannah/06_featureCounts_Zm.R "$species"
#else
#    Rscript /rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez/hannah/06_featureCounts.R "$species"
#fi

# After counting, run the summary statistics script
#bash /rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez/hannah/07_generate_summary_statistics.sh "$species"
