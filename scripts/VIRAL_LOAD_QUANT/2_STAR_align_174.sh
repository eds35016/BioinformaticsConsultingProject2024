#!/bin/bash

# Define variables
RUN_NAME="PRRSV2_174_QUANT"
FASTA_DIR_NAME="PRRSV2"
RUN_ID="CKDL210021614-1a-AK7207-AK7213_HMH2WDSX2"
GENOME_DIR="/scratch/bioconsult/Ethan_James/BINF_Consulting_2024/ongoing/mapping/STAR_idx/${FASTA_DIR_NAME}/pig_PRRSV2_174"
IN_DIR="/scratch/bioconsult/Ethan_James/BINF_Consulting_2024/raw_data/${FASTA_DIR_NAME}"
OUT_DIR="/scratch/bioconsult/Ethan_James/BINF_Consulting_2024/ongoing/mapping/${RUN_NAME}"
METADATA_FILE="/scratch/bioconsult/Ethan_James/BINF_Consulting_2024/metadata/PRRSV2_metadata.csv"

# Create output directory if needed
mkdir -p $OUT_DIR

# Extract sample IDs
SAMPLE_IDS=($(awk -F ',' 'NR>1 {print $2}' ${METADATA_FILE}))

# Check if there are any samples to process
if [ ${#SAMPLE_IDS[@]} -eq 0 ]; then
    echo "No samples to process. Exiting."
    exit 1
fi

# Print sample IDs to terminal
echo "Processing the following samples:"
for id in "${SAMPLE_IDS[@]}"; do
  echo $id
done

# Exporting necessary variables for the Slurm script
export RUN_NAME
export RUN_ID
export GENOME_DIR
export IN_DIR
export OUT_DIR
export SAMPLE_IDS
export METADATA_FILE

# Create a string with comma-separated sample IDs for exporting
export SAMPLE_IDS_STR="${SAMPLE_IDS[*]}"

# Submitting the slurm array job. 
sbatch --array=1-${#SAMPLE_IDS[@]}%12 --export=RUN_NAME=${RUN_NAME},RUN_ID=${RUN_ID},GENOME_DIR=${GENOME_DIR},SAMPLE_IDS_STR="${SAMPLE_IDS_STR}",IN_DIR=${IN_DIR},OUT_DIR=${OUT_DIR},METADATA_FILE=${METADATA_FILE} STAR_align_174.sbatch

echo "Slurm array job submitted."

