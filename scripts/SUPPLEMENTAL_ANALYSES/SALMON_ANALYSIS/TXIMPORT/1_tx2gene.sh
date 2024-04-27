#!/bin/bash
#SBATCH --job-name=tximport
#SBATCH --output=//scratch/bioconsult/Ethan_James/BINF_Consulting_2024/ongoing/mapping/salmon_out/logs/tximport_out_%A_%a.log
#SBATCH --error=/scratch/bioconsult/Ethan_James/BINF_Consulting_2024/ongoing/mapping/salmon_out/logs/tximport_err_%A_%a.log
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --array=0-22%11

# Define variables
RUN_NAME="PRRSV2"
RUN_ID="CKDL210021614-1a-AK7207-AK7213_HMH2WDSX2"
GENOME_DIR="/scratch/bioconsult/Ethan_James/BINF_Consulting_2024/ongoing/mapping/salmon_index"
IN_DIR="/scratch/bioconsult/Ethan_James/BINF_Consulting_2024/ongoing/mapping/salmon_out/${RUN_NAME}"
OUT_DIR="/scratch/bioconsult/Ethan_James/BINF_Consulting_2024/ongoing/mapping/salmon_out/${RUN_NAME}"
METADATA_FILE="/scratch/bioconsult/Ethan_James/BINF_Consulting_2024/metadata/PRRSV2_metadata.csv"
MAP_FILE="/scratch/bioconsult/Ethan_James/BINF_Consulting_2024/raw_data/ref_genomes/pig/tx2gene.tsv"

# Create output directory if needed
mkdir -p $OUT_DIR

# Extract sample IDs
SAMPLE_IDS=($(awk -F ',' 'NR>1 {print $2}' ${METADATA_FILE}))

# Set Sample ID based on the array task ID
SAMPLE_ID=${SAMPLE_IDS[$SLURM_ARRAY_TASK_ID]}
IN_FILE="${IN_DIR}/${RUN_NAME}/${SAMPLE_ID}/quant.sf"
QUANT_OUT_DIR="${OUT_DIR}/${RUN_NAME}/${SAMPLE_ID}/"
mkdir -p ${QUANT_OUT_DIR}

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
export IN_DIR
export OUT_DIR
export SAMPLE_IDS
export METADATA_FILE
export QUANT_OUT_DIR
export MAP_FILE


# Create a string with comma-separated sample IDs for exporting
export SAMPLE_IDS_STR="${SAMPLE_IDS[*]}"

# load conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate R


# Running R with input and output DIRs and the transcripts-to-map file
Rscript ./tx2gene.R $IN_FILE $QUANT_OUT_DIR $MAP_FILE

echo "Slurm array job submitted."

