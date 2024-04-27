#!/bin/bash
#SBATCH --job-name=salmon_align
#SBATCH --output=/scratch/bioconsult/Ethan_James/BINF_Consulting_2024/ongoing/mapping/salmon_out/PRRSV2/logs/salmon_align_out_%A_%a.log
#SBATCH --error=/scratch/bioconsult/Ethan_James/BINF_Consulting_2024/ongoing/mapping/salmon_out/PRRSV2/logs/salmon_align_err_%A_%a.log
#SBATCH --cpus-per-task=4
#SBATCH --mem=40G
#SBATCH -p standard

# Get the sample ID from the array
SAMPLE_ID=($(echo $SAMPLE_IDS_STR | tr " " "\n"))
SAMPLE_ID=${SAMPLE_ID[$SLURM_ARRAY_TASK_ID - 1]}
echo "Processing ${SAMPLE_ID}..."

# Define the directory containing the fastq files
FASTQ_DIR="$IN_DIR/${SAMPLE_ID}"

# Construct the fastq filenames
READ1_FILES=(${FASTQ_DIR}/${SAMPLE_ID}*"_1.fq.gz")
READ2_FILES=(${FASTQ_DIR}/${SAMPLE_ID}*"_2.fq.gz")

# Check if exactly one file is matched for each read
if [ ${#READ1_FILES[@]} -ne 1 ] || [ ${#READ2_FILES[@]} -ne 1 ]; then
    echo "Input files for $SAMPLE_ID do not exist or multiple matches found. Exiting."
    exit 1
fi

# If exactly one match is found for each, assign to READ1 and READ2
READ1=${READ1_FILES[0]}
READ2=${READ2_FILES[0]}

echo $READ1
echo $READ2

# Check if the input files exist
if [ ! -f "$READ1" ] || [ ! -f "$READ2" ]; then
    echo "Input files for $SAMPLE_ID do not exist. Exiting."
    exit 1
fi

# Create output directory for this sample
mkdir -p "$OUT_DIR/$RUN_NAME/$SAMPLE_ID"

# Run salmon alignment
/scratch/bioconsult/Ethan_James/BINF_Consulting_2024/tools/salmon-latest_linux_x86_64/bin/salmon quant -i $GENOME_DIR -l A -1 $READ1 -2 $READ2 --validateMappings -o $OUT_DIR/${SAMPLE_ID}

echo "Alignment done for $SAMPLE_ID"
