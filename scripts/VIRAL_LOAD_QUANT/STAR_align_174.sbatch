#!/bin/bash
#SBATCH --job-name=STAR_align_174
#SBATCH --output=/scratch/bioconsult/Ethan_James/BINF_Consulting_2024/ongoing/mapping/PRRSV2_174_QUANT/logs/STAR_align_out_%A_%a.log
#SBATCH --error=/scratch/bioconsult/Ethan_James/BINF_Consulting_2024/ongoing/mapping/PRRSV2_174_QUANT/logs/STAR_align_err_%A_%a.log
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

# Run STAR alignment. Specified peOverlap parameter.
/scratch/bioconsult/Ethan_James/BINF_Consulting_2024/tools/STAR-2.7.11b/bin/Linux_x86_64_static/STAR --quantMode GeneCounts \
    --genomeDir $GENOME_DIR \
    --readFilesIn $READ1 $READ2 \
    --readFilesCommand zcat \
    --runThreadN $SLURM_CPUS_PER_TASK \
    --peOverlapNbasesMin 10 \
    --outFileNamePrefix ${OUT_DIR}/${RUN_NAME}/${SAMPLE_ID}/${RUN_NAME}_${SAMPLE_ID}_

echo "Alignment done for $SAMPLE_ID"
