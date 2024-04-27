#!/bin/bash
#SBATCH --partition=standard
#SBACTH --output=/scratch/bioconsult/Ethan_James/BINF_Consulting_2024/ongoing/mapping/PRRSV2/logs/trimalign_out_%A.log
#SBATCH --error=/scratch/bioconsult/Ethan_James/BINF_Consulting_2024/ongoing/mapping/PRRSV2/logs/trimalign_err_%A.log
#SBATCH --cpus-per-task 4
#SBATCH --mem=50G

mkdir -p "/scratch/bioconsult/Ethan_James/BINF_Consulting_2024/ongoing/mapping/PRRSV2/PRRSV2/trimmed_data/cDC_Mock_4/cDC_Mock_4_peOverlap"

# Define variables
INPUT_DIR="/scratch/bioconsult/Ethan_James/BINF_Consulting_2024/raw_data/PRRSV2/cDC_Mock_4"
OUTPUT_DIR="/scratch/bioconsult/Ethan_James/BINF_Consulting_2024/ongoing/mapping/PRRSV2/PRRSV2/trimmed_data/cDC_Mock_4"
ALIGN_INPUT_DIR=${OUTPUT_DIR}

SAMPLE_ID="cDC_Mock_4"
RUN_NAME="PRRSV2"
GENOME_DIR="/scratch/bioconsult/Ethan_James/BINF_Consulting_2024/ongoing/mapping/STAR_idx/${RUN_NAME}/pig_only"

# Call conda environment. Change environment name as needed.
source ~/miniconda3/etc/profile.d/conda.sh
conda activate DSS

# Adapters unknown. 
# F_ADAPTER=""
# R_ADAPTER=""

# Pair-end sequencing files
R1=$(find ${INPUT_DIR} -name "${SAMPLE_ID}*_1.fq.gz" | head -n 1)
R2=$(find ${INPUT_DIR} -name "${SAMPLE_ID}*_2.fq.gz" | head -n 1)

# Run Trim Galore. -t n.  Do not have adapter sequence - automatically detect adapters.
trim_galore --paired ${R1} ${R2} -o ${OUTPUT_DIR} --adapter "AAGCAGTGGTATCAACGCAGAGTAC" \
--clip_r1 1 \
--clip_r2 1 \
--three_prime_clip_r1 30 \
--three_prime_clip_r2 30 \
--dont_gzip

# Set trimmed files for alignment
READ1=$(find ${ALIGN_INPUT_DIR} -name "${SAMPLE_ID}*_1_val_1.fq" | head -n 1)
READ2=$(find ${ALIGN_INPUT_DIR} -name "${SAMPLE_ID}*_2_val_2.fq" | head -n 1)

# Perform alignment on the trimmed reads.
/scratch/bioconsult/Ethan_James/BINF_Consulting_2024/tools/STAR-2.7.11b/bin/Linux_x86_64_static/STAR \
    --quantMode GeneCounts \
    --genomeDir $GENOME_DIR \
    --readFilesIn $READ1 $READ2 \
    --runThreadN $SLURM_CPUS_PER_TASK \
    --outFileNamePrefix ${OUTPUT_DIR}/${SAMPLE_ID}/${SAMPLE_ID}_clip30

# Perform alignment on the trimmed reads. with peOverlap
/scratch/bioconsult/Ethan_James/BINF_Consulting_2024/tools/STAR-2.7.11b/bin/Linux_x86_64_static/STAR \
    --quantMode GeneCounts \
    --genomeDir $GENOME_DIR \
    --readFilesIn $READ1 $READ2 \
    --runThreadN $SLURM_CPUS_PER_TASK \
    --peOverlapNbasesMin 10 \
    --outFileNamePrefix ${OUTPUT_DIR}/${SAMPLE_ID}_peOverlap/${SAMPLE_ID}_peOverlap_

date
exit 0
