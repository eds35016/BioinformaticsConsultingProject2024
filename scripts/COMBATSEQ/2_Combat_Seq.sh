#!/bin/bash
#SBATCH --job-name=combatseq
#SBATCH --output=/scratch/bioconsult/Ethan_James/BINF_Consulting_2024/ongoing/mapping/salmon_out/logs/COMBATSEQ_out_%A_%a.log
#SBATCH --error=/scratch/bioconsult/Ethan_James/BINF_Consulting_2024/ongoing/mapping/salmon_out/logs/COMBATSEQ_err_%A_%a.log
#SBATCH --mem=10G
#SBATCH --cpus-per-task=2

# Define variables
RUN_NAME="PRRSV2"
IN_DIR="/scratch/bioconsult/Ethan_James/BINF_Consulting_2024/ongoing/mapping/salmon_out/${RUN_NAME}/${RUN_NAME}" # For Salmon input
OUT_DIR="/scratch/bioconsult/Ethan_James/BINF_Consulting_2024/ongoing/mapping/salmon_out/${RUN_NAME}/${RUN_NAME}" # For Salmon output

#IN_DIR="/scratch/bioconsult/Ethan_James/BINF_Consulting_2024/ongoing/mapping/${RUN_NAME}/${RUN_NAME}" #For STAR input
#OUT_DIR="/scratch/bioconsult/Ethan_James/BINF_Consulting_2024/ongoing/mapping/${RUN_NAME}/${RUN_NAME}" #For STAR input

COUNT_FILE="${IN_DIR}/combined_gene_counts.csv"

# Create output directory if needed
mkdir -p ${OUT_DIR}


# Check if count file exists
if [ ! -f ${COUNT_FILE} ]; then
	echo "Combined Counts File does not exist. Exiting."
	exit 1
fi

echo "Processing Combined Counts File"
export IN_DIR
export OUT_DIR
export COUNT_FILE

# load conda environment
#source ~/miniconda3/etc/profile.d/conda.sh
conda activate /scratch/bioconsult/Ethan_James/BINF_Consulting_2024/tools/r_env

# Run Rscript to execute combat-seq
/scratch/bioconsult/Ethan_James/BINF_Consulting_2024/tools/R-4.3.3/bin/Rscript combatseq.R

echo "SLURM job sumitted."
