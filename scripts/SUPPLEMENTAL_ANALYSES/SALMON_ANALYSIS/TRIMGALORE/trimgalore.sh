#!/bin/bash
#SBATCH --partition=standard
#SBACTH --output=/scratch/bioconsult/Ethan_James/BINF_Consulting_2024/ongoing/mapping/PRRSV2/logs/trimmed_data_out_%A.log
#SBATCH --error=/scratch/bioconsult/Ethan_James/BINF_Consulting_2024/ongoing/mapping/PRRSV2/logs/trimmed_data_err_%A.log
#SBATCH --cpus-per-task 4

mkdir -p "/scratch/bioconsult/Ethan_James/BINF_Consulting_2024/ongoing/mapping/PRRSV2/PRRSV2/trimmed_data/cDC_134_3"
# Define variables
INPUT_DIR="/scratch/bioconsult/Ethan_James/BINF_Consulting_2024/raw_data/PRRSV2/cDC_134_3"
OUTPUT_DIR="/scratch/bioconsult/Ethan_James/BINF_Consulting_2024/ongoing/mapping/PRRSV2/PRRSV2/trimmed_data"

# Call conda environment. Change environment name as needed.
source ~/miniconda3/etc/profile.d/conda.sh
conda activate DSS

# Adapters unknown. 
# F_ADAPTER=""
# R_ADAPTER=""

# Pair-end sequencing files
R1="${INPUT_DIR}/cDC_134_3_CKDL210021614-1a-AK7209-AK7215_HMH2WDSX2_L4_1.fq.gz"
R2="${INPUT_DIR}/cDC_134_3_CKDL210021614-1a-AK7209-AK7215_HMH2WDSX2_L4_2.fq.gz"

# Run Trim Galore. -t n.  Do not have adapter sequence - automatically detect adapters  # -j sets parallelization. It  requiers 'pigz' installed.
trim_galore --paired ${R1} ${R2} --dont_gzip -o ${OUTPUT_DIR}  #-j 4  #/ --adapter ${F_ADAPTER} --adapter2 ${R_ADAPTER} 

date
exit 0
