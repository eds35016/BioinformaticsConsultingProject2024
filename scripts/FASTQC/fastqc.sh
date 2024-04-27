#!/bin/bash
#SBATCH --partition=standard
#SBACTH --output=fastqc_%A.out
#SBATCH --cpus-per-task 4




SAMPLE="moDC_Mock_6"

INPUT_DIR="/scratch/bioconsult/Ethan_James/BINF_Consulting_2024/raw_data/PRRSV2/${SAMPLE}"
OUTPUT_DIR="/scratch/bioconsult/Ethan_James/BINF_Consulting_2024/ongoing/QC/PRRSV2/${SAMPLE}"
mkdir -p ${OUTPUT_DIR}



# Call conda environment. Change environment name as needed.
conda activate /scratch/bioconsult/Ethan_James/BINF_Consulting_2024/tools/r_env

R1=$(find ${INPUT_DIR} -name ${SAMPLE}*_1.fq.gz |head -n 1)
R2=$(find ${INPUT_DIR} -name ${SAMPLE}*_2.fq.gz |head -n 1)
echo "R1 is ${R1}"
echo "R2 is ${R2}"
# run fastqc. -o output directory. -t number of cores, -noextract for no GUI
fastqc -o ${OUTPUT_DIR} -t 4 -noextract ${R1}
fastqc -o ${OUTPUT_DIR} -t 4 -noextract ${R2}

exit 0
