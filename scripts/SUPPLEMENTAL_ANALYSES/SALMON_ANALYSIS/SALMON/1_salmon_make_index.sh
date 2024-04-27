#!/bin/bash
#SBATCH -p standard
#SBATCH -o salmon-%j.out
#SBATCH --exclusive
#SBATCH --mem=0

/scratch/bioconsult/Ethan_James/BINF_Consulting_2024/tools/salmon-latest_linux_x86_64/bin/salmon index -t /scratch/bioconsult/Ethan_James/BINF_Consulting_2024/raw_data/ref_genomes/pig/gentrome.fna -d /scratch/bioconsult/Ethan_James/BINF_Consulting_2024/raw_data/ref_genomes/pig/decoys.txt -p 12 -i /scratch/bioconsult/Ethan_James/BINF_Consulting_2024/ongoing/mapping/salmon_index
