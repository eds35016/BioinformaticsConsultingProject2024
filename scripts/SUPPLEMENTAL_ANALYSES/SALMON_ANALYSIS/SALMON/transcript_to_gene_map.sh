#!/bin/bash
#SBATCH -p standard
#SBATCH -c 1
#SBATCH --mem=10G
#SBATCH --job-name=tx2gene
#SBATCH -o tx2gene_output_%j.txt
#SBATCH -e tx2gene_error_%j.txt


# Take a GTF file and outputs unique combination of transcripts ID - gene ID pair
zless -S genomic.gtf | grep -v "#" | awk '$3=="transcript"' | cut -f9 | tr -s ";" " " | awk '{print$4"\t"$2}' | sort | uniq |  sed 's/\"//g' | tee tx2gene.tsv
