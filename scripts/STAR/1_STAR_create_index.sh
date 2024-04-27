#!/bin/bash
#SBATCH --job-name=STAR_create_index
#SBATCH --output=/scratch/bioconsult/Ethan_James/BINF_Consulting_2024/ongoing/mapping/STAR_idx/PRRSV2/logs/STAR_create_index_out_%A_%a.log
#SBATCH --error=/scratch/bioconsult/Ethan_James/BINF_Consulting_2024/ongoing/mapping/STAR_idx/PRRSV2/logs/STAR_create_index_err_%A_%a.log
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH -p standard

# set up reference genome index

        star="/scratch/bioconsult/Ethan_James/BINF_Consulting_2024/tools/STAR-2.7.11b/bin/Linux_x86_64_static/STAR"

        refseqdir="/scratch/bioconsult/Ethan_James/BINF_Consulting_2024/raw_data/ref_genomes/pig"

        idxdir="/scratch/bioconsult/Ethan_James/BINF_Consulting_2024/ongoing/mapping/STAR_idx/PRRSV2/pig_only"

        $star --runThreadN 8 --runMode genomeGenerate --genomeDir $idxdir \
                --genomeFastaFiles $refseqdir/GCF_000003025.6_Sscrofa11.1_genomic.fna \
                --sjdbGTFfile $refseqdir/genomic.gtf \
                --sjdbOverhang 149
