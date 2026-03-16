#!/bin/bash
#SBATCH --job-name=bedtools_intersect_count
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=compute9

#Bedtools intersect command that will indicate how many cpg sites are present in the promoter windows of the respective promoter sources
DIR="/home/UNIXHOME/echia/SV_Anno/intersect_output/files_Venn/intersect_final"
COUNT_DIR="${DIR}/count"
#HG03579 cell line used
CELL_LINE="/pbi/vast-collections/appslabht/human_germline_wgs/hprc/HG03579/HG03579.GRCh38.cpg_pileup.combined.bed.gz"

#rm -r "$COUNT_DIR"
mkdir -p "$COUNT_DIR"
#All intersect
bedtools intersect \
        -a "${DIR}/venn_all(1v2v3).bed.gz" \
        -b "$CELL_LINE" \
	-c | gzip > "${COUNT_DIR}/venn_all.bed.gz" 

#GENCODE and ENCODE intersect regions
bedtools intersect \
        -a "${DIR}/GENCODE_encode(1v2).bed.gz" \
        -b "$CELL_LINE" \
        -c | gzip > "${COUNT_DIR}/GEN_EN.bed.gz"

#GENCODE and Ensembl intersect regions
bedtools intersect \
        -a "${DIR}/GENCODE_ensembl(1v3).bed.gz" \
        -b "$CELL_LINE" \
        -c | gzip > "${COUNT_DIR}/GEN_Ensem.bed.gz"

#ENCODE and Ensembl intersect regions
bedtools intersect \
        -a "${DIR}/ENCODE_ensembl(2v3).bed.gz" \
        -b "$CELL_LINE" \
        -c | gzip > "${COUNT_DIR}/EN_Ensem.bed.gz"

#GENCODE only regions
bedtools intersect \
        -a "${DIR}/GENCODE_only.bed.gz" \
        -b "$CELL_LINE" \
        -c | gzip > "${COUNT_DIR}/GEN_only.bed.gz"

#ENCODE only regions
bedtools intersect \
        -a "${DIR}/ENCODE_only.bed.gz" \
        -b "$CELL_LINE" \
        -c | gzip > "${COUNT_DIR}/EN_only.bed.gz"

#Ensembl only regions
bedtools intersect \
        -a "${DIR}/Ensembl_only.bed.gz" \
        -b "$CELL_LINE" \
        -c | gzip > "${COUNT_DIR}/Ensem_only.bed.gz"
