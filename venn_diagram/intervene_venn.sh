#!/bin/bash
#SBATCH --job-name=bedtools_intersect
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=compute9

HOME_DIR="/home/UNIXHOME/echia/SV_Anno/"
GENCODE="${HOME_DIR}/promoter/GENCODE/GENCODE_clean_original.gff3.gz"
ENCODE="${HOME_DIR}/promoter/ENCODE/ENCODE_original.bed.gz"
Ensembl="${HOME_DIR}/promoter/Ensembl/Ensembl_promoters_clean.gff3.gz"
OUTPUT_DIR="${HOME_DIR}/intersect_output/diagrams/venn_results"
mkdir -p "$OUTPUT_DIR"
PERCENT="${1}"

intervene venn -i \
	"$GENCODE" \
	"$ENCODE" \
	"$Ensembl" \
	--names=GENCODE,ENCODE,Ensembl \
	--bedtools-options f="$PERCENT" \
	--output "${OUTPUT_DIR}/Venn_${PERCENT}"

