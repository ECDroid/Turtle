#!/bin/bash
#SBATCH --job-name=bedtools_intersect
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=compute9

INTERSECT_DIR="/home/UNIXHOME/echia/SV_Anno/intersect_output"
OUTPUT_DIR="${INTERSECT_DIR}/Venn_diagrams/Venn_overall"
rm -r "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR"

GENCODE="${INTERSECT_DIR}/count/GENCODE_HG03579_inter_count.bed.gz"
ENCODE="${INTERSECT_DIR}/count/ENCODE_HG03579_inter_count.bed.gz"
Ensembl="${INTERSECT_DIR}/count/Ensembl_HG03579_inter_count.bed.gz"
#Intersect commands where 90% of a GENCODE segment must overlap with Ensembl/ENCODE segment

#GENCODE vs ENCODE (1v2)
bedtools intersect \
	-a "$GENCODE" \
	-b "$ENCODE" \
	-wa | \
bedtools intersect \
	-a stdin \
	-b "$Ensembl" \
	-v > "${OUTPUT_DIR}/GENCODE_ENCODE(1v2).bed" \
&& gzip "${OUTPUT_DIR}/GENCODE_ENCODE(1v2).bed" \
&& zcat "${OUTPUT_DIR}/GENCODE_ENCODE(1v2).bed.gz" | wc -l

#GENCODE vs Ensembl (1v3)
bedtools intersect \
	-a "$GENCODE" \
	-b "$Ensembl" \
	-wa | \
bedtools intersect \
	-a stdin \
	-b "$ENCODE" \
	-v > "${OUTPUT_DIR}/GENCODE_Ensembl(1v3).bed" \
&& gzip "${OUTPUT_DIR}/GENCODE_Ensembl(1v3).bed" \
&& zcat "${OUTPUT_DIR}/GENCODE_Ensembl(1v3).bed.gz" | wc -l

#ENCODE vs Ensembl (2v3)
bedtools intersect \
	-a "$ENCODE" \
	-b "$Ensembl" \
	-wa | \
bedtools intersect \
	-a stdin \
	-b "$GENCODE" \
	-v > "${OUTPUT_DIR}/ENCODE_Ensembl(2v3).bed" \
&& gzip "${OUTPUT_DIR}/ENCODE_Ensembl(2v3).bed" \
&& zcat "${OUTPUT_DIR}/ENCODE_Ensembl(2v3).bed.gz" | wc -l

#GENCODE vs Ensembl vs ENCODE (1v2v3)
bedtools intersect \
        -a "$GENCODE" \
        -b "$ENCODE" \
	-wa | \
bedtools intersect \
	-a stdin \
	-b "$Ensembl" \
        -wa > "${OUTPUT_DIR}/all_overlap(1v2v3).bed" \
&& gzip "${OUTPUT_DIR}/all_overlap(1v2v3).bed" \
&& zcat "${OUTPUT_DIR}/all_overlap(1v2v3).bed.gz" | wc -l

#Test Ensembl vs GENCODE vs ENCODE (1v2v3)
bedtools intersect \
        -a "$Ensembl" \
        -b "$ENCODE" \
        -wa | \
bedtools intersect \
        -a stdin \
        -b "$GENCODE" \
        -wa > "${OUTPUT_DIR}/all_overlap(1v2v3)_test.bed" \
&& gzip "${OUTPUT_DIR}/all_overlap(1v2v3)_test.bed" \
&& zcat "${OUTPUT_DIR}/all_overlap(1v2v3)_test.bed.gz" | wc -l

#Only in GENCODE (1-2-3)
bedtools intersect \
        -a "$GENCODE" \
        -b "$ENCODE" \
        "$Ensembl" \
	-wa -v > "${OUTPUT_DIR}/GENCODE(only).bed" \
&& gzip "${OUTPUT_DIR}/GENCODE(only).bed" \
&& zcat "${OUTPUT_DIR}/GENCODE(only).bed.gz" | wc -l

#Only in ENCODE (2-1-3)
bedtools intersect \
        -a "$ENCODE" \
        -b "$GENCODE" \
        "$Ensembl" \
        -wa -v > "${OUTPUT_DIR}/ENCODE(only).bed" \
&& gzip "${OUTPUT_DIR}/ENCODE(only).bed" \
&& zcat "${OUTPUT_DIR}/ENCODE(only).bed.gz" | wc -l

#Only in Ensembl (3-1-2)
bedtools intersect \
        -a "$Ensembl" \
        -b "$GENCODE" \
        "$ENCODE" \
        -wa -v > "${OUTPUT_DIR}/Ensembl(only).bed" \
&& gzip "${OUTPUT_DIR}/Ensembl(only).bed" \
&& zcat "${OUTPUT_DIR}/Ensembl(only).bed.gz" | wc -l

