#!/bin/bash
#SBATCH --job-name=venn_diagram_results
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=compute9

INTERSECT_DIR="/home/UNIXHOME/echia/SV_Anno/intersect_output"
OUTPUT_DIR="${INTERSECT_DIR}/files_Venn/intersect_final"
PROMOTER_DIR="/home/UNIXHOME/echia/SV_Anno/promoter"
rm -r "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR"

# Code used to create files from GENCODE, ENCODE and Ensembl that intersect with each other
# Files are first merged together using bedtools merge and then bedtools intersect is used to differentiate which promoter windows intersect with which regions

GENCODE="${PROMOTER_DIR}/GENCODE/GENCODE_clean_sorted.bed"
ENCODE="${PROMOTER_DIR}/ENCODE/ENCODE_clean_sorted.bed"
Ensembl="${PROMOTER_DIR}/Ensembl/Ensembl_clean_sorted.bed"
Overall_bed="${OUTPUT_DIR}/all_merged_regions.bed.gz"
COMBINED_TMP="${OUTPUT_DIR}/combined_temp.bed"

cat \
	<(grep -v '^#' "$GENCODE" | grep -v '^track' | cut -f1-6) \
	<(grep -v '^#' "$ENCODE" | grep -v '^track' | cut -f1-6) \
	<(grep -v '^#' "$Ensembl" | grep -v '^track'| cut -f1-6) \
| sort -k1,1 -k2,2n > "$COMBINED_TMP"

echo "Lines in combined file: $(wc -l < "$COMBINED_TMP")"

bedtools merge -c 4,1 -o collapse,count \
	-i "$COMBINED_TMP" | gzip > "$Overall_bed"
echo -n "Total segments:" && zcat "$Overall_bed" | wc -l

#Overall intersecting parts (GENCODE overlap with ENCODE with Ensembl)
bedtools intersect -a "$Overall_bed" -b "$GENCODE" -u | \
bedtools intersect -a stdin -b "$Ensembl" -u | \
bedtools intersect -a stdin -b "$ENCODE" -u > "${OUTPUT_DIR}/venn_all(1v2v3).bed" \
&& gzip "${OUTPUT_DIR}/venn_all(1v2v3).bed" 
echo -n "All intersecting parts:" && zcat "${OUTPUT_DIR}/venn_all(1v2v3).bed.gz" | wc -l

#GENCODE intersect with ENCODE only
bedtools intersect -a "$Overall_bed" -b "$GENCODE" -u | \
bedtools intersect -a stdin -b "$ENCODE" -u | \
bedtools intersect -a stdin -b "$Ensembl" -v > "${OUTPUT_DIR}/GENCODE_encode(1v2).bed" \
&& gzip "${OUTPUT_DIR}/GENCODE_encode(1v2).bed" 
echo -n "GENCODE intersect with ENCODE only:" && zcat "${OUTPUT_DIR}/GENCODE_encode(1v2).bed.gz" | wc -l

#GENCODE intersect with Ensembl only
bedtools intersect -a "$Overall_bed" -b "$GENCODE" -u | \
bedtools intersect -a stdin -b "$Ensembl" -u | \
bedtools intersect -a stdin -b "$ENCODE" -v > "${OUTPUT_DIR}/GENCODE_ensembl(1v3).bed" \
&& gzip "${OUTPUT_DIR}/GENCODE_ensembl(1v3).bed"
echo -n "GENCODE intersect with Ensembl only:" && zcat "${OUTPUT_DIR}/GENCODE_ensembl(1v3).bed.gz" | wc -l

#ENCODE intersect with Ensembl only
bedtools intersect -a "$Overall_bed" -b "$ENCODE" -u | \
bedtools intersect -a stdin -b "$Ensembl" -u | \
bedtools intersect -a stdin -b "$GENCODE" -v > "${OUTPUT_DIR}/ENCODE_ensembl(2v3).bed" \
&& gzip "${OUTPUT_DIR}/ENCODE_ensembl(2v3).bed" 
echo -n "ENCODE intersect with Ensembl only:" && zcat "${OUTPUT_DIR}/ENCODE_ensembl(2v3).bed.gz" | wc -l

#GENCODE only
bedtools intersect -a "$Overall_bed" -b "$GENCODE" -u | \
bedtools intersect -a stdin -b "$Ensembl" -v | \
bedtools intersect -a stdin -b "$ENCODE" -v > "${OUTPUT_DIR}/GENCODE_only.bed" \
&& gzip "${OUTPUT_DIR}/GENCODE_only.bed" 
echo -n "GENCODE only:" && zcat "${OUTPUT_DIR}/GENCODE_only.bed.gz" | wc -l

#ENCODE only
bedtools intersect -a "$Overall_bed" -b "$ENCODE" -u | \
bedtools intersect -a stdin -b "$Ensembl" -v | \
bedtools intersect -a stdin -b "$GENCODE" -v > "${OUTPUT_DIR}/ENCODE_only.bed" \
&& gzip "${OUTPUT_DIR}/ENCODE_only.bed"
echo -n "ENCODE only:" && zcat "${OUTPUT_DIR}/ENCODE_only.bed.gz" | wc -l

#Ensembl only
bedtools intersect -a "$Overall_bed" -b "$Ensembl" -u | \
bedtools intersect -a stdin -b "$ENCODE" -v | \
bedtools intersect -a stdin -b "$GENCODE" -v > "${OUTPUT_DIR}/Ensembl_only.bed" \
&& gzip "${OUTPUT_DIR}/Ensembl_only.bed" 
echo -n "Ensembl only:" && zcat "${OUTPUT_DIR}/Ensembl_only.bed.gz" | wc -l




