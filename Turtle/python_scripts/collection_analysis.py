#!/usr/bin/env python3

import pandas as pd
import numpy as np
import argparse
import sys
import os
import glob
from typing import Dict, List, Tuple
from collections import defaultdict
import gc

def process_single_chromosome(collection_file: str, output_file: str, chromosome: str):
    print(f"=== Processing chromosome {chromosome} ===")

    # Read collection
    df_collection = pd.read_csv(collection_file, sep='\t')
    file_collection = {}

    for _, row in df_collection.iterrows():
        identifier = row['identifier']
        base_path = row['filename'].strip()

        bed_files = glob.glob(f"{base_path}.combined.bed.gz")
        if bed_files:
            file_collection[identifier] = bed_files[0]
            print(f"Found file for {identifier}: {os.path.basename(bed_files[0])}")

    if not file_collection:
        print("❌ No valid files found in collection")
        return

    print(f"Processing {len(file_collection)} samples for chromosome {chromosome}")

    # Aggregate data for this chromosome
    all_position_data = defaultdict(list)
    processed_samples = 0

    for identifier, bed_file in file_collection.items():
        try:
            processed_samples += 1
            print(f"[{processed_samples}/{len(file_collection)}] Processing {identifier}...")

            # Read with chromosome filter
            dtype = {'chrom': 'category', 'start': 'int32', 'end': 'int32', 'score': 'float32'}
            df = pd.read_csv(bed_file, sep='\t', header=None,
                           usecols=[0, 1, 2, 3],
                           names=['chrom', 'start', 'end', 'score'],
                           dtype=dtype)

            # Filter by chromosome
            df_chrom = df[df['chrom'] == chromosome].copy()

            if len(df_chrom) == 0:
                print(f"  No data for {chromosome} in {identifier}")
                continue

            print(f"  Found {len(df_chrom)} CpG sites")

            # Aggregate scores
            sites_processed = 0
            for _, row in df_chrom.iterrows():
                pos_key = (row['chrom'], row['start'], row['end'])
                all_position_data[pos_key].append(row['score'])
                sites_processed += 1

            print(f"  Processed {sites_processed} sites")

        except Exception as e:
            print(f"❌ Error processing {bed_file}: {e}")

    print(f"\nFound {len(all_position_data)} unique CpG positions on {chromosome}")

    # Compute statistics
    if not all_position_data:
        print("❌ No data to process")
        return

    print("Computing statistics...")
    results = []

    batch_size = 10000
    position_items = list(all_position_data.items())
    total_batches = (len(position_items) + batch_size - 1) // batch_size

    for batch_idx in range(total_batches):
        start_idx = batch_idx * batch_size
        end_idx = min((batch_idx + 1) * batch_size, len(position_items))

        print(f"  Batch {batch_idx + 1}/{total_batches}...")

        for (chrom, start, end), scores in position_items[start_idx:end_idx]:
            scores_array = np.array(scores, dtype=np.float32)

            mean_score = np.mean(scores_array)
            std_score = np.std(scores_array)
            variance_score = np.var(scores_array)
            median = np.median(scores_array)
            range_score = np.max(scores_array) - np.min(scores_array)

            if std_score > 0:
                z_scores = (scores_array - mean_score) / std_score
                max_zscore = np.max(np.abs(z_scores))
            else:
                max_zscore = 0.0

            results.append({
                'chrom': chrom,
                'start': start,
                'end': end,
                'mean': float(mean_score),
                'std': float(std_score),
                'variance': float(variance_score),
                'max_zscore': float(max_zscore),
                'median': float(median),
                'range': float(range_score),
                'no_of_samples': len(scores)
            })

    # Save results
    results_df = pd.DataFrame(results)
    if len(results_df) > 0:
        results_df.to_csv(output_file, sep='\t', index=False)
        print(f"Saved {len(results_df)} sites to {output_file}")

        # Summary
        avg_coverage = results_df['no_of_samples'].mean()
        print(f"=== Summary for {chromosome} === ")
        print(f"CpG sites: {len(results_df):,}")
        print(f"Average coverage: {avg_coverage:.1f} samples/site")
    else:
        print(f"ERROR! No results to save for {chromosome}")

def main():
    parser = argparse.ArgumentParser(description='Process CpG data for specific chromosomes')
    parser.add_argument('-c', '--collection', required=True, help='Collection TSV file')
    parser.add_argument('-o', '--output', required=True, help='Output file')
    parser.add_argument('--chromosome', required=True, help='Chromosome to process (e.g., chr1, chrX)')
    args = parser.parse_args()

    process_single_chromosome(args.collection, args.output, args.chromosome)

if __name__ == "__main__":
    main()



