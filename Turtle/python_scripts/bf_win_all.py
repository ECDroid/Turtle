#!/usr/bin/env python3

import pandas as pd
import numpy as np
import argparse
import sys
import gzip
import os
from typing import List, Dict, Tuple

def count_GEN_genes(gene_info: str) -> int:
    """Count number of GEN genes in the gene_info string"""
    if pd.isna(gene_info):
        return 0
    count = 0
    for item in gene_info.split(','):
        if item.strip().startswith('GEN('):
            count += 1
    return count

# Function to read the CpG file input
def read_cpg_sites(input_file: str, total_samples: int = 195) -> pd.DataFrame:
    """Read CpG sites with additional statistical columns"""
    df = pd.read_csv(input_file, sep='\t', header=None,
                    usecols=[0, 1, 2, 3, 4, 5, 6, 9],  # Added std, variance, no_of_samples
                    names=['chrom', 'start', 'end', 'mean', 'median', 'std', 'variance', 'no_of_samples'])
 
    # Normalise mean and median values and store them  as lamda values
    df['mean_lamda'] = df['mean'].apply(lambda x: 100 - x if x > 50 else x)
    df['median_lamda'] = df['median'].apply(lambda x: 100 - x if x > 50 else x)

    # Calculate trustworthiness index for each CpG
    df['trustworthiness'] = (df['no_of_samples'] / total_samples) * (1 / (1 + df['variance']))

    df = df.sort_values(['chrom', 'start']).reset_index(drop=True)
    return df

# Function to read promoter windows from input file
def promoter_windows(input_file: str) -> pd.DataFrame:
    df = pd.read_csv(input_file, sep='\t', header=None,
                    usecols=[0, 1, 2, 3],
                    names=['chrom', 'start', 'end', 'gene_info'])
    df = df.sort_values(['chrom', 'start']).reset_index(drop=True)
    return df

# Function to read individual promoters and their strand information (GENCODE file)
def read_individual_promoters(input_file: str) -> pd.DataFrame:
    # Read only the relevant columns, skip header rows (starting with 'track' or '#')
    df = pd.read_csv(input_file, sep='\t',
                    comment='#',  # Skip lines starting with #
                    skiprows=1,   # Skip the track line
                    usecols=[0, 1, 2, 3, 5],  # chrom, start, end, gene_name, strand
                    names=['chrom', 'start', 'end', 'gene_name', 'strand'])

    df = df.sort_values(['chrom', 'start']).reset_index(drop=True)
    return df

# Function to call TSS location from GENCODE promoter file
def get_tss_position(promoter_row: pd.Series) -> int:
    if promoter_row['strand'] == '+':
        return promoter_row['end']    # TSS at end for positive strand
    else:
        return promoter_row['start']  # TSS at start for negative strand

# Function to call promoter regions that overlap with a large promoter region
def find_individual_promoters_in_merged(merged_promoter: pd.Series, individual_promoters_df: pd.DataFrame) -> List[pd.Series]:
    chrom, start, end = merged_promoter['chrom'], merged_promoter['start'], merged_promoter['end']

    overlapping_individuals = individual_promoters_df[
        (individual_promoters_df['chrom'] == chrom) &
        (individual_promoters_df['start'] <= end) &
        (individual_promoters_df['end'] >= start)
    ]

    return [row for _, row in overlapping_individuals.iterrows()]

# Function to identify overlappping regions between individual promoters
def find_overlap_region(individual_promoters: List[pd.Series]) -> Tuple[int, int]:
    if len(individual_promoters) < 2:
        return None, None

    # Get TSS positions for each promoter
    tss_positions = [get_tss_position(p) for p in individual_promoters]

    # Find maximum start and minimum end among all promoters
    overlap_start = min(tss_positions)
    overlap_end = max(tss_positions)

    if overlap_start < overlap_end:
        return overlap_start, overlap_end
    else:
        return None, None

# Function to identify if a promoter window from the cost function has at least 1 boundary at either the start or end of the promtoer that is within the overlap region
def window_has_boundary_in_overlap(window: pd.Series, overlap_start: int, overlap_end: int) -> bool:
    if overlap_start is None or overlap_end is None:
        return True  # No overlap region defined, accept all windows
    
    start_in_overlap = overlap_start <= window['start'] <= overlap_end
    end_in_overlap = overlap_start <= window['end'] <= overlap_end

    return start_in_overlap or end_in_overlap 

# Function to generate statistics for the promoter windows
def window_stats_calculation(original_scores: List[float], lamda_scores: List[float], trustworthiness_scores: List[float] = None) -> Dict[str, float]:
    if len(original_scores) == 0:
        return {'std': 0.0, 'variance': 0.0, 'mean': 0.0, 'median': 0.0, 'weighted_std': 0.0, 'weighted_variance': 0.0, 'weighted_mean': 0.0, 'lamda_mean': 0.0}

    scores_array = np.array(original_scores)
    lamda_array = np.array(lamda_scores)

    window_mean = np.mean(scores_array)
    if window_mean > 50:
        lamda_window_based = 100 - window_mean
    else:
        lamda_window_based = window_mean

    # Calculate average trustworthiness if provided
    if trustworthiness_scores and len(trustworthiness_scores) > 0:
        weights_array = np.array(trustworthiness_scores)

        # Calculate weighted statistics
        weighted_mean_original = np.average(scores_array, weights=weights_array)
        weighted_variance = np.average((scores_array - weighted_mean_original)**2, weights=weights_array)
        weighted_std = np.sqrt(weighted_variance)
        lamda_mean = np.mean(lamda_array)

        return {
            'mean': np.mean(scores_array),
            'std': np.std(scores_array),
            'variance': np.var(scores_array),
            'weighted_std': weighted_std,
            'weighted_variance': weighted_variance,
            'weighted_mean': weighted_mean_original,
            'lamda_mean': lamda_mean,
            'lamda_window_based': lamda_window_based,
            'median': np.median(scores_array)
        }
    else:
        return {
            'std': np.std(scores_array),
            'variance': np.var(scores_array),
            'mean': np.mean(scores_array),
            'median': np.median(scores_array),
            'lamda_mean': np.mean(lamda_array),
            'lamda_window_based': lamda_window_based,
            'weighted_std': 0.0,
            'weighted_variance': 0.0,
            'weighted_mean': 0.0
        }

# Function to call specific cost functions for comparison
def calculate_cost_function(stats: Dict, genomic_length: int, window_size_cpgs: int, cost_function: str) -> float:
    """Calculate specific cost function based on user choice"""
    if cost_function == 'cost_1':
        return (stats['mean'] + stats['std']) / genomic_length
    elif cost_function == 'cost_2':
        return (stats['mean'] + stats['variance']) / genomic_length
    elif cost_function == 'cost_3':
        return stats['variance'] / genomic_length
    elif cost_function == 'cost_4':
        return (stats['mean'] + stats['weighted_std']) / genomic_length
    elif cost_function == 'cost_5':
        return (stats['mean'] + stats['weighted_variance']) / genomic_length
    elif cost_function == 'cost_6':
        return (stats['lamda_window_based'] + stats['weighted_variance']) / window_size_cpgs
    elif cost_function == 'cost_7':
        return (stats['lamda_mean'] + stats['weighted_variance']) / window_size_cpgs
    elif cost_function == 'cost_8':
        return (stats['mean'] + stats['weighted_variance']) / window_size_cpgs
    elif cost_function == 'cost_9':
        return (stats['lamda_mean'] + stats['weighted_variance']) / window_size_cpgs
    elif cost_function == 'cost_10':
        return (stats['lamda_mean'] + stats['weighted_variance']) / window_size_cpgs
    else:
        raise ValueError(f"Unknown cost function: {cost_function}")

# Function to perform sliding window analysis on a single promoter region
def analyze_single_promoter(promoter_row: pd.Series, cpg_df: pd.DataFrame, cost_function: str,
                           step_size: int = 1, individual_promoters_df: pd.DataFrame = None) -> pd.DataFrame:
    """Perform sliding window analysis for a single promoter region"""
    results = []

    # Extract CpG sites within this promoter region
    chrom, start, end, gene_info = promoter_row['chrom'], promoter_row['start'], promoter_row['end'], promoter_row['gene_info']

    # Count GEN genes in this promoter
    num_GEN_genes = count_GEN_genes(gene_info)

    promoter_cpgs = cpg_df[
        (cpg_df['chrom'] == chrom) &
        (cpg_df['start'] >= start) &
        (cpg_df['start'] <= end)
    ].copy()

    if len(promoter_cpgs) == 0:
        return pd.DataFrame()

    promoter_cpgs = promoter_cpgs.sort_values('start').reset_index(drop=True)
    n_cpg_sites = len(promoter_cpgs)
    cpg_positions = promoter_cpgs['start'].tolist()

    # Find individual promoters and overlap region for boundary checking (if provided)
    individual_promoters = []
    overlap_start, overlap_end = None, None
    
    if individual_promoters_df is not None and num_GEN_genes >= 2:
        individual_promoters = find_individual_promoters_in_merged(promoter_row, individual_promoters_df)
        overlap_start, overlap_end = find_overlap_region(individual_promoters)

    # Generate all possible windows within this promoter
    for window_size in range(1, n_cpg_sites + 1):
        for start_index in range(0, n_cpg_sites - window_size + 1, step_size):
            end_index = start_index + window_size

            start_pos = cpg_positions[start_index]
            end_pos = cpg_positions[end_index - 1]
            genomic_length = (end_pos + 1) - start_pos

            # Get CpG sites in this sub-window
            window_cpgs = promoter_cpgs.iloc[start_index:end_index]
            window_scores = window_cpgs['mean'].tolist()
            window_lamda_scores = window_cpgs['mean_lamda'].tolist()
            window_trust = window_cpgs['trustworthiness'].tolist()

            # Pass trustworthiness scores to stats calculation
            stats = window_stats_calculation(window_scores, window_lamda_scores, window_trust)
            cost_value = calculate_cost_function(stats, genomic_length, window_size, cost_function)
            # Create window record
            window_record = {
                'chrom': chrom,
                'start': start_pos,
                'end': end_pos,
                'gene_info': gene_info,
                'window_size_cpgs': window_size,
                'cost_value': cost_value,
            }

            # For multi-gene promoters with individual promoter data, check boundary
            if num_GEN_genes >= 2 and individual_promoters_df is not None:
                has_boundary = window_has_boundary_in_overlap(window_record, overlap_start, overlap_end)
                window_record['has_boundary_in_overlap'] = has_boundary

            results.append(window_record)

    return pd.DataFrame(results)

# Function to identify the best window
def find_best_window(df: pd.DataFrame) -> pd.Series:
    """Find the single best window (lowest cost)"""
    if len(df) == 0:
        return None

    # Filter out zero cost values and find minimum
    non_zero_df = df[df['cost_value'] != 0]
    if len(non_zero_df) == 0:
        return None

    min_index = non_zero_df['cost_value'].idxmin()
    return non_zero_df.loc[min_index]

# Function to find the optimal, non-overlapping windows for merged promoters with boundary constraints
def find_optimal_windows_for_merged_promoter(all_windows: pd.DataFrame, num_GEN_genes: int,
                                           individual_promoters: List[pd.Series]) -> List[pd.Series]:
    best_windows = []

    if len(all_windows) == 0 or num_GEN_genes < 2:
        return best_windows

    if len(individual_promoters) < 2:
        # Fall back to original method if we don't have individual promoter info
        return find_non_overlapping_best_windows(all_windows, num_GEN_genes)

    # Find overlap region for boundary checking
    overlap_start, overlap_end = find_overlap_region(individual_promoters)

    # Filter windows that have boundaries in the overlap region
    valid_windows = all_windows[all_windows['has_boundary_in_overlap'] == True]

    if len(valid_windows) == 0:
        # If no windows meet boundary criteria, use all windows
        valid_windows = all_windows

    # Find first best window from valid windows
    first_best = find_best_window(valid_windows)
    if first_best is None:
        return best_windows

    best_windows.append(first_best)

    # Find second non-overlapping window that also meets boundary criteria
    non_overlapping = valid_windows[
        (valid_windows['start'] >= first_best['end']) |
        (valid_windows['end'] <= first_best['start'])
    ]

    second_best = find_best_window(non_overlapping)
    if second_best is not None:
        best_windows.append(second_best)

    return best_windows

# Function to find the best non-overlapping windows if there are 2 promoters in an overlapping region
def find_non_overlapping_best_windows(all_windows: pd.DataFrame, num_GEN_genes: int) -> List[pd.Series]:
    """Find best non-overlapping windows based on number of GEN genes"""
    best_windows = []

    if len(all_windows) == 0:
        return best_windows

    # Always find the first best window
    first_best = find_best_window(all_windows)
    if first_best is not None:
        best_windows.append(first_best)

    # If we have 2 GEN genes, find a second non-overlapping best window
    if num_GEN_genes >= 2 and first_best is not None:
        # Filter out windows that overlap with the first best window
        non_overlapping = all_windows[
            (all_windows['start'] >= first_best['end']) |
            (all_windows['end'] <= first_best['start'])
        ]

        second_best = find_best_window(non_overlapping)
        if second_best is not None:
            best_windows.append(second_best)

    return best_windows

# Function to process one chromosome worth of promoters
def process_single_chromosome(cpg_df: pd.DataFrame, promoter_df: pd.DataFrame,
                            output_file: str, step_size: int,
                            target_chrom: str, cost_function: str,
                            individual_promoters_df: pd.DataFrame = None):
    if target_chrom is None:
        print("Error: No target chromosome specified")
        return

    print(f"Processing chromosome {target_chrom} with cost function {cost_function}...")

    # Filter for ONLY the target chromosome
    chrom_promoters = promoter_df[promoter_df['chrom'] == target_chrom]
    chrom_cpgs = cpg_df[cpg_df['chrom'] == target_chrom]
    
    if individual_promoters_df is not None:
        chrom_individual_promoters = individual_promoters_df[individual_promoters_df['chrom'] == target_chrom]
    else:
        chrom_individual_promoters = None

    print(f"Found {len(chrom_promoters)} promoters and {len(chrom_cpgs)} CpG sites for chromosome {target_chrom}")

    if len(chrom_promoters) == 0 or len(chrom_cpgs) == 0:
        print(f"No data found for chromosome {target_chrom}")
        return

    best_windows = []
    processed_count = 0

    for idx, promoter in chrom_promoters.iterrows():
        # Count GEN genes for this promoter
        num_GEN_genes = count_GEN_genes(promoter['gene_info'])

        # Analyze this promoter
        promoter_results = analyze_single_promoter(promoter, chrom_cpgs, cost_function, step_size, chrom_individual_promoters)

        if len(promoter_results) > 0:
            # Choose appropriate window selection strategy
            if num_GEN_genes >= 2 and chrom_individual_promoters is not None:
                # Get individual promoters for this region
                individual_promoters_for_region = find_individual_promoters_in_merged(promoter, chrom_individual_promoters)
                if len(individual_promoters_for_region) >= 2:
                    current_best_windows = find_optimal_windows_for_merged_promoter(
                        promoter_results, num_GEN_genes, individual_promoters_for_region)
                else:
                    current_best_windows = find_non_overlapping_best_windows(promoter_results, num_GEN_genes)
            else:
                # Use original method for single-gene promoters
                current_best_windows = find_non_overlapping_best_windows(promoter_results, num_GEN_genes)
            
            best_windows.extend(current_best_windows)

        processed_count += 1
        if processed_count % 100 == 0:
            print(f"Processed {processed_count}/{len(chrom_promoters)} promoters on chromosome {target_chrom}")

    # Save chromosome results as clean BED file
    if best_windows:
        # Create DataFrame with only essential BED columns
        bed_columns = ['chrom', 'start', 'end', 'gene_info']
        bed_df = pd.DataFrame(best_windows)[bed_columns]

        # Ensure output directory exists
        os.makedirs(os.path.dirname(output_file), exist_ok=True)

        bed_df.to_csv(output_file, sep='\t', index=False, header=False)
        print(f"Saved {len(bed_df)} best windows for chromosome {target_chrom} to {output_file}")
    else:
        print(f"No valid windows found for chromosome {target_chrom}")

# Function to define parameters
def main():
    parser = argparse.ArgumentParser(
        description='Sliding window analysis of CpG methylation in promoter regions with TSS-aware segmentation'
    )
    parser.add_argument('-c', '--cpg_file', required=True,
                       help='Input BED file with CpG sites (can be gzipped)')
    parser.add_argument('-p', '--promoter', required=True,
                        help='Input BED file with promoter windows')
    parser.add_argument('-i', '--individual_promoters', required=False,
                        help='Input BED file with individual promoter regions (TSS information)')
    parser.add_argument('-o', '--output', required=True,
                       help='Output prefix')
    parser.add_argument('-s', '--step', type=int, default=1,
                       help='Step size for sliding window (default: 1)')
    parser.add_argument('--chromosome', required=True,
                       help='Target chromosome to process (e.g., chr1, 1, X, etc.)')
    parser.add_argument('--cost_function', required=True,
                       choices=['cost_1', 'cost_2', 'cost_3', 'cost_4', 'cost_5', 'cost_6', 'cost_7', 'cost_8', 'cost_9', 'cost_10'],
                       help='Cost function to use for optimization')
    parser.add_argument('--total_samples', type=int, default=195,
                       help='Total number of samples for trustworthiness calculation (default: 195)')
    args = parser.parse_args()

    print("Reading CpG sites...")
    cpg_df = read_cpg_sites(args.cpg_file, args.total_samples)
    print(f"Loaded {len(cpg_df)} CpG sites")

    print("Reading promoter windows...")
    promoter_df = promoter_windows(args.promoter)
    print(f"Loaded {len(promoter_df)} promoter windows")

    # Read individual promoters if provided
    individual_promoters_df = None
    if args.individual_promoters:
        print("Reading individual promoter regions...")
        individual_promoters_df = read_individual_promoters(args.individual_promoters)
        print(f"Loaded {len(individual_promoters_df)} individual promoter regions")

    # Create output path: output_dir/cost_function/chromosome.bed
    output_file = f"{args.output}/{args.cost_function}/{args.chromosome}.bed"

    process_single_chromosome(
        cpg_df, promoter_df, output_file, args.step, args.chromosome, args.cost_function, individual_promoters_df
    )

    print(f"Analysis complete for chromosome {args.chromosome}")

if __name__ == "__main__":
    main()


