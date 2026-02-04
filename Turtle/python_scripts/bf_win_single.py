#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import sys
import gzip
from typing import List, Dict

# Global variable for total samples if needed - you might want to pass this as an argument
total_samples = None # You should set this from command line or calculate it

def read_cpg_sites(input_file: str) -> pd.DataFrame:
    """Read CpG sites from input file with the specified format."""
    try:
        df = pd.read_csv(input_file, sep='\t', header=None,
                        usecols=[0, 1, 2, 3, 4, 5, 6, 9],
                        names=['chrom', 'start', 'end', 'mean', 'median', 'std', 'variance', 'no_of_samples'])

        # NEW: Keep original mean and create transformed versions
        df['original_mean'] = df['mean'].copy()  # Store original 0-100 values
        
        # Create per-site lambda (distance from nearest extreme)
        df['mean_lambda'] = df['mean'].apply(lambda x: 100 - x if x > 50 else x)
        
        # For compatibility, keep 'mean' as the lambda value (as before)
        df['mean'] = df['mean_lambda'].copy()
        
        # Fix: You need to define total_samples or pass it as argument
        if total_samples is None:
            total_samples_val = df['no_of_samples'].max()
        else:
            total_samples_val = total_samples

        df['trustworthiness'] = (df['no_of_samples'] / total_samples_val) * (1 / (1 + df['variance']))

        df = df.sort_values(['chrom', 'start']).reset_index(drop=True)
        return df
    except Exception as e:
        print(f"Error reading input file: {e}")
        sys.exit(1)

def window_stats_calculation(scores: List[float], trustworthiness_scores: List[float] = None, 
                           original_scores: List[float] = None) -> Dict[str, float]:
    """Calculate window statistics including weighted stats if trustworthiness is provided."""
    if len(scores) == 0:
        return {
            'std': 0.0, 'variance': 0.0, 'mean': 0.0, 'median': 0.0,
            'weighted_std': 0.0, 'weighted_variance': 0.0, 'weighted_mean': 0.0,
            'lamda_window_based': 0.0  # NEW
        }

    scores_array = np.array(scores)

    # NEW: Calculate window-based lambda
    lamda_window_based = 0.0
    if original_scores is not None and len(original_scores) > 0:
        window_mean = np.mean(original_scores)
        if window_mean > 50:
            lamda_window_based = (100 - window_mean)  # 0-50 range
        else:
            lamda_window_based = window_mean  # 0-50 range

    # Calculate basic statistics
    result = {
        'mean': float(np.mean(scores_array)),
        'std': float(np.std(scores_array)),
        'variance': float(np.var(scores_array)),
        'median': float(np.median(scores_array)),
        'lamda_window_based': float(lamda_window_based)  # NEW
    }

    # Calculate weighted statistics if trustworthiness scores are provided
    if trustworthiness_scores and len(trustworthiness_scores) > 0:
        weights_array = np.array(trustworthiness_scores)
        weighted_mean = np.average(scores_array, weights=weights_array)
        weighted_variance = np.average((scores_array - weighted_mean)**2, weights=weights_array)
        weighted_std = np.sqrt(weighted_variance)

        result.update({
            'weighted_std': float(weighted_std),
            'weighted_variance': float(weighted_variance),
            'weighted_mean': float(weighted_mean)
        })
    else:
        result.update({
            'weighted_std': float(np.std(scores_array)),
            'weighted_variance': float(np.var(scores_array)),
            'weighted_mean': float(np.mean(scores_array))
        })

    return result

def calculate_cost_function(stats: Dict, genomic_length: int, cost_function: str, window_size: int = 1) -> float:
    """Calculate specific cost function based on user choice."""
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
        # NEW: Your proposed window-based lambda with weighted variance
        return (stats['lamda_window_based'] + stats['weighted_variance']) / window_size
    elif cost_function == 'cost_7':
        return (stats['mean'] + stats['weighted_variance'])
    elif cost_function == 'cost_8':
        return (stats['mean'] + stats['weighted_variance']) / window_size
    elif cost_function == 'cost_9':
        # NEW: Scaled version - lambda 0-50, variance <20, so scale lambda down
        scaled_lambda = stats['lamda_window_based'] / 2.5  # Now 0-20 range
        return (scaled_lambda + stats['weighted_variance']) / window_size
    elif cost_function == 'cost_10':
        # NEW: No window size division
        scaled_lambda = stats['lamda_window_based'] / 2.5
        return scaled_lambda + stats['weighted_variance']
    elif cost_function == 'cost_11':
        # NEW: Square the variance to penalize inconsistency
        scaled_lambda = stats['lamda_window_based'] / 2.5
        return (scaled_lambda + stats['weighted_variance']**2) / window_size
    else:
        raise ValueError(f"Unknown cost function: {cost_function}")

def sliding_window_analysis(cpg_df: pd.DataFrame, cost_function: str, step_size: int = 1) -> pd.DataFrame:
    """Perform sliding window analysis and calculate cost values."""
    results = []
    n_cpg_sites = len(cpg_df)

    if n_cpg_sites < 1:
        print("Warning: No CpG sites for analysis")
        return pd.DataFrame()

    chrom = cpg_df['chrom'].iloc[0]
    print(f"Analyzing {n_cpg_sites} CpG sites on {chrom}...")

    cpg_positions = sorted(cpg_df['start'].unique())

    for window_size in range(3, n_cpg_sites + 1):
        for start_index in range(0, n_cpg_sites - window_size + 1, step_size):
            end_index = start_index + window_size

            start_pos = cpg_positions[start_index]
            end_pos = cpg_positions[end_index - 1]
            genomic_length = (end_pos + 1) - start_pos

            window_df = cpg_df[
                (cpg_df['start'] >= start_pos) &
                (cpg_df['start'] <= end_pos)
            ]

            if len(window_df) > 0:
                # Use lambda values (per-site transformed)
                window_scores = window_df['mean'].tolist()  # These are lambda values
                window_trust = window_df['trustworthiness'].tolist()
                
                # NEW: Get original methylation scores for window-based lambda
                window_original_scores = window_df['original_mean'].tolist()

                stats = window_stats_calculation(window_scores, window_trust, window_original_scores)
                window_size_cpg = len(window_df)
                cost_value = calculate_cost_function(stats, genomic_length, cost_function, window_size_cpg)

                results.append({
                    'chrom': chrom,
                    'start': start_pos,
                    'end': end_pos,
                    'window_size': window_size,
                    'std_dev': stats['std'],
                    'variance': stats['variance'],
                    'mean': stats['mean'],
                    'median': stats['median'],
                    'weighted_std': stats['weighted_std'],
                    'weighted_variance': stats['weighted_variance'],
                    'weighted_mean': stats['weighted_mean'],
                    'lamda_window_based': stats['lamda_window_based'],  # NEW
                    'genomic_length': genomic_length,
                    'cost_value': cost_value,
                })

    return pd.DataFrame(results)

def transform_value(value: float, transform: str = None) -> float:
    """Apply transformation to value for better visualization."""
    if transform == 'log':
        return np.log1p(value) # log(1 + x) to handle zeros
    elif transform == 'log10':
        return np.log10(value + 1) # log10(1 + x)
    elif transform == 'sqrt':
        return np.sqrt(value) # square root
    elif transform == 'inverse':
        return 1 / (value + 1e-10) # inverse (add small epsilon to avoid division by zero)
    elif transform == 'neg_log':
        return -np.log(value + 1e-10) # negative log for emphasizing low values
    else:
        return value # no transformation

def write_tsv_output(results_df: pd.DataFrame, output_file: str, cost_function: str):
    """Write results to TSV file."""
    if len(results_df) == 0:
        print("No results to write")
        return

    # Create a copy to avoid modifying original
    output_df = results_df.copy()

    # Format all numeric columns to 4 decimal places
    numeric_cols = output_df.select_dtypes(include=[np.number]).columns
    for col in numeric_cols:
        output_df[col] = output_df[col].round(4)

    # Add a column indicating which cost function was used
    output_df['cost_function'] = cost_function

    # Write to TSV
    output_df.to_csv(output_file, sep='\t', index=False)
    print(f"Results written to {output_file}")

def position_matrix(results_df: pd.DataFrame, matrix_type: str, cost_function: str = None) -> pd.DataFrame:
    """Create matrix of results based on selected statistic or cost value."""
    # If cost_function is specified, use cost_value, otherwise use the specified stat
    if cost_function:
        stat_column = 'cost_value'
        matrix_type = cost_function # Use cost function name for labeling
    else:
        # Map matrix_type to column names
        type_to_column = {
            'mean': 'mean',
            'variance': 'variance',
            'std': 'std_dev',
            'median': 'median',
            'weighted_mean': 'weighted_mean',
            'weighted_std': 'weighted_std',
            'weighted_variance': 'weighted_variance',
            'lamda_window_based': 'lamda_window_based'  # NEW
        }

        if matrix_type in type_to_column:
            stat_column = type_to_column[matrix_type]
        else:
            raise ValueError(f"Unknown matrix type: {matrix_type}")

    # Extract unique start positions
    all_starts = sorted(results_df['start'].unique())

    # Create an empty matrix with NaN values
    matrix_data = pd.DataFrame(
        index=all_starts,
        columns=all_starts,
        dtype=float
    )

    # Fill matrix with values
    for _, row in results_df.iterrows():
        i = row['start']
        j = row['end']
        value = row[stat_column]
        matrix_data.loc[i, j] = value

    matrix_data = matrix_data.iloc[::-1] # Reverse y-axis
    return matrix_data

def heatmap(matrix_df: pd.DataFrame, output_png: str, matrix_type: str, cost_function: str = None, transform: str = None):
    """Generate heatmap from matrix."""
    plt.figure(figsize=(12, 10))
    mask = matrix_df.isna()

    # Apply transformation if specified
    if transform:
        matrix_df_transformed = matrix_df.copy()
        matrix_df_transformed = matrix_df_transformed.applymap(lambda x: transform_value(x, transform) if not pd.isna(x) else x)
    else:
        matrix_df_transformed = matrix_df

    # Set appropriate label
    if cost_function:
        label = f'Cost Value ({cost_function})'
        if transform:
            label = f'{transform.capitalize()}({label})'
    else:
        label = matrix_type.upper()
        if transform:
            label = f'{transform.capitalize()}({label})'

    if cost_function and not matrix_df.empty:
        min_val = matrix_df.min().min()
        if not pd.isna(min_val):
            # Create a mask for the minimum value(s)
            min_mask = matrix_df == min_val
            if min_mask.any().any():
                # Create a custom colormap or annotation
                plt.text(0.02, 0.98, f'Min cost: {min_val:.6f}',
                        transform=plt.gca().transAxes,
                        fontsize=10, verticalalignment='top',
                        bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.5))

    sns.heatmap(matrix_df_transformed,
                mask=mask,
                cmap='viridis',
                cbar_kws={'label': label},
                square=True)

    title = f'Genomic Position Matrix - {label}'
    plt.title(title)
    plt.xlabel('End Position (j)')
    plt.ylabel('Start Position (i)')
    plt.tight_layout()
    plt.savefig(output_png, dpi=300, bbox_inches='tight')
    plt.close()

def main():
    parser = argparse.ArgumentParser(
        description='Sliding window analysis of CpG methylation in promoter regions'
    )
    parser.add_argument('-i', '--input', required=True,
                       help='Input BED file with CpG sites (can be gzipped)')
    parser.add_argument('-o', '--output', required=True,
                       help='Output TSV file')
    parser.add_argument('-c', '--cost_function', default='cost_1',
                       choices=['cost_1', 'cost_2', 'cost_3', 'cost_4',
                               'cost_5', 'cost_6', 'cost_7', 'cost_8',
                               'cost_9', 'cost_10', 'cost_11'],
                       help='Cost function to use (default: cost_1)')
    parser.add_argument('-s', '--step', type=int, default=1,
                       help='Step size for sliding window (default: 1)')
    parser.add_argument('-t', '--total_samples', type=int, default=195,
                       help='Total number of samples for trustworthiness calculation ')
    parser.add_argument('-m', '--matrix',
                       type=str,
                       nargs='?',
                       const='mean',
                       choices=['mean', 'variance', 'std', 'median',
                               'weighted_mean', 'weighted_std', 'weighted_variance',
                               'lamda_window_based', 'cost', 'all'],
                       help='Creates a matrix PNG output - choose statistic or "cost" for cost value (default: mean)')
    parser.add_argument('--transform',
                       choices=['log', 'log10', 'sqrt', 'inverse', 'neg_log', 'none'],
                       default='none',
                       help='Transformation to apply to cost values for visualization (default: none)')

    parser.add_argument('--highlight_min',
                       action='store_true',
                       help='Highlight the minimum cost value in the heatmap')
    args = parser.parse_args()

    # Set global total_samples if provided
    if args.total_samples:
        global total_samples
        total_samples = args.total_samples

    print("Reading CpG sites...")
    cpg_df = read_cpg_sites(args.input)
    print(f"Loaded {len(cpg_df)} CpG sites")

    print(f"Performing sliding window analysis with cost function: {args.cost_function}...")
    results_df = sliding_window_analysis(cpg_df, args.cost_function, args.step)
    print(f"Generated {len(results_df)} window analyses")

    print("Writing results to TSV file...")
    write_tsv_output(results_df, args.output, args.cost_function)

    if args.matrix is not None:
        matrix_types = []

        if args.matrix == 'all':
            # Include all statistic types plus cost value
            matrix_types = ['mean', 'variance', 'std', 'median',
                           'weighted_mean', 'weighted_std', 'weighted_variance',
                           'lamda_window_based']  # ADDED
            # Add cost matrix separately
            cost_matrix = True
        elif args.matrix == 'cost':
            matrix_types = [] # Will handle cost separately
            cost_matrix = True
        else:
            matrix_types = [args.matrix]
            cost_matrix = False

        base_name = args.output.replace('.tsv', '')

        # Generate matrices for regular statistics
        for matrix_type in matrix_types:
            print(f"Creating {matrix_type} matrix")
            matrix_df = position_matrix(results_df, matrix_type, None)
            output_png = f"{base_name}_{matrix_type}.png"
            output_tsv = f"{base_name}_{matrix_type}_matrix.tsv"

            heatmap(matrix_df, output_png, matrix_type, None)
            matrix_df.to_csv(output_tsv, sep='\t', na_rep='NA')
            print(f"Matrix and heatmap data saved as {output_png}")

        # Generate cost matrix if requested
        if args.matrix == 'cost' or (args.matrix == 'all' and cost_matrix):
            print(f"Creating cost matrix for {args.cost_function}")
            matrix_df = position_matrix(results_df, None, args.cost_function)
            output_png = f"{base_name}_{args.cost_function}_cost.png"
            if args.transform != 'none':
                output_png = f"{base_name}_{args.cost_function}_cost_{args.transform}.png"
            output_tsv = f"{base_name}_{args.cost_function}_cost_matrix.tsv"

            heatmap(matrix_df, output_png, args.cost_function, args.cost_function, args.transform)
            matrix_df.to_csv(output_tsv, sep='\t', na_rep='NA')
            print(f"Cost matrix and heatmap data saved as {output_png}")

    print("All processes completed!")

if __name__ == "__main__":
    main()


