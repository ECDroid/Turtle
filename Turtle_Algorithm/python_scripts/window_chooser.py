#!/usr/bin/env python3

import pandas as pd
import numpy as np
import argparse
import sys
import os
import gzip
from typing import List, Dict

# Function to read TSV file input
def read_tsv_file(input_file: str) -> pd.DataFrame:
    df = pd.read_csv(input_file, sep = '\t')
    cost_columns = [col for col in df.columns if col.startswith('cost_')]

    if not cost_columns:
        print("Error: No columns that start with cost_")
        sys.exit(1)
    print(f"Cost columns: {cost_columns}")
    return df

def find_min_rows(df: pd.DataFrame) -> dict:
    min_rows = {}
    cost_columns = [col for col in df.columns if col.startswith('cost_')]

    for cost_col in cost_columns:
        non_zero_df = df[df[cost_col] != 0]
        if len(non_zero_df) == 0:
            print(f"Warning: All values for {cost_col} are 0. Skipping {cost_col}")
            continue
        min_value = non_zero_df[cost_col].min()
        min_rows[cost_col] = {
            'rows': non_zero_df[non_zero_df[cost_col] == min_value],
            'min_value': min_value
        }
    return min_rows

def write_output(min_rows: dict, output_pattern: str):
    output_files = []
    base_name = os.path.splitext(output_pattern)[0]
    extension = os.path.splitext(output_pattern)[1] or '.tsv'

    for cost_col, data in min_rows.items():
        output_file=f"{base_name}_{cost_col}{extension}"
        data['rows'].to_csv(output_file, sep='\t', index=False)
        output_files.append(output_file)
        window_sizes = data['rows']['window_size'].unique()

        print(f"{len(data['rows'])} row(s), window_sizes: {', '.join(map(str, window_sizes))} for {cost_col} = {data['min_value']}")
    return output_files


def main():
    parser = argparse.ArgumentParser(
        description='Find minimum cost sliding windows'
        )
    parser.add_argument('-i', '--input', required=True,
                        help='Input TSV file from brute_force_window.py output')
    parser.add_argument('-o', '--output', required=True,
                        help='Output TSV file naming convention')

    args = parser.parse_args()

    print("Reading files")
    cost_df = read_tsv_file(args.input)
    print("Finding min rows")
    dict_min = find_min_rows(cost_df)
    print("Writing output files")
    write_output(dict_min, args.output)

if __name__ == "__main__":
    main()





