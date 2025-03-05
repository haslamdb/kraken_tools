"""
This script demonstrates how to:
1) Filter taxa present in <20% of samples.
2) Filter taxa that make up <0.01% of total assigned reads.
3) Remove samples that have very low total reads.
4) Rarefy ("rrarefy") the table to a minimum count.

We rely on the merged table that has row=sample, columns=metadata + all species.
We'll define some utility functions.

Usage:
python normalization_and_filtering.py --input merged.csv --output filtered.csv

"""

import argparse
import pandas as pd
import numpy as np
import random


def prevalence_filter(df_counts: pd.DataFrame, min_prevalence_fraction=0.2) -> pd.DataFrame:
    """
    Remove species (columns) that are present in fewer than 20% of the samples.
    df_counts: row=samples, columns=species
    """
    n_samples = df_counts.shape[0]
    min_nonzero = int(np.floor(n_samples * min_prevalence_fraction))

    keep_columns = []
    for col in df_counts.columns:
        # number of nonzero in this column
        nonzero_count = (df_counts[col] > 0).sum()
        if nonzero_count >= min_nonzero:
            keep_columns.append(col)
    return df_counts[keep_columns]


def abundance_filter(df_counts: pd.DataFrame, min_fraction=1e-4) -> pd.DataFrame:
    """
    Remove species that make up <0.01% of total reads.
    """
    total_assigned = df_counts.sum().sum()
    # compute sum of each species across all samples
    col_sums = df_counts.sum(axis=0)
    keep_columns = col_sums[col_sums / total_assigned >= min_fraction].index
    return df_counts[keep_columns]


def rrarefy(df_counts: pd.DataFrame, min_count=None, random_seed=42) -> pd.DataFrame:
    """
    Rarefies each sample to the same number of reads (min_count) by random subsampling.
    If min_count not provided, we find the smallest sample sum.
    This replicates the R vegan::rrarefy approach.
    """
    random.seed(random_seed)
    if min_count is None:
        sample_sums = df_counts.sum(axis=1)
        min_count = sample_sums.min()

    # We'll do it row by row.

    rarefied_data = []
    for i in range(df_counts.shape[0]):
        row = df_counts.iloc[i].values
        total = row.sum()
        if total <= min_count:
            # no changes, the row has less or equal to min_count.
            rarefied_data.append(row)
        else:
            # sample from the distribution of reads
            # approach: create a list of indices repeated by read counts, then pick min_count.
            # Then sum up how many times each index was picked.
            indices = np.repeat(np.arange(len(row)), row)
            chosen = np.random.choice(indices, size=min_count, replace=False)
            # count how many times each index was chosen
            new_counts = np.bincount(chosen, minlength=len(row))
            rarefied_data.append(new_counts)
    rarefied_data = np.array(rarefied_data)
    return pd.DataFrame(rarefied_data, columns=df_counts.columns, index=df_counts.index)


def main():
    parser = argparse.ArgumentParser(description="Filter and rarefy.")
    parser.add_argument("--input", required=True, help="Merged table with row=sample, columns=metadata + species.")
    parser.add_argument("--output", default="filtered.csv", help="Output CSV.")
    parser.add_argument("--metadata_cols", type=int, default=5,
                        help="Number of left-most columns that are metadata (e.g., sample ID, location, etc.)")
    parser.add_argument("--min_prevalence", type=float, default=0.2,
                        help="Minimum fraction of samples in which a species must appear.")
    parser.add_argument("--min_fraction", type=float, default=1e-4,
                        help="Minimum fraction of total reads for a species.")
    parser.add_argument("--min_reads", type=int, default=100000,
                        help="Min reads for a sample.")

    args = parser.parse_args()

    df = pd.read_csv(args.input)

    # We assume the first 'metadata_cols' columns are metadata.
    meta = df.iloc[:, :args.metadata_cols].copy()
    counts = df.iloc[:, args.metadata_cols:].copy()

    # 1) remove samples with < min_reads
    sample_sums = counts.sum(axis=1)
    good_indices = sample_sums[sample_sums > args.min_reads].index
    meta = meta.loc[good_indices]
    counts = counts.loc[good_indices]

    # 2) filter species by prevalence
    counts = prevalence_filter(counts, min_prevalence_fraction=args.min_prevalence)

    # 3) filter species by total fraction
    counts = abundance_filter(counts, min_fraction=args.min_fraction)

    # 4) rarefy
    counts_rarefied = rrarefy(counts)

    # Combine back with metadata
    out = pd.concat([meta, counts_rarefied], axis=1)
    out.to_csv(args.output, index=False)

    print("Normalization & filtering complete.")

if __name__ == "__main__":
    main()
