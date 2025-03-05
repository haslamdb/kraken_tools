"""
This Python script covers reading Kraken/Bracken output tables,
merging them into a single dataframe, merging with sample metadata,
and preliminary data checks.

Please ensure that:
1) You have all your Kraken2/Bracken abundance files in a folder.
2) Each file is read and the species column is made consistent (no spaces, only dots, etc.).
3) You have a separate metadata CSV with sample IDs and relevant columns (e.g. location, week, antibiotic usage).

Usage:
python data_preprocessing.py --input_dir path_to_kraken_outputs --metadata path_to_metadata.csv --output raw_merged.csv

"""

import os
import sys
import argparse
import pandas as pd


def standardize_species_name(species_name: str) -> str:
    # Replace various characters with '.' or remove them,
    # akin to the R script approach.
    name = species_name.replace(' ', '.')
    name = name.replace('_', '.')
    name = name.replace('-', '.')
    name = name.replace('/', '.')
    name = name.replace('[', '')
    name = name.replace(']', '')
    name = name.replace('..', '.')  # repeated a few times in the R code
    return name


def load_kraken_abundance(filepath: str) -> pd.DataFrame:
    # Reads a single Kraken or Bracken output (TSV) that has columns:
    #   1) name of the taxon
    #   2) ...
    #   6) column with the estimated abundance or reads
    # We assume something like: 'Species', 'something', 'something', 'something', 'something', 'num_reads'
    # We'll read all columns, rename them, keep only species column and read count.

    df = pd.read_csv(filepath, sep='\t', header=0)

    # If we follow the typical Bracken style, the columns might be:
    # ["name", "taxonomy_id", "taxonomy_lvl", "kraken_assigned_reads",
    #  "added_reads", "new_est_reads", "fraction", "rank", "taxon_id"]
    # Or they might be slightly different. We'll assume the last col is abundance.

    # Let's just rename the columns generically
    df.columns = [f"col_{i}" for i in range(len(df.columns))]

    # We interpret col_0 as species name, col_5 (or last) as read count.
    # But let's just pick the first and 6th column if it exists.
    if len(df.columns) >= 6:
        df = df[["col_0", "col_5"]].copy()
        df.columns = ["Species", "Count"]
    else:
        raise ValueError(f"Not enough columns in file {filepath}")

    # Standardize species name.
    df["Species"] = df["Species"].apply(standardize_species_name)

    return df


def merge_kraken_outputs(input_dir: str) -> pd.DataFrame:
    # For each file in input_dir that ends with _species_abundance.txt or something,
    # read it in, store in a dictionary keyed by sample name.

    all_files = [f for f in os.listdir(input_dir) if f.endswith("_species_abundance.txt")]

    merged = None

    for f in all_files:
        # Sample name derived from the filename by removing the suffix
        sample_name = f.replace("_species_abundance.txt", "")
        filepath = os.path.join(input_dir, f)
        df_temp = load_kraken_abundance(filepath)
        # rename "Count" column to this sample
        df_temp.rename(columns={"Count": sample_name}, inplace=True)

        if merged is None:
            merged = df_temp
        else:
            merged = pd.merge(merged, df_temp, on="Species", how="outer")

    # Fill NA with 0
    if merged is not None:
        merged.fillna(0, inplace=True)
    else:
        raise ValueError("No abundance files found or empty data.")

    return merged


def main():
    parser = argparse.ArgumentParser(description="Combine Kraken/Bracken species abundance tables.")
    parser.add_argument("--input_dir", help="Path to directory with all the _species_abundance.txt files.")
    parser.add_argument("--metadata", help="Path to CSV with sample metadata.")
    parser.add_argument("--output", default="raw_merged.csv", help="Output CSV for merged table.")

    args = parser.parse_args()

    # 1) Merge all kraken/bracken outputs
    merged = merge_kraken_outputs(args.input_dir)

    # 2) Merged table is wide: row=species, columns=sample. We want to transpose to match the R approach.
    # Actually, let's keep it in wide form for now, or let's do both.

    # Write wide form
    merged.to_csv("wide_" + args.output, index=False)

    # Transpose: row=sample, columns=species
    merged_t = merged.set_index("Species").T.reset_index()
    merged_t.rename(columns={"index": "SampleID"}, inplace=True)
    merged_t.to_csv("transposed_" + args.output, index=False)

    # 3) Merge with sample metadata.
    if args.metadata:
        meta = pd.read_csv(args.metadata)
        # We'll assume meta has a column "SampleID"
        # Then we do a left or inner join on SampleID.
        # If the user wants to do it differently, adapt.

        merged_meta = pd.merge(meta, merged_t, on="SampleID", how="inner")
        merged_meta.to_csv(args.output, index=False)

    print("Preprocessing complete.")

if __name__ == "__main__":
    main()
