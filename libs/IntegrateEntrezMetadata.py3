"""
The purpose of this script is to take the results objects from more than one query to ObtainEntrezMetadata.py3 and integrate them, anding completeness to as many records as possible.
for instance, calls to sra, biosample, bioproject, and geo with regard to gene expression data in humans produces non-redundant datasets.

Integration of this data is aimed at 3 main goals:
  1) leverage non-redundant data streams to create more comprehensive metadata possible using 1 alone.
  2) streamline the process of manual study identification and curation / allow more rapid prioritization of candidate studies.
  3) identify / draw attention to studies that may not have been identified without an agnostic scan.
"""

import pandas as pd
import sys


""" Step 1: Identify candidate primary keys """
def identify_candidate_primary_keys(df):
    num_rows = len(df)
    column_names = df.columns.tolist()
    percentage_data = {}
    primary_key_candidates = []

    for column in column_names:
        not_null = df[column].notnull()
        not_empty = df[column].replace('', pd.NA).notna()
        valid_data = not_null & not_empty
        total_valid = valid_data.sum()
        empty_strings = (not_null & ~not_empty).sum()
        nans = num_rows - not_null.sum()

        percentage_data[column] = {
            'Valid data %': 100 * total_valid / num_rows,
            'Empty strings %': 100 * empty_strings / num_rows,
            'NaNs %': 100 * nans / num_rows
        }

        if nans == 0 and df[column].is_unique:
            primary_key_candidates.append(column)

    return primary_key_candidates, percentage_data

""" Step 2: Match the contents of all candidate primary keys. Return the amount of matching content per candidate pair
    Ideally, there should be 100% mapping of df with fewer rows onto the df with more rows if the keys are ideal. """
def match_key_candidate_contents(df1, primary_keys1, df2, primary_keys2):
    match_percentages = {}

    for key1 in primary_keys1:
        for key2 in primary_keys2:
            if key1 in df1.columns and key2 in df2.columns:
                unique_values_df1 = set(df1[key1].dropna().unique())
                unique_values_df2 = set(df2[key2].dropna().unique())
                common_values = unique_values_df1.intersection(unique_values_df2)
                total_values_df2 = len(unique_values_df2)

                if total_values_df2 > 0:
                    match_percentage = 100 * len(common_values) / total_values_df2
                else:
                    match_percentage = 0

                match_percentages[(key1, key2)] = match_percentage

    return match_percentages

""" Run Steps 1 and 2 """
def main(file1, file2):
    # Load the data from the input files into pandas DataFrames
    df1 = pd.read_csv(file1)
    df2 = pd.read_csv(file2)

    # Identify primary key candidates and percentage data for both DataFrames
    primary_keys1, percentages1 = identify_candidate_primary_keys(df1)
    primary_keys2, percentages2 = identify_candidate_primary_keys(df2)

    # Match the contents of the potential primary keys between both DataFrames
    matches = match_key_candidate_contents(df1, primary_keys1, df2, primary_keys2)

    # Print the match percentages
    print("Match Percentages:")
    for keys, percentage in matches.items():
        print(f"{keys}: {percentage}%")

""" File names, which are also input args """
inDir='/data/user/vlaufer/DeathStar/Data/merge'
file1=f'{inDir}/biosamples.proc.tsv'
file2=f'{inDir}/bioprojects.proc.tsv'
outfile=f'{inDir}/merge.proc.tsv'

""" Examples of ways to call the main function: """
if len(sys.argv) > 2:
	main(sys.argv[1], sys.argv[2], sys.argv[3])
else:
	main(file1, file2, outfile)


