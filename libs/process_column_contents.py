from collections import Counter
import pandas as pd

def remove_sparse_columns(df, threshold=0.95):
    """ Removes columns that have a percentage of missingness above a specified threshold. """
    original_columns = set(df.columns)
    min_count = int((1 - threshold) * len(df))
    df = df.dropna(axis=1, thresh=min_count)
    removed_columns = original_columns - set(df.columns)
    return df, list(removed_columns)

def remove_identical_columns(df, threshold=0.9):
    """ Removes columns that have a percentage of identity above a specified threshold. """
    columns = df.columns.tolist()
    redundant_cols = []
    for i in range(len(columns)):
        for j in range(i + 1, len(columns)):
            if (df[columns[i]] == df[columns[j]]).mean() > threshold:
                redundant_cols.append(columns[j])
    df = df.drop(columns=set(redundant_cols))
    return df, redundant_cols

def remove_substrings_from_colnames(df, threshold_percentage=0.20):
    """ Splits column names and counts occurrences of each part. """
    original_col_names = df.columns.tolist()
    df.columns = [col.lstrip('_') for col in df.columns]
    col_names = df.columns.tolist()
    parts_list = [name.split('_') for name in col_names]
    flat_list = [part for sublist in parts_list for part in sublist]
    part_counts = Counter(flat_list)
    threshold = len(col_names) * threshold_percentage
    common_leading_parts = {part for part, count in part_counts.items() if count > threshold and any(name.startswith(part + '_') for name in col_names)}
    new_col_names = []
    for name in col_names:
        for part in common_leading_parts:
            if name.startswith(part + '_'):
                name = name.replace(part + '_', '', 1)
                break
        new_col_names.append(name)
    removed_col_names = [original for original in original_col_names if original not in new_col_names]
    df.columns = new_col_names
    return df, removed_col_names

def main():
    dbIn = '/data/user/vlaufer/DeathStar/Data/proc/sra/sra_db_pre_proc.tsv'
    dbOut = '/data/user/vlaufer/DeathStar/Data/proc/sra/sra_db_processed.tsv'
    try:
        sra_db = pd.read_csv(dbIn, sep='\t')
        print('Data read into DF.')

        sraDbSlim, removed_columns = remove_sparse_columns(sra_db)
        print(f'Sparse columns removed: {len(removed_columns)}', removed_columns)
        print(sraDbSlim.head())

        sraDbSlim, redundant_cols = remove_identical_columns(sraDbSlim)
        print(f'Highly similar columns removed: {len(redundant_cols)}', redundant_cols)
        print(sraDbSlim.head())

        sraDbSlim.to_csv(dbOut, index=False, sep='\t')
        print(f'Processed data saved to {dbOut}.')
    except Exception as e:
        print(f"An error occurred: {e}")

main()

