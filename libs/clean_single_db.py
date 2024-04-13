import pandas as pd
from collections import Counter
import sys

def remove_sparse_columns(df):
	""" Step 1.1 - Removes columns that have more than 95% missing values """
	threshold = 0.95
	return df.dropna(thresh=(1 - threshold) * len(df), axis=1)

def standardize_data(df):
	""" Step 1.2 - Creates a copy of the DataFrame and standardizes all text to lower case and strips whitespace """
	for col in df.columns:
		df[col] = df[col].apply(lambda x: x.lower().strip() if isinstance(x, str) else x)
	return df

def remove_empty_rows(df):
    """ drop rows that are entirely empty """ 
    df = df.dropna(axis=1, how='all')
    return(df)

def remove_sparse_columns(df, threshold=0.95):
    """ Removes columns that have a percentage of missingness above a specified threshold. """
    original_columns = set(df.columns)
    min_count = int((1 - threshold) * len(df))
    df = df.dropna(axis=1, thresh=min_count)
    removed_columns = original_columns - set(df.columns)
    return df, list(removed_columns)

def identify_primary_key_candidates(df):
	""" Step 1.4 - Removes columns from both the original and its copy that are completely identical """
	candidates = []
	for col in df.columns:
		if df[col].is_unique and df[col].notna().all():
			candidates.append(col)
	return candidates

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

def remove_identical_columns(df, df_copy):
	"""Identifies columns that have a one-to-one relationship, useful for deduplication and data integrity."""
	identical_cols = set()  # Use a set to prevent duplicates
	columns = df.columns
	identity_check = {}
	for i in range(len(columns)):
		for j in range(i + 1, len(columns)):
			col1 = columns[i]
			col2 = columns[j]
			# Check if the columns are identical
			identical = (df[col1] == df[col2]).all()
			# Store the result in the dictionary
			identity_check[(col1, col2)] = identical
			# If columns are identical, add them to the set of columns to drop
			if identical:
				identical_cols.add(col2)  # Typically add the second column, or choose strategy

	# Drop the identical columns from both dataframes
	df.drop(columns=list(identical_cols), inplace=True)
	df_copy.drop(columns=list(identical_cols), inplace=True)
	return df, df_copy


def identify_bijective_relationships(df_copy):
	""" Step 1.5 - Identify columns that map onto one another in a one-to-one fashion, even if the strings differ."""
	num_cols = len(df_copy.columns)
	cols_to_drop = []
	mapping_tables = []

	# Precompute unique value counts and non-null unique value counts for each column
	unique_counts = {col: df_copy[col].dropna().unique().size for col in df_copy.columns}
	unique_nunique_counts = {col: df_copy[col].dropna().nunique() for col in df_copy.columns}

	for i in range(num_cols):
		for j in range(i + 1, num_cols):
			col_i_name = df_copy.columns[i]
			col_j_name = df_copy.columns[j]
			col_i = df_copy.iloc[:, i].dropna()
			col_j = df_copy.iloc[:, j].dropna()

			if (unique_counts[col_i_name] == unique_counts[col_j_name] ==
				unique_nunique_counts[col_i_name] == unique_nunique_counts[col_j_name] and
				col_i.isin(col_j).all() and col_j.isin(col_i).all()):
				cols_to_drop.append(col_j_name)
				mapping_table = pd.DataFrame({
					col_i_name: col_i.unique(),
					col_j_name: pd.Series(col_i.unique()).map(pd.Series(col_j.unique(), index=col_i.unique())).values
				})
				mapping_tables.append(mapping_table)
	return cols_to_drop, mapping_tables

def surjective_cols_to_lookup_table(df, df_copy, colToKeep, colToDrop):
	""" Step 1.6 - Following identification of a surjective relationship, this function will write the one-to-one 
	correspondence to a separate DataFrame and update the original DataFrame by dropping one column."""
	lookup_df = df[[colToKeep, colToDrop]]       # Create a new DF with only the two specified columns
	df_copy = df_copy.dropna().drop_duplicates() # Remove NaN values and drop duplicates from the new DataFrame
	df_copy = df_copy.drop(columns=[colToDrop])  # Drop the specified column to drop from the original DataFrame
	df = df.drop(columns=[colToDrop])            # Now finally drop the specified columns to drop from the **original** DataFrame
	return df, df_copy, lookup_df                # Return both the modified original DataFrame and the new DataFrame

def manage_bijective_and_lookup(df, df_copy):
	""" Steps 1.5 and 1.6 - Manage bijective relationships and create lookup tables """
	dropped_columns, mapping_tables = identify_bijective_relationships(df_copy)
	for table in mapping_tables:
		colToKeep, colToDrop = table.columns[0], table.columns[1]
		df, df_copy, new_lookup_table = surjective_cols_to_lookup_table(df, df_copy, colToKeep, colToDrop)
		print(f"Lookup table created for {colToKeep} and {colToDrop}")
		print(new_lookup_table.head())  # Display a preview of the lookup table
	return df

def write_outfiles(df_Out, outfile):
	outfileXlsx=outfile + '.xlsx'
	df_Out, column_mapping_df = remove_substrings_from_colnames(df_Out, 0.20)
	try:
		with pd.ExcelWriter(outfileXlsx, engine='openpyxl') as writer:
			df_Out.to_excel(writer, sheet_name='MetaData', index=False)
			column_mapping_df.to_excel(writer, sheet_name='ColumnMapping', index=False)
#
	except Exception as e:
		print(f"Failed to write using openpyxl due to {e}. Trying xlsxwriter...")
		try:
			with pd.ExcelWriter(outfileXlsx, engine='xlsxwriter') as writer:
				df_Out.to_excel(writer, sheet_name='MetaData', index=False)
				column_mapping_df.to_excel(writer, sheet_name='ColumnMapping', index=False)
#
		except Exception as e:
			print(f"Failed to write using xlsxwriter as well due to {e}. Writing .tsv instead.")   # Additional error handling or fallback logic can be implemented here
			outfileTsv=outfile + '.tsv'
			df_Out.to_csv(outfileTsv, sep='\t', index=False) # Write output to tsv
			outfileColTsv=outfile + '.col_mapping.tsv'
			column_mapping_df.to_csv(outfileColTsv, sep='\t', index=False) # Write output to tsv

def main(infile, outfile_stem):
	""" Step 0: Main function to load data, apply cleaning transformations, and print out the final DataFrame """
	try:
		df = pd.read_csv(infile, delimiter='\t', dtype=str, encoding='utf-8')
	except UnicodeDecodeError:
		try:
			df = pd.read_csv(infile, delimiter='\t', dtype=str, encoding='ISO-8859-1')
		except UnicodeDecodeError:
			df = pd.read_csv(infile, delimiter='\t', dtype=str)  # Fallback to default encoding

	df = remove_sparse_columns(df)  # Step 1.1
	df = standardize_data(df)  # Step 1.2
	df_copy=df.copy()
	df, df_copy = remove_identical_columns(df, df_copy)  # Step 1.3
	primary_keys = identify_primary_key_candidates(df)  # Step 1.4
	print("Primary Key Candidates:", primary_keys)
	final_df = manage_bijective_and_lookup(df, df_copy)  # Step 1.5 and 1.6
	write_outfiles(final_df, outfile_stem)

	return final_df

if __name__ == "__main__":
	import sys
	if len(sys.argv) > 1:
		main(sys.argv[1])  # Pass the input file as an argument
	else:
		main('../path/to/input/Biosample_Metadata_Formatted.txt', '../path/to/output/biosample.0.95.proc')
"""
requires integration into main() above. formerly from dedicated column name handling script.
def main():
    # Read from stdin, assume each input line is a single string of tab-delimited fields
    input_lines = [line.strip() for line in sys.stdin if line.strip()]
    fields = [item for line in input_lines for item in line.split('\t')]

    # Process the fields to refine and make them unique
    refined_fields = refine_segments(fields)

    # Output the final unique fields
    for field in refined_fields:
        print(field)

if __name__ == "__main__":
    main()

"""
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
"""
