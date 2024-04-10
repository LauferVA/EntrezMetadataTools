import pandas as pd

def Duplicate_And_PreProcess_Db(df):
	"""Preprocesses the DataFrame by stripping whitespace and converting to lowercase."""
	df_processed = df.copy()
	for col in df_processed.columns:
		df_processed[col] = df_processed[col].str.strip().str.lower()
	return df_processed

def Identify_and_Remove_Duplicate_Cols(df):
	"""Identifies completely identical columns in a DataFrame and removes them."""
	cols_to_drop = []
	num_cols = len(df.columns)
	for i in range(num_cols):
		for j in range(i+1, num_cols):
			if df.iloc[:, i].equals(df.iloc[:, j]):
				cols_to_drop.append(df.columns[j])
	df_cleaned = df.drop(columns=cols_to_drop) # Drop identified columns
	return df_cleaned, cols_to_drop

def Identify_Primary_Key_Cols(df):
  """ Identify the primary key by proving a column has no empty fields and is bijective unto itself. We do this simply by: 
  1. 
  2. Showing the number of unique values in the column considered is equivalent to the number of rows of the df 
  - 
  """

def Surjective_Cols_To_LookUp(df, colToKeep, colToDrop):
	""" Following the identification of a surjective relationship, this function will write the one-to-one correspondence to a separate Df """
	new_df = df[[colToKeep, colToDrop]]			# Create a new DF with only the two specified columns
	new_df = new_df.fillna('')					# Convert NaN values in the new DataFrame to 0-length strings
	new_df = new_df.drop_duplicates()			# Drop duplicates from the new DataFrame
	modified_df = df.drop(columns=[colToDrop])	# Drop the specified column to drop from the original DataFrame
	return modified_df, new_df					# Return both the modified original DataFrame and the new DataFrame

def identify_bijective_relationships(df):
	""" Identify columns that map onto one another in a one-to-one fashion, even if the string differ """
	cols_to_drop = []
	mapping_tables=[]
	num_cols = len(df.columns)
	for i in range(num_cols):
		for j in range(i+1, num_cols):
			col_i = df.iloc[:, i]
			col_j = df.iloc[:, j]
			if col_i.dropna().unique().size == col_j.dropna().unique().size == col_i.dropna().nunique() == col_j.dropna().nunique() and col_i.isin(col_j).all() and col_j.isin(col_i).all():
				cols_to_drop.append(df.columns[j])
				mapping_tables=ColToMetaTable(cols_to_drop)
	return cols_to_drop, mapping_tables

# Check if every department is represented
def identify_surjective_relationships(df):
    # Get all column names
    columns = df.columns
    # Prepare a dictionary to store results for each pair
    results = {}
    
    # Iterate over all pairs of columns
    for i in range(len(columns)):
        for j in range(len(columns)):
            if i != j:
                # Get column names
                domain_col = columns[i]
                codomain_col = columns[j]
                # Check if all values in the codomain are covered by the domain
                surjective = all(item in df[domain_col].unique() for item in df[codomain_col].unique())
                # Store the result in the dictionary
                results[(domain_col, codomain_col)] = surjective
    
    return results

def main():
	# Sample data setup
	df = pd.DataFrame({
		'A': [' Apple ', ' Banana ', 'Citrus ', 'Apple'],
		'B': ['apple', 'banana', 'citrus', 'apple'],
		'C': ['1', '2', '3', '1'],
		'D': ['one', 'two', 'three', 'one']
	})

	# Identify candidates that could possibly serve as a primary key
	dedup_df, dropped_columns = Identify_Primary_Key_Cols(df)
	print("Dropped Columns:", dropped_columns, '\n', "De-duplicated DataFrame:\n", dedup)

	# Preprocess the DataFrame
	duplicate_df = Duplicate_And_PreProcess_Db(df)
	print("duplicate DataFrame in all lower case; all leading and trailing white spaces removed:\n", duplicate_df)

	# Identify and remove duplicate columns
	dedup_df, dropped_columns = Identify_and_Remove_Duplicate_Cols(dedup_df)
	print("Dropped Columns:", dropped_columns, '\n', "De-duplicated DataFrame:\n", dedup)

	# Like Identical Columns, one of the two of a pair of bijectives can be removed as long as the mapping is stored
	proc_df, oneToOne_columns = Bijectives_to_LookUp(processed_df)
	print("One-to-one Columns:", oneToOne_columns, '\n', "Bijective Relationships removed to increase LookUp Speed (data footprint unchanged).\n", final_df)

  # Surjective columns to LookUp Table
  proc_df, oneToOne_columns = Surjective_Cols_To_LookUp(processed_df)
  print("One-to-one Columns:", oneToOne_columns, '\n', "Surjective Relationships removed to decrease data footprint.\n", final_df)

	# Now that the 

if __name__ == "__main__":
	main()

