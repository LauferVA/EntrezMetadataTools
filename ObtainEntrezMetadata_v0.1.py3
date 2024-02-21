import csv
import xml.etree.ElementTree as ET
from Bio import Entrez
import json
import logging
import time
import pandas as pd
import numpy as np
import os
import subprocess

# Set up logging
logging.basicConfig(filename='entrez_errors.log', level=logging.ERROR, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')

with open('../EntrezCredentials.txt', 'r') as UserData:
    EntrezInfo = UserData.readline().split()
    Entrez.email = EntrezInfo[0]
    Entrez.api_key = EntrezInfo[1]
    Entrez.tool = EntrezInfo[2]

def check_and_write_output_file(output_file, record, total_records):
    """
    Checks if the output file already has the expected number of records.
    Writes to the file only if the existing number of records is less than total_records.
    """
    # Check if the file exists and has the expected number of lines/records
    if os.path.exists(output_file):
        with open(output_file, "r") as file:
            lines = file.readlines()
            # Assuming each record is written in a separate line
            if len(lines) >= total_records: 
                print(f"File {output_file} already has the expected number of records or more.")
                return False  # Indicates no need to write again
    # If file does not exist or has fewer records, proceed to write
    with open(output_file, "a") as file:
        file.write(json.dumps(record, indent=2))
        file.write("\n\n")  # Add spacing between entries for readability
    return True  # Indicates that writing was performed

def clean_and_write_df(df, output_file_path, threshold_percentage=5):
    """
    Removes columns with non-null row less than % threshold from DF; writes the cleaned DF to output_file_path.
    - df: The pandas DataFrame to clean and write.
    - output_file_path: Path to the output tab-delimited file.
    - threshold_percentage: Minimum percentage of non-null entries required to keep a column.
    """
    df = df.replace("", np.nan)     # Replace empty strings with NaN to treat them as missing values
    percentage_non_null = df.notnull().mean() * 100 # Calculate the percentage of non-null entries for each column
#
    columns_to_drop = percentage_non_null[percentage_non_null < threshold_percentage].index      # Identify columns with % non-null entries < the threshold
    df_cleaned = df.drop(columns=columns_to_drop)     # Drop these columns from the DataFrame
    df_cleaned.to_csv(output_file_path, sep='\t', index=False)     # Write the cleaned DataFrame to a tab-delimited file
    print(f"Cleaned DataFrame written to {output_file_path}. Dropped columns: {list(columns_to_drop)}")

    return df

# def parse_xml_to_dict(xml_string):
#     try:
#         # Parse the XML string
#         root = ET.fromstring(xml_string)
#         # Initialize a dictionary to hold the parsed data
#         parsed_data = {}
#         # Iterate through all elements in the XML
#         for element in root.iter():
#             # Use the element tag as the key and the element text as the value
#             # If the element has attributes, add them as keys with their respective values
#             if element.text:
#                 parsed_data[element.tag] = element.text.strip()
#             for name, value in element.attrib.items():
#                 parsed_data[name] = value
#         return parsed_data
#     except ET.ParseError:
#         return {}  # Return an empty dictionary if parsing fails

# def parse_xml_columns(df):
#     for column in df.columns:
#         # Check if the column content looks like XML (e.g., starts with '<' and ends with '>')
#         if df[column].dtype == object and df[column].str.startswith('<').all():
#             # Parse the XML content in the column
#             parsed_columns = df[column].apply(parse_xml_to_dict).apply(pd.Series)
#             # Drop the original XML column
#             df = df.drop(column, axis=1)
#             # Concatenate the parsed columns to the original DataFrame
#             df = pd.concat([df, parsed_columns], axis=1)
#     return df

def fetch_ids_in_batches_and_write_record(db, term, batch_size=10000):
    print(f"attempting to fetch IDs from {db}.")
    handle = Entrez.esearch(db=db, term=term, retmax=0)
    record = Entrez.read(handle)
    total_records = int(record["Count"])
    output_file = f"../Data/{db}_summary.json"

    ids = []
    for start in range(0, total_records, batch_size):
        handle = Entrez.esearch(db=db, term=term, retstart=start, retmax=batch_size)
        record = Entrez.read(handle)
        ids.extend(record["IdList"])
        handle.close()

        # Before writing, check if the output file already meets the expected record count
        if not check_and_write_output_file(output_file, record, total_records):
            break  # Skip writing if the file already has the expected number of records or more

    print(f"done fetching {db} IDs")
    return ids


def fetch_summaries_in_batches_to_df(db, id_list, batch_size=10000, max_batches_per_cycle=25):
    total_batches = (len(id_list) - 1) // batch_size + 1
    cycles = (total_batches - 1) // max_batches_per_cycle + 1
    
    for cycle in range(cycles):
        print(f"Starting cycle {cycle + 1} of {cycles}")
        all_summaries = []  # Reinitialize for each cycle
        
        batch_start = cycle * max_batches_per_cycle
        batch_end = min((cycle + 1) * max_batches_per_cycle, total_batches)
        
        for batch_index in range(batch_start, batch_end):
            start = batch_index * batch_size
            batch_ids = id_list[start:start + batch_size]
            handle = Entrez.esummary(db=db, id=",".join(batch_ids))
            summaries = Entrez.read(handle, validate=False)
            handle.close()

            if 'DocumentSummarySet' in summaries:
                all_summaries.extend(summaries['DocumentSummarySet']['DocumentSummary'])
            else:
                all_summaries.extend(summaries)

            print(f"Processed batch {batch_index + 1} of {total_batches} in cycle {cycle + 1}")
            
        df = pd.DataFrame(all_summaries)
#        df = parse_xml_columns(df)  # Here, insert the new function to check and parse XML-like content
        output_file_path = f"../Data/{db}_metadata_cycle_{cycle + 1}.tsv"
        df_cleaned = clean_and_write_df(df, output_file_path, threshold_percentage=5)

        print(f"Cycle {cycle + 1} completed. Data written to {output_file_path}.")



# Set search terms
DbBioProject = "bioproject"; DbBioSample = "biosample"; DbSra='sra'
query='("Homo sapiens"[Organism] OR "Mus musculus"[Organism]) AND RNA-Seq[All Fields]'

# Usage for BioProject
BioProjectIds = fetch_ids_in_batches_and_write_record(DbBioProject, query, batch_size=10000)
BioprojectDf = fetch_summaries_in_batches_to_df(DbBioProject, BioProjectIds)

# Usage for BioSample
BioSampleIds = fetch_ids_in_batches_and_write_record(DbBioSample, query, batch_size=10000)
BioSampleDf = fetch_summaries_in_batches_to_df(DbBioSample, BioSampleIds)

# Usage for Sra
SraIds = fetch_ids_in_batches_and_write_record(DbSra, query, batch_size=10000)
SraDf = fetch_summaries_in_batches_to_df(DbSra, SraIds)


# Define the command to be executed
command = "awk 'FNR==1 && NR!=1{next;}{print}' sra_metadata_cycle_{1..12}.tsv > sra_metadata.tsv"

# Use the subprocess.run method to execute the command
result = subprocess.run(command, shell=True, capture_output=True, text=True)
print(result)

