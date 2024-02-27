# import packages needed for file I/O
import os
import subprocess
import csv
import json

# import packages for data storage in memory
import pandas as pd
import numpy as np

# import packages needed to parallelize
import logging
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

# error handling
from urllib.error import HTTPError

# import Bio and XML parsing packages to actually do the work
from Bio import Entrez
import xml.etree.ElementTree as ET

# Define the function needed by this workflow.
def get_worker_count(api_key):
    """As per https://www.ncbi.nlm.nih.gov/books/NBK25497/"""
    return 9 if api_key else 3

def check_and_write_output_file(output_file, record, total_records):
    """Check and write output file if needed."""
    if os.path.exists(output_file):
        with open(output_file, "r") as file:
            lines = file.readlines()
            if len(lines) >= total_records:
                print(f"File {output_file} already has the expected number of records or more.")
                return False
    with open(output_file, "a") as file:
        file.write(json.dumps(record, indent=2))
        file.write("\n\n")
    return True

def clean_and_write_df(df, output_file_path, threshold_percentage=5):
    """Clean and write DataFrame to file."""
    df = df.replace("", np.nan)
    percentage_non_null = df.notnull().mean() * 100
    columns_to_drop = percentage_non_null[percentage_non_null < threshold_percentage].index
    df_cleaned = df.drop(columns=columns_to_drop)
    df_cleaned.to_csv(output_file_path, sep='\t', index=False)
    print(f"Cleaned DataFrame written to {output_file_path}. Dropped columns: {list(columns_to_drop)}")
    return df_cleaned

def fetch_ids(db, term, batch_size=10000):
    print(f"Attempting to fetch IDs from {db}.")
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

    print(f"Done fetching {db} IDs")
    return ids

def fetch_summary_batch(db, id_list, retries=5, backoff=1.5):
    """Fetch a batch of summaries from NCBI Entrez database with error handling."""
    attempt = 0
    while attempt < retries:
        try:
            handle = Entrez.esummary(db=db, id=",".join(id_list))
            records = Entrez.read(handle)
            handle.close()
            return records
        except HTTPError as e:
            print(f"HTTP Error encountered: {e.code} {e.reason}. Attempt {attempt+1} of {retries}.")
            time.sleep(backoff ** attempt)  # Exponential backoff
            attempt += 1
        except Exception as e:
            print(f"An error occurred: {e}. Attempt {attempt+1} of {retries}.")
            time.sleep(backoff ** attempt)  # Exponential backoff
            attempt += 1
    raise Exception(f"Failed to fetch summaries after {retries} attempts.")

def dump_dataframe_to_file(df):
    """Dump the DataFrame to a CSV file."""
    global current_file_index
    current_file_index += 1
    file_name = f'summary_batch_{current_file_index}.csv'
    df.to_csv(file_name, index=False)
    print(f"Dumped {len(df)} records to {file_name}.")

def parallel_fetch_summaries(db, id_list, batch_size=10000):
    """Fetch summaries in parallel and manage dumping to files at specified thresholds."""
    global records_processed
    worker_count = get_worker_count(Entrez.api_key)
    batches = [id_list[i:i + batch_size] for i in range(0, len(id_list), batch_size)]
    params = [(db, batch) for batch in batches]

    df_accumulator = pd.DataFrame()  # Initialize an empty DataFrame to accumulate fetched summaries

    with ThreadPoolExecutor(max_workers=worker_count) as executor:
        future_to_batch = {executor.submit(fetch_summary_batch, db, batch): batch for batch in batches}
        for future in as_completed(future_to_batch):
            records = future.result()
            df_batch = pd.DataFrame(records)  # Convert list of dictionaries to DataFrame
            df_accumulator = pd.concat([df_accumulator, df_batch], ignore_index=True)
            records_processed += len(df_batch)

            # Check if the threshold is reached to dump data to file and reset DataFrame
            if records_processed >= dump_threshold:
                dump_dataframe_to_file(df_accumulator)
                df_accumulator = pd.DataFrame()  # Reset DataFrame after dumping
                records_processed = 0  # Reset records processed counter

    # Dump any remaining records that didn't meet the threshold
    if not df_accumulator.empty:
        dump_dataframe_to_file(df_accumulator)

# Initialize global variables to manage the process
current_file_index = 0          # To keep track of the file number being written
records_processed = 0           # To track the number of records processed across batches
batch_size = 10000              # Number of records pulled in an any single API call.
dump_threshold = 100000         # Threshold to dump data to file

# Main logic to use the functions
def main():
    # Set your Entrez email, tool, and api_key
    # Define your databases and query
    # Fetch IDs, summaries, and process data as needed
    with open('../EntrezCredentials.txt', 'r') as UserData:
        EntrezInfo = UserData.readline().split()
        Entrez.email = EntrezInfo[0]
        Entrez.tool = EntrezInfo[1]
        Entrez.api_key = EntrezInfo[2]

    # Set search terms
    db = 'sra'
    query = '("Homo sapiens"[Organism] OR "Mus musculus"[Organism]) AND RNA-Seq[All Fields]'

    # Fetch IDs based on the query
    id_list = fetch_ids(db, query)

    # Fetch summaries in parallel and manage dumping to files
    parallel_fetch_summaries(db, id_list, batch_size=10000)

if __name__ == "__main__":
    # Configure logging
    logging.basicConfig(filename='entrez_errors.log', level=logging.ERROR, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    try:
        main()
    except Exception as e:
        logging.error("An error occurred during execution", exc_info=True)
