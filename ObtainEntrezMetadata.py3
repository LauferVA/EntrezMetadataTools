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
import shortuuid
import threading

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
    """Check and write to the output file if it doesn't have enough records."""
    try:
        # Attempt to open for reading and writing; this assumes the file exists.
        with open(output_file, "r+") as file:
            lines = file.readlines()
            if len(lines) >= total_records:
                print(f"File {output_file} already has the expected number of records or more.")
                return False
            else:
                # Move the cursor to the end of the file to append
                file.seek(0, os.SEEK_END)
                file.write(json.dumps(record, indent=2) + "\n\n")
    except FileNotFoundError:
        # File doesn't exist, so create it and write the record
        with open(output_file, "w") as file:
            file.write(json.dumps(record, indent=2) + "\n\n")
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
    print(f"Fetching the query records from {db}.")
    handle = Entrez.esearch(db=db, term=term, retmax=0)
    record = Entrez.read(handle)
    # print(f"An error occurred while fetching IDs: {e}")
    # logging.error(f"Failed to fetch IDs for {db} with term {term}: {e}", exc_info=True)
    handle.close()
    total_records = int(record["Count"])
    output_file = f"../Data/{db}_summary.json"

    ids = []
    for start in range(0, total_records, batch_size):
        print(f"Initial annotation of {start} records is complete.")
        try:
            handle = Entrez.esearch(db=db, term=term, retstart=start, retmax=batch_size)
            record = Entrez.read(handle)
        except Exception as e:
            print(f"An error occurred in the while reading batched records: {e}")
            logging.error(f"Failed to read IDs for {db} with term {term}: {e}", exc_info=True)
        ids.extend(record["IdList"])
        handle.close()
        if not check_and_write_output_file(output_file, record, total_records):
            break
    print(f"Done fetching {db} IDs")
    return ids

def fetch_summary_batch(db, id_list, retries=5, backoff=1.5):
    """Fetch a batch of summaries from NCBI Entrez database."""
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

def dump_dataframe_to_file(df, cycle_number):
    """Dump the DataFrame to a CSV file."""
    unique_id = shortuuid.ShortUUID().random(length=10) 
    file_name = f"../Data/{db}_metadata_cycle_{cycle_number}_{unique_id}.tsv"
    df.to_csv(file_name, index=False)
    print(f"A group of batched queries of length {len(df)} records has been written to {file_name}.")


def parallel_fetch_summaries(db, id_list, api_key, batch_size=10000, max_batches_per_cycle=16):
    total_batches = (len(id_list) - 1) // batch_size + 1
    cycles = (total_batches - 1) // max_batches_per_cycle + 1
    worker_count = get_worker_count(api_key)  # Assuming this function returns the number of workers based on the API key
    
    for cycle in range(cycles):
        print(f"\nStarting cycle {cycle + 1} of {cycles} at {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}")
        df_accumulator = pd.DataFrame()  # Initialize an empty DataFrame for each cycle
        
        batch_start_index = cycle * max_batches_per_cycle * batch_size
        batch_end_index = min((cycle + 1) * max_batches_per_cycle * batch_size, len(id_list))
        current_cycle_batches = [id_list[i:i + batch_size] for i in range(batch_start_index, batch_end_index, batch_size)]
        
        with ThreadPoolExecutor(max_workers=worker_count) as executor:
            future_to_batch = {executor.submit(fetch_summary_batch, db, batch): batch for batch in current_cycle_batches}
            for future in as_completed(future_to_batch):
                try:
                    records = future.result()
                    df_batch = pd.DataFrame(records)  # Assuming records is a list of dictionaries
                    df_accumulator = pd.concat([df_accumulator, df_batch], ignore_index=True)
                    print(f"Processed a batch within cycle {cycle + 1}")
                except Exception as e:
                    logging.error(f"An error occurred while processing a batch: {e}")

        # After all batches in a cycle are processed, dump the accumulated DataFrame to file
        if not df_accumulator.empty:
            dump_dataframe_to_file(df_accumulator, cycle + 1)  # Assuming this function takes a DataFrame and a file path

# Initialize global variables to manage the process
current_file_index = 0          # To keep track of the file number being written
records_processed = 0           # To track the number of records processed across batches
batch_size = 10000              # Number of records pulled in an any single API call.
batches_per_cycle = 10         # Threshold to dump data to file
db = 'sra'
query = '("Homo sapiens"[Organism] OR "Mus musculus"[Organism]) AND RNA-Seq[All Fields]'

# Configure logging
logging.basicConfig(filename=f'../Out/{db}_metadata_extraction_errors.log', level=logging.ERROR, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')

# Main logic to use the functions
def main():
    # get an initial start time for benchmarking purposes
    print(f'Begin extracting metadata from {db}')
    print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    logging.error(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), exc_info=True)

    # Set your Entrez email, tool, and api_key
    # Define your databases and query
    # Fetch IDs, summaries, and process data as needed
    with open('../EntrezCredentials.txt', 'r') as UserData:
        EntrezInfo = UserData.readline().split()
        Entrez.email = EntrezInfo[0]
        Entrez.tool = EntrezInfo[1]
        Entrez.api_key = EntrezInfo[2]

    # Fetch IDs based on the query
    queryIds = fetch_ids(db, query)

    # Fetch summaries in parallel and manage dumping to files
    parallel_fetch_summaries(db, queryIds, Entrez.api_key, batch_size, batches_per_cycle)

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logging.error("An error occurred during execution - no more specific process identified if not listed above.", exc_info=True)

    print('All processes complete')
    print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    logging.error(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), exc_info=True)
