import json
import pandas as pd
from Bio import Entrez
from urllib.error import HTTPError
import time
import os
import concurrent.futures
import certifi
import ssl
import logging
from datetime import datetime
from http.client import IncompleteRead

# Configure logging to write to a file and print to console.
log_directory = "../Logs"
os.makedirs(log_directory, exist_ok=True)
log_file_path = os.path.join(log_directory, 'ncbi_fetch_log.txt')

# Configure logging to write to a file and print to console.
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(levelname)s: %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    handlers=[logging.FileHandler(log_file_path),
                              logging.StreamHandler()])

print('\n\n')

def printAndLog(message, level="info"):
    print(message)  # Print to console
    if level.lower() == "info":
        logging.info(message)
    elif level.lower() == "warning":
        logging.warning(message)
    elif level.lower() == "error":
        logging.error(message)

# Override the default HTTPS context
def create_https_context():
  context = ssl.create_default_context(cafile=certifi.where())
  return context

ssl._create_default_https_context = create_https_context

def handle_http_errors(func):
    def wrapper(*args, **kwargs):
        retries = 4
        delay = 1.5  # Initial delay in seconds
        for attempt in range(retries):
            try:
                return func(*args, **kwargs)
            except HTTPError as e:
                if e.code == 429:
                    printAndLog(f"429 Too Many Requests: Retrying in {delay} seconds...", "warning")
                    time.sleep(delay)
                    delay *= 2  # Double the delay for the next attempt
                else:
                    printAndLog(f"HTTP error occurred: {e}", "error")
                    if e.code in [400, 500, 502, 503, 504]:
                        time.sleep(delay)  # Wait before retrying
                        delay *= 2  # Increase delay for the next attempt
                    else:
                        break  # Don't retry for other HTTP errors
            except IncompleteRead as e:
                printAndLog(f"Incomplete read error: {e}, bytes read: {e.partial}", "warning")
                if attempt < retries - 1:
                    printAndLog(f"Retrying due to incomplete read... Attempt {attempt + 1} of {retries}", "info")
                    time.sleep(delay)
                else:
                    printAndLog("Final attempt failed due to incomplete read.", "error")
                    raise
    return wrapper

def read_credentials(credentialFile):
    with open(credentialFile, 'r') as file:
        credentials = json.load(file)
    return credentials

@handle_http_errors
def authenticate_with_entrez(credentials):
    Entrez.email = credentials['email']
    Entrez.tool = credentials['projectName']
    Entrez.api_key = credentials['apiKey']
    try:
        handle = Entrez.esearch(db="pubmed", term="NCBI", retmax=1)
        handle.close()
        return 9
    except HTTPError as e:
        if e.code == 401:
            return 3
        else:
            raise

def read_primary_ids(json_file_path):
    with open(json_file_path, 'r') as file:
        data = json.load(file)
    return data

@handle_http_errors
def fetch_metadata_for_batch(db, batch_ids, index, total_batches, output_dir):
    try:
        # Attempt fetching with the original batch size
        records = Entrez.efetch(db=db, id=','.join(batch_ids), rettype="xml").read()
        if isinstance(records, bytes):
            records = records.decode('utf-8')  # Decode bytes to string using utf-8 or appropriate encoding
        with open(os.path.join(output_dir, f'batch_{index+1}.xml'), 'w') as file:
            file.write(records)
        return f"Batch {index+1} saved."
    except HTTPError as e:
        if e.code == 400:
            # If Bad Request, attempt to handle by splitting the batch
            mid = len(batch_ids) // 2
            if mid > 0:
                printAndLog(f"Splitting batch {index+1}/{total_batches} due to HTTP Error 400.", "warning")
                first_half = fetch_metadata_for_batch(db, batch_ids[:mid], index, total_batches, output_dir)
                second_half = fetch_metadata_for_batch(db, batch_ids[mid:], index, total_batches, output_dir)
                return first_half + second_half
            else:
                printAndLog(f"Failed to fetch batch {index+1}/{total_batches}: {e}", "error")
                return ""
        else:
            raise

def parallel_fetch_metadata(primary_ids, db, batch_size, output_dir):
    total_records = len(primary_ids)
    total_batches = (total_records + batch_size - 1) // batch_size  # Ceiling division
    printAndLog(f"Total number of records to process: {total_records}\n")

    with concurrent.futures.ThreadPoolExecutor(max_workers=25) as executor:
        futures = []
        for index, start in enumerate(range(0, len(primary_ids), batch_size)):  # Index usage
            time.sleep(3)
            end = min(start + batch_size, len(primary_ids))
            batch_ids = primary_ids[start:end]
            futures.append(executor.submit(fetch_metadata_for_batch, db, batch_ids, index, total_batches, output_dir)) # Pass index and total_batches to the fetch_metadata_for_batch()
        
        results = []
        for future in concurrent.futures.as_completed(futures):
            result = future.result()
            if result:
                results.append(result)
    return results

def main():
    start_time = datetime.now()
    printAndLog(f"Script started at {start_time.strftime('%Y-%m-%d %H:%M:%S')}")

    db = 'sra'
    credentialFile = '../inputFiles/Auth/EntrezCredentials.json'
    file_path_prefix = '/data/user/vlaufer/DeathStar'
    json_file_path = f"{file_path_prefix}/Data/esearch/{db}/{db}_Hs_Ids.json"
    output_dir = f"{file_path_prefix}/Data/efetch/{db}"
    os.makedirs(output_dir, exist_ok=True)
    batch_size = 5000  # You might adjust this based on your API limitations and testing

    credentials = read_credentials(credentialFile)
    maxWorkers = authenticate_with_entrez(credentials)
    primary_ids = read_primary_ids(json_file_path)
    results = parallel_fetch_metadata(primary_ids, db, batch_size, output_dir)

    end_time = datetime.now()
    printAndLog(f"Script ended at {end_time.strftime('%Y-%m-%d %H:%M:%S')}")
    elapsed_time = end_time - start_time
    printAndLog(f"Total execution time: {elapsed_time}")

if __name__ == "__main__":
    main()
