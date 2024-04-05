import json
import pandas as pd
from Bio import Entrez
from bs4 import BeautifulSoup
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

print('\n\n'

def printAndLog(message, level="info"):
    """
    Prints and logs a message at the specified level.
    :param message: The message to print and log.
    :param level: The logging level ('info', 'warning', 'error').
    """
    print(message)  # Print to console
    if level.lower() == "info":
        logging.info(message)
    elif level.lower() == "warning":
        logging.warning(message)
    elif level.lower() == "error":
        logging.error(message)
    else:
        logging.info(message)  # Default to info level if unspecified/incorrect level

# Override the default HTTPS context
def create_https_context():
  context = ssl.create_default_context(cafile=certifi.where())
  return context

ssl._create_default_https_context = create_https_context

def handle_http_errors(func):
    def wrapper(*args, **kwargs):
        retries = 3  # Number of retries for IncompleteRead
        for attempt in range(retries):
            try:
                return func(*args, **kwargs)
            except HTTPError as e:
                printAndLog(f"HTTP error occurred: {e}", "error")
                if e.code == 400:
                    printAndLog("Bad request. Please check the parameters or request syntax.", "error")
                    break  # No retry for client errors
                elif e.code == 429:
                    printAndLog("Too many requests. Sleeping for 10 seconds.", "warning")
                    time.sleep(10)  # Retry after waiting
                elif 500 <= e.code < 600:
                    printAndLog("Server error. Retrying the request.", "warning")
                    if attempt < retries - 1:
                        time.sleep(2 ** (attempt + 1))  # Exponential backoff for server errors
                    else:
                        printAndLog("Final attempt failed.", "error")
                        raise
                else:
                    raise
            except IncompleteRead as e:
                # Handle IncompleteRead separately to potentially retry
                printAndLog(f"Incomplete read error: {e}, bytes read: {e.partial}", "warning")
                if attempt < retries - 1:
                    printAndLog(f"Retrying due to incomplete read... Attempt {attempt + 1} of {retries}", "info")
                    time.sleep(2 ** (attempt + 1))  # Exponential backoff
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
def fetch_metadata_for_batch(db, batch_ids, index, total_batches):
    try:
        # Attempt fetching with the original batch size
        records = Entrez.efetch(db=db, id=','.join(batch_ids), rettype="xml").read()
        return records
    except HTTPError as e:
        if e.code == 400:
            # If Bad Request, attempt to handle by splitting the batch
            mid = len(batch_ids) // 2
            if mid > 0:
                printAndLog(f"Splitting batch {index+1}/{total_batches} due to HTTP Error 400.", "warning")
                first_half = fetch_metadata_for_batch(db, batch_ids[:mid], index, total_batches)
                second_half = fetch_metadata_for_batch(db, batch_ids[mid:], index, total_batches)
                return first_half + second_half
            else:
                printAndLog(f"Failed to fetch batch {index+1}/{total_batches}: {e}", "error")
                return ""
        else:
            raise

def parallel_fetch_metadata(primary_ids, db, batch_size=500):
    total_records = len(primary_ids)
    total_batches = (total_records + batch_size - 1) // batch_size  # Ceiling division
    printAndLog(f"Total number of records to process: {total_records}\n")

    with concurrent.futures.ThreadPoolExecutor(max_workers=25) as executor:
        futures = []
        for index, start in enumerate(range(0, len(primary_ids), batch_size)):  # Index usage
            end = min(start + batch_size, len(primary_ids))
            batch_ids = primary_ids[start:end]
            futures.append(executor.submit(fetch_metadata_for_batch, db, batch_ids, index, total_batches)) # Pass index and total_batches to the fetch_metadata_for_batch()
        
        results = []
        for future in concurrent.futures.as_completed(futures):
            result = future.result()
            if result:
                results.append(result)
    return results

def flatten_xml(element, path="", flat_dict=None):
    if flat_dict is None:
        flat_dict = {}
    if len(element.find_all(recursive=False)) > 0:
        for child in element.find_all(recursive=False):
            child_path = f"{path}/{child.name}" if path else child.name
            flatten_xml(child, path=child_path, flat_dict=flat_dict)
    else:
        flat_dict[path] = element.text.strip()
    return flat_dict

def parse_and_flatten_xml(xml_records):
    all_records = []
    for record in xml_records:
        try:
            soup = BeautifulSoup(record, 'lxml')            # First attempt to parse with 'lxml'
        except Exception as lxml_error:
            printAndLog(f"lxml parsing failed: {lxml_error}, trying 'html5lib'", "warning")
            try:
                soup = BeautifulSoup(record, 'html5lib')                # Then fall back to 'html5lib' if 'lxml' fails
            except Exception as html5lib_error:
                printAndLog(f"html5lib parsing failed: {html5lib_error}, falling back to 'html.parser'", "warning")
                soup = BeautifulSoup(record, 'html.parser')                # Last resort: use Python’s built-in HTML parser
        flat_dict = flatten_xml(soup)
        all_records.append(flat_dict)
    return all_records

def write_to_dataframe(records, output_file):
    df = pd.DataFrame(records)
    common_prefix = os.path.commonprefix(df.columns.tolist())
    df.rename(columns=lambda x: x.replace(common_prefix, '').replace('/', '_'), inplace=True)
    df.dropna(thresh=len(df) * 0.95, axis=1, inplace=True)
    directory = os.path.dirname(output_file)
    if not os.path.exists(directory):
        os.makedirs(directory)
    df.to_csv(output_file, index=False, sep='\t')

def main():
    start_time = datetime.now()
    printAndLog(f"Script started at {start_time.strftime('%Y-%m-%d %H:%M:%S')}")

    db = 'bioproject'
    credentialFile = '../inputFiles/Auth/EntrezCredentials.json'
    file_path_prefix = '/data/user/vlaufer/DeathStar'
    json_file_path = f"{file_path_prefix}/Data/esearch/{db}/{db}_Hs_Ids.json"
    output_file = f"{file_path_prefix}/Data/efetch/{db}/{db}_hS_metadata.tsv"
    #### start timer
    credentials = read_credentials(credentialFile)
    maxWorkers = authenticate_with_entrez(credentials)
    primary_ids = read_primary_ids(json_file_path)
    xml_records = parallel_fetch_metadata(primary_ids, db)
    records_dict = parse_and_flatten_xml(xml_records)
    df = pd.concat([pd.DataFrame([record]) for record in records_dict], ignore_index=True)
    write_to_dataframe(df, output_file)
    #### end timer
    end_time = datetime.now()
    printAndLog(f"Script ended at {end_time.strftime('%Y-%m-%d %H:%M:%S')}")
    elapsed_time = end_time - start_time
    printAndLog(f"Total execution time: {elapsed_time}")

if __name__ == "__main__":
    main()
