"""
This script is designed to interface with a biomedical database API, offering several functionalities:


efetch: Retrieves records in the requested format from a list of one or more primary IDs or from the user's environment. Useful for fetching detailed information about specific entries.
epost: Posts a file containing a list of primary IDs for future use in the user's environment, enabling batch operations on these IDs with subsequent search strategies.
esearch: Searches and retrieves primary IDs for use in other operations like EFetch, ELink, and ESummary. It can also retain results for future use in the user's environment, aiding in complex search strategies.
elink: Checks for the existence of external or Related Articles links from a list of one or more primary IDs. It can retrieve primary IDs and relevancy scores for links to Entrez databases or Related Articles, facilitating the discovery of related content.
einfo: Provides field index term counts, last update, and available links for each database, which is helpful for understanding database structure and content coverage.
esummary: Retrieves document summaries from a list of primary IDs or from the user's environment, offering a quick overview of the entries without fetching complete records.
egquery: Provides Entrez database counts in XML for a single search using Global Query, allowing for cross-database searches to understand the spread of information across the system.
espell: Retrieves spelling suggestions, improving search accuracy by suggesting corrections to potentially misspelled search terms.
ecitmatch: Retrieves PubMed IDs (PMIDs) that correspond to a set of input citation strings, enabling the linking of bibliographic references to database entries.

The script requires further implementation to make API calls, parse responses, and handle data according to specific requirements.
"""

from Bio import Entrez
import json
import os
from urllib.error import HTTPError  # If using Bio.Entrez which uses urllib under the hood
from functools import wraps
import time
import certifi
import ssl

# Override the default HTTPS context
def create_https_context():
  context = ssl.create_default_context(cafile=certifi.where())
  return context

ssl._create_default_https_context = create_https_context

def retry_on_exception(retries=5, backoff=1.5, exceptions=(Exception,)):
    """ A decorator to retry a function call with exponential backoff. """
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            attempts = 0
            while attempts < retries:
                try:
                    return func(*args, **kwargs)
                except exceptions as e:
                    attempts += 1
                    print(f"An error occurred: {e}. Attempt {attempts} of {retries}.")
                    time.sleep(backoff ** attempts)
            raise Exception(f"Failed after {retries} attempts.")
        return wrapper
    return decorator

# Define the function needed by this workflow.
def authenticate_with_ncbi(credentialInfo):
    """As per https://www.ncbi.nlm.nih.gov/books/NBK25497/"""
    Entrez.email = credentialInfo['email']
    Entrez.tool = credentialInfo['projectName']
    Entrez.api_key = credentialInfo['apiKey']
    try:
        handle = Entrez.einfo()  # A simple request to test the API key
        record = Entrez.read(handle)
        handle.close()
        print("API key is valid.")
        return 9
    except RuntimeError as e:
        print(f"Error: {e} API key may be invalid incorrectly entered, etc.")
        return 3

with open('../inputFiles/Auth/EntrezCredentials.json', 'r') as auth_file:
    credentials = json.load(auth_file) # {'email': 'vlaufer@med.umich.edu', 'projectName': 'The100KTranscriptomesProject', 'apiKey': 'a572ea77509b246c85fc9c0a59e48dc70908'}

worker_count=authenticate_with_ncbi(credentials)

@retry_on_exception(retries=5, backoff=1.5, exceptions=(HTTPError, Exception))
def dbIds2json(dbs, term, file_path_prefix, retmax=2000000):
    if isinstance(dbs, str):
        dbs = [dbs]  # Convert a single database name into a list for uniform handling
    result_paths = []  # Initialize an empty list to store file paths of the results
    for db in dbs:
        file_path = f"{file_path_prefix}/{db}/{db}_Hs_Ids.json"  # Prepare file path, ensuring unique paths for different dbs
        
        if os.path.exists(file_path):
            print(f"File already exists for {db}, skipping to next database.")
            continue  # Skip to the next database
        handle = Entrez.esearch(db=db, term=term, retmax=retmax)  # Fetch IDs
        record = Entrez.read(handle)
        handle.close()
        ids = record['IdList']
        os.makedirs(os.path.dirname(file_path), exist_ok=True)        # Create directories if they do not exist
        with open(file_path, 'w') as file:
            json.dump(ids, file, indent=2)         # Write data to JSON file
        print(f"Data successfully written to {file_path}.")
        result_paths.append(file_path)  # Append the file path to the results list
    return result_paths  # Return the list of paths where data has been written

# Run time
db = ['sra', 'gds', 'biosample', 'bioproject']
term = '"Homo sapiens"[Organism] AND RNA-Seq[All Fields]'
file_path_prefix = "../Data/esearch"
db_Ids=dbIds2json(db, term, file_path_prefix)
