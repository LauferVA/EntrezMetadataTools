import json
from Bio import Entrez

def fetch_entrez_info(credentials, project_info):
    """
    Fetches Entrez database information using the einfo tool and writes it to a JSON file.

    Parameters:
    - databases: List of strings, each an Entrez database name.
    - output_file: String, the filename to which the JSON data will be saved.
    """
    Entrez.email = "your_email@example.com"  # Replace with your email address.
    Entrez.api_key = "Your_Entrez_API_Key"   # Replace with your actual API key if you have one.

    info_dict = {}
    for db in databases:
        with Entrez.einfo(db=db) as handle:
            db_info = Entrez.read(handle)
            info_dict[db] = db_info

    output_file='../Logs/QueryID_{hs}.json'
    with open(output_file, 'w') as f:
        json.dump(info_dict, f, indent=2)

    print(f"Information written to {output_file}")

# Example usage:
databases = ["biosamples", "bioproject", "gds", "pubmed", "sra"]
fetch_entrez_info(databases, "entrez_info.json")
