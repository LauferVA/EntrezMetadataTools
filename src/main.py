
from entrez_handler import EntrezHandler
from config import api_key

""" Place holder for main.py to come """

def main():
    # Create an instance of the EntrezHandler class
    handler = EntrezHandler(api_key=api_key)
    
    # Example function calls
    raw_data = handler.download_metadata('pubmed', ['12345', '67890'])
    parsed_data = handler.parse_xml(raw_data)
    cleaned_data = handler.clean_data(parsed_data)
    merged_data = handler.merge_data([cleaned_data])
    results = handler.interrogate_data(merged_data)
    
    print(results)  # or handle results in a different way

if __name__ == "__main__":
    main()
