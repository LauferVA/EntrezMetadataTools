import config
import libs
import utils

from Bio import Entrez

import certifi
import ssl

# Override the default HTTPS context
def create_https_context():
    context = ssl.create_default_context(cafile=certifi.where())
    return context

ssl._create_default_https_context = create_https_context # Get secure connection

class NCBI_Metadata_Aggregator:
    def __init__(self, query_params):
        self.email=settings.credentials['email']
        self.project=settings.credentials['project']
        self.apikey=settings.credentials['apikey']

        self.organism = settings.project_info['organism']
        self.taxid = settings.project_info['taxid']
        self.org_abbrev = settings.project_info['org_abbrev']
        self.query = settings.project_info['query']
        self.db = settings.project_info['db']

    @log_decorator
    @handle_errors(default_value=[])
    def download_ids(self):
        return fetch_ids(self.db, self.query)

    @log_decorator
    @handle_errors(default_value={})
    def download_id_metadata(self, ids):
        return ids_to_records(ids)

    @log_decorator
    @handle_errors(default_value="")
    def xml_to_fi(self, data):
        return xml_to_tsv(data)  # Assuming xml_to_tsv is correctly named and should be xml_to_data or similar

    @log_decorator
    @handle_errors(default_value="")
    def convert_to_tsv(self, data):
        return dictionary_to_tsv(data)

    @log_decorator
    @handle_errors(default_value="")
    def aggregate_and_process(self, ids):
        records = self.get_records(ids)
        merged_data = merge_metadata(records)
        cleaned_data = clean_metadata(merged_data)
        return cleaned_data


def main():
    aggregator = NCBI_Metadata_Aggregator()
    ids = aggregator.get_ids()
    processed_data = aggregator.aggregate_and_process(ids)
    tsv_data = aggregator.convert_to_tsv(processed_data)
    xml_data = aggregator.convert_to_xml(processed_data)

    print("TSV Data:\n", tsv_data)
    print("XML Data:\n", xml_data)

if __name__ == "__main__":
    main()
