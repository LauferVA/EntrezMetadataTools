import pandas as pd
import subprocess


# input file paths
dD='../Data/'
SraMd=f"{dD}sra_metadata.tsv"
BpMd=f"{dD}bioproject_metadata_cycle_1.tsv"
BsMd=f"{dD}biosample_metadata_cycle_1.tsv"

# intermediate file paths
SraMdTemp=f"{dD}sra_temp.tsv"
BpMdTemp=f"{dD}bioproject_temp.tsv"
BsMdTemp=f"{dD}biosample_temp.tsv"

# output file paths:
SraMdOut=f"{dD}sra_metadata_proc.tsv"
BpMdOut=f"{dD}bioproject_metadata_proc.tsv"
BsMdOut=f"{dD}biosample_metadata_proc.tsv"

# Execute the awk command
xml_column_index=10
awk_cmd = f"awk -F'\\t' '{{print $xml_column_index}}' {BsMd} > {BsMdTemp}"
subprocess.run(awk_cmd, shell=True, check=True)


def parse_xml(xml_string):
    # Placeholder for the function that parses a single XML record
    # and returns a dictionary of the data
    pass

# Define the output CSV file for the processed data
processed_data_path = 'path/to/your/processed_data.csv'

# Open the output file in append mode
with open(processed_data_path, 'a') as output_file:
    # Read the XML data line by line or in small batches
    with open(output_xml_file_path, 'r') as xml_file:
        for xml_record in xml_file:
            # Parse the XML record
            parsed_data = parse_xml(xml_record)
            # Convert the parsed data to a DataFrame and write to the CSV file
            df = pd.DataFrame([parsed_data])
            df.to_csv(output_file, header=False, index=False)
