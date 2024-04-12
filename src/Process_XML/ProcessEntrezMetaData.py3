import pandas as pd
from bs4 import BeautifulSoup
from collections import defaultdict
import os
import csv

def xml_to_dict(xml_string):
    """ convert XML to dict using beautiful soup. name dict keys using the XML element and subelement 
    names in order to preserve knowledge of incoming XML """
    soup = BeautifulSoup(xml_string, 'html.parser')  # lxml is default for bs4, but this selection more likely to work for varety of users (html.parser as an alternative parser).
    xml_dict = defaultdict(list)
    for tag in soup.find_all(recursive=False): # iterate over tags
        if tag.string and tag.string.strip():
            xml_dict[tag.name].append(tag.string.strip()) # extract the tags
        for attr_name, attr_value in tag.attrs.items():
            xml_dict[f"{tag.name}@{attr_name}"].append(attr_value) # extract the values of each tag.
    return xml_dict

def process_xml_data(xml_row_data):
    """ handles cases in which the xml data are poorly formed and breaks standard xml parsers an initial attempt to parse the 
    xml is made. if this fails, the index where the non-standard formatting (that caused the error) is recorded then used to 
    isolate a manageable section of the xml. this process is repeated until the xml data are fully parsed. """
    final_xml_df = pd.DataFrame()
    xml_data = xml_row_data
    while xml_data:
        try:
            xml_dict = xml_to_dict(xml_data)
            final_xml_df = pd.concat([final_xml_df, pd.DataFrame(xml_dict)], ignore_index=True)
            break  # Process only the first chunk of XML data
        except Exception as e:
            error_message = str(e)
            print(f"XML parsing error encountered: {error_message}")
            # Find the position of the parsing error and truncate XML data
            error_index = int(error_message.split(':')[-1].strip())
            xml_data = xml_data[:error_index]
    """ I recall adding recursion to this - need to enable serial flattening of nested structures using recursion. """
    return final_xml_df

def process_line_with_xml(row_data):
    """ shunts various parts of each row of the data being read in into a different processing routine """
    row_combined_data = {}
    non_xml_fields = [key for key, value in row_data.items() if not key.endswith('Xml')]

    for key in non_xml_fields:
        row_combined_data[key] = row_data[key]     # Copy non-XML fields
    xml_fields = [key for key in row_data if key.endswith('Xml')]
    for xml_field in xml_fields:
        xml_data = row_data.get(xml_field, '')     # Process XML fields
        if xml_data:
            # Process XML column
            final_xml_df = process_xml_data(xml_data)
            if not final_xml_df.empty:
                # Merge flattened XML data into combined data
                row_combined_data.update(final_xml_df.iloc[0])
    return row_combined_data

def process_metadata_file_with_xml(file_path):
    """ at present, this script takes flat files created by ObtainEntrezMetadata.py3 as input via this function in the 
    future, will process the xml at the same time as Md are obtained from Entrez, so this function will be unnecessary. """
    combined_data = []
    with open(file_path, mode='r', encoding='utf-8') as file:
        reader = csv.DictReader(file, delimiter=",")
        for row in reader:
            row_combined_data = process_line_with_xml(row)         # Process each line with embedded XML
            combined_data.append(row_combined_data)
    final_df = pd.DataFrame(combined_data)                         # Convert combined data to DataFrame
    return final_df

def process_all_files_in_folder(folder_path):
    """ get all .tsv fnames in target dir; serially process stored metadata-containg flatfiles. as with
    process_metadata_file_with_xml(), will be phased out once this script is integrated with ObtainEntrezMetadata.py3. """
    dfs = []
    for filename in os.listdir(folder_path):
        if filename.endswith('.tsv'):
            file_path = os.path.join(folder_path, filename)
            df = process_metadata_file_with_xml(file_path)
            dfs.append(df)
    combined_df = pd.concat(dfs, ignore_index=True)
    return combined_df

# actually set target file / folder paths. replace this with command line arguments or draw the info from flat file.
folder_path = '../Data/'
final_combined_df = process_all_files_in_folder(folder_path)
del final_combined_df['Item']
final_combined_df.to_excel('../Out/combined_sra_rnaseq_metadata.xlsx', index=False) # Save the combined DataFrame as an Excel file
final_combined_df.to_csv(file_path, sep='\t', index=False) # can be helpful to save as flatfile instead ... example: Excel chokes after > 2^20 rows.

