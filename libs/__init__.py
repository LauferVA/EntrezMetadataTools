# libs/__init__.py

# Import specific functions or classes from each module
from .cleanMd import clean_metadata
from .dict2Tsv import dictionary_to_tsv
from .getIds import fetch_ids
from .Ids2Records.py import ids_to_records
from .mergeMd import merge_metadata
from .xml2Tsv import xml_to_tsv


# After this runs, can run from libs import clean_metadata, dictionary_to_tsv, fetch_ids from any part of project
def setup():
    """ Initialize some settings or perform additional checks here if needed """
    print("Libs package initialized.")

setup()


# 1) Consider defining a function or class here that utilizes these imports,
# 2) Do we want to provide a composite utility that combines several of these functionalities?  # (question for later after discussion with team)
