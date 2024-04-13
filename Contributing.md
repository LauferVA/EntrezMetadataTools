
Purpose:

1. rapid, optimized metadata download using efetch, esummary, etc as implemented in Bio.Entrez
2. implementation of recursive, robust xml parsing to handle the variety of xml structures used by entrez to export them as flat files
3. data characterization and cleaning functions for the (now tab delimited) metadata that came from xml
4. once data are cleaned, merge functions to join entrez metadata from multiple dbs to build more and more complete records
5. once data are merged, interrogation scripts for use in prioritizing studies


Contributions welcome!

0. Suggestions, improvements, bug reports, etc always appreciated. 

In addition, plans for future functionality include:

1. helping improve the existing class to be robust enough to be offered as a PyPI package.
2. create parallel functionality using entrez cloud resources.

help with those latter two goals would be very helpful
