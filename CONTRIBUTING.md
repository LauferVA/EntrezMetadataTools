
First to lay out the 5 core functions:
1. rapid, optimized metadata download using efetch, esummary, etc as implemented in Bio.Entrez
2. implementation of recursive, robust xml parsing to handle the variety of xml structures used by entrez to export them as flat files
3. data characterization and cleaning functions for the (now tab delimited) metadata that came from xml
4. once data are cleaned, merge functions to join entrez metadata from multiple dbs to build more and more complete records
5. once data are merged, interrogation scripts for use in prioritizing studies

Contributions and improvements, bug reports, etc are appreciated. In addition, plans for future functionality include:
1. recasting the functionality described above as a python class
2. create parallel functionality using entrez cloud resources.
