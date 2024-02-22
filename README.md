# EntrezMetadataTools

ObtainEntrezMetadata.py3:
Obtain Entrez Metadata based on an arbitrary query. 
Works irrespective of target Db (Bioproject, Biosample, Sra, etc).

Built for high-throughput via batched queries. 
Will pull IDs and metadata for 3M records in >2hrs with low overhead, **serially.**

At present there are 5 planned improvements to these functions, see the Issues posted on 2/22/2024.

ProcessEntrezMetadata.py3: 
Explodes the xml into an expanded tab delimited format.


