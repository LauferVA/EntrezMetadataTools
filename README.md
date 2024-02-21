# EntrezMetadataTools

ObtainEntrezMetadata.py3:
Obtain Entrez Metadata based on an arbitrary query. 
Works irrespective of target Db (Bioproject, Biosample, Sra, etc).

Built for high-throughput via batched queries. 
Will pull IDs and metadata for 3M records in ~ 1 hr with low overhead.

> improvements needed: enable parallel processing
> enable error catching (e.g. due to http 500) and fluidly continue
> re-combine separate records into one afterwards.

ProcessEntrezMetadata.py3: 
Explodes the xml into an expanded tab delimited format.


