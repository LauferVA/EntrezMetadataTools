# EntrezMetadataTools

ObtainEntrezMetadata.py3:
Capabilities:
- Obtain Entrez Metadata based on an arbitrary query. 
- Works irrespective of target Db (Bioproject, Biosample, Sra, etc).
- Built for high-throughput via batched queries. 
- Will pull IDs and metadata for 3M records in >2hrs with low overhead, **without parallelization.**
- Parallelization was added 2/27/2024 using `concurrent.futures`  `ThreadPoolExecutor` and `as_completed`. Script will now pull all Rna-Seq metadata from SRA in ~30 minutes

3 tasks remain before versioning as 1.0:
1. Add XML parsing (v0.5)
2. Add Data Integration (using ProcessEntrezMetadata.py3) for any number of queries run using ObtainEntrezMetadata.py3 (v0.9)
3. Add prioritization script to be able to deliver to Matt rank-ordered studies that are recommended.

Once the script provides helpful prioritization of studies, version as 1.0.

  
