# EntrezMetadataTools

ObtainEntrezMetadata.py3:

Capabilities:
- Obtain Entrez Metadata based on an arbitrary query. 
- Works irrespective of target Db (Bioproject, Biosample, Sra, etc).
- Built for high-throughput via batched queries. 
- Will pull IDs and metadata for 3M records in >2hrs with low overhead, **without parallelization.**
- Parallelization was added 2/27/2024 using `concurrent.futures`  `ThreadPoolExecutor` and `as_completed`. Script will now pull all Rna-Seq metadata from SRA in ~30 minutes
- XML parsing capability added 3/21/2024.
    - XML parser capable of serial flattening of compound nested XML structures.
    - XML parser capable of handling broken / incorrect XML formatted files of certain kinds.

2 tasks remain before versioning as 1.0:
1. Add Data Integration across multiple queries run using ObtainEntrezMetadata.py3. This will be with IntegrateEntrezMetadata.py3.
2. Add prioritization script to be able to deliver to Matt rank-ordered studies that are recommended.

Once the script provides helpful prioritization of studies, version as 1.0.

  
