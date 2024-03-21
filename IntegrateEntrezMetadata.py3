
# the purpose of this script is to take the results objects from more than one query to ObtainEntrezMetadata.py3 and integrate them, anding completeness to as many records as possible.
# for instance, calls to sra, biosample, bioproject, and geo with regard to gene expression data in humans produces non-redundant datasets.

# integration of this data is aimed at 3 main goals:
  1) leverage non-redundant data streams to create more comprehensive metadata possible using 1 alone.
  2) streamline the process of manual study identification and curation / allow more rapid prioritization of candidate studies.
  3) identify / draw attention to studies that may not have been identified without an agnostic scan.

