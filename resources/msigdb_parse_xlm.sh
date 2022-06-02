#!/bin/bash
set -eo pipefail

# Example:
#   STANDARD_NAME="IBRAHIM_NRF1_UP
#   PMID="33096892
#   DESCRIPTION_BRIEF="Genes up-regulated in HEK293T cells overexpressing FLAG-NRF1
#   DESCRIPTION_FULL="The NRF transcription factors NRF1, NRF2, and NRF3, are a subset of Cap'n'collar transcriptional regulators which modulate the expression of genes harboring antioxidant-response element (ARE) sequences within their genomic loci. Despite the emerging physiological importance of NRF family members, the repertoire of their genetic targets remains incompletely defined. Here we use RNA-sequencing-based transcriptional profiling and quantitative proteomics to delineate the overlapping and differential genetic programs effected by the three NRF transcription factors. Comparing our data to recent profiling analyses, we create consensus target gene sets regulated by NRF1, NRF2, and NRF3, genetic programs which we determine to be differentially regulated in human tissues. Together, our data provide a quantitative assessment of how NRF family members sculpt proteomes and transcriptomes, essential information for future studies evaluating the role of NRF factors in normal physiology and disease.


echo -e "pathway\tdescription"

grep "STANDARD_NAME" $1 \
  | sed -r 's/.*STANDARD_NAME="([^"]*)".*DESCRIPTION_BRIEF="([^"]*)".*/\1\t\2/'

