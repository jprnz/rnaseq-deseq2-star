# Changelog

## 2023-06-29
* Switched to use snakemake's slurm module
* Added default command-line options via --workflow-profiles 
* Refactored fastq-combine and fastp to enable job-grouping
* Added option to restrict heatmaps to use sample found in contrasts
* Added rules for delivery
* Added `setup` script

## 2023-04-26
* Sort DESeq results by p-value
* Truncate long names in GSEA enrichment plots 

## 2022-08-10
* Re-work of analysis.yaml and sample-sheet definitions (see [here](config/README.md) for more information)
* Set ranking of results tables / heatmaps to use shrunken LFC (as recommended in the DESeq2 documentation)
* Added ability to define models with interaction terms or utilize the LRT analysis
* Added ability to specify strandedness for each sample separately
* Re-organized `workflow/scripts` 
* Simplified DAG
