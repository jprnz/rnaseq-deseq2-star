# General settings
To configure this workflow, modify `config/config.yaml` according to your needs paying special attention to the refrence genome, output folders, and strandedness options. 

# Samples for preprocessing
FASTQ files should be listed in `fastqs.tsv` file. For each sample, the columns `sample`, `pair`, and `path` must be defined. `sample` indicates a unique unit that will be processed, `pair` indicates the sequence pair of the FASTQ file, and `path` should describe the location of the FASTQ file (as an absolute path or relative to the project directory). Multiple values for the same `sample` and `pair` will merged together before processing (useful when a sample is sequenced across lanes or if FASTQ files are split into separate parts).     

# Differential expression
Each analysis that will be performed is defined in `analysis.yaml`. Here, `analyses` is a list that can contain any number of models to analyse. Values for `model` are used specify the design for DESeq, and `contrasts` indicate which conditions to include in the analysis. A special condition `all` may be used to facilitate the analysis of every pairwise comparison for the condition of interest. See the DESeq2 [manual](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html) for more information about setting up a design and using contrasts.

In addition to `samples`, every variable in your model will need to have a corresponding column in `samples.tsv`.

An example set of configuration files can be found in the `.test/config` directory.
