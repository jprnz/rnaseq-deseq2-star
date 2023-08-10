# General settings
To configure this workflow, modify `config/config.yaml` according to your needs paying special attention to the refrence genome, output folders, and strandedness options.

# Samples for preprocessing
FASTQ files should be listed in `fastqs.tsv` file. For each sample, the columns `sample`, `pair`, and `path` must be defined. `sample` indicates a unique unit that will be processed, `pair` indicates the sequence pair of the FASTQ file, and `path` should describe the location of the FASTQ file (as an absolute path or relative to the project directory). Multiple values for the same `sample` and `pair` will merged together before processing (useful when a sample is sequenced across lanes or if FASTQ files are split into separate parts).

# Differential expression
## Strandedness
The `strandedness` option in `config.yaml` defines the strandedness for all samples in the project. However, per-sample strandedness can be defined in a seperate file. To do so, include a tab-separated table with the columns `sample` and `strandedness`. Values for this table can be `none`, `yes`, or `reverse` (see `config.yaml` for a description of each of these).

## Samples for analysis
Samples to be included in an analysis need to be defined in a tab-separated file that contain `sample` (as defined in `fastqs.tsv`) and at least one condition to test. Any configuration of samples can be included. An optional column `alias` can be included to modify the name of a sample. Additional columns will be converted to factors with levels ordered according to the input. When doing an all-vs-all comparison using `[<condtion>, all]` as a value for `contrasts`, the sample order will determine which comparisons will be made. Generally it is best to include WT / baseline samples first, followed by the most interesting conditions.

## Analysis configuration
See the DESeq2 [manual](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html) for more information about setting up a design and using contrasts.

Each analysis that will be performed is defined in `analysis.yaml`. Here, `analyses` is a list that can contain any number of analyses:

```yaml
analysese:
    - analysis: name of the analysis
        design: ~ analysis formula
        samplesheet: file containing sample information (ex. samples.tsv)
        results:
            contrasts:
                - [condition, all] and / or...
                - [condition, level, level] and / or...
                - [display name, coefficent, coefficent]
            name:
                - [display name, coefficent] 
            <parameter>: <value> 
        deseq:
            <parameter>: <value>
        pca:
            pcs: list of two PCs to plot
            # Following are optional 
            color_by: variable used to determine color
            shape_by: variable used to determine shape
            samples: include sample names in plot (true/false)
            size: size of points
            alpha: transparency of points, as a percent
            color_legend: name to use for color legend
            shape_legend: name to use for shape legend
        heatmap:
            only_contrasts: Only use samples defined in contrasts (true/false)

    # Addtional analyses can be included... 
    - analysis: name of second analysis
        samplesheet: ... 
        design: ...
        # etc ...
```
### Notes
* Multiple analyses can be defined within this document
* At least one value for `names` or `contrasts` under `results` needs to be specified. 
* Both `results` and `deseq` can be used to pass any arbitrary parameter to either `results()` or `DESeq()`. 

### Tips for LRT analysis and more complex analyses
* For LRT, use 'name' to determine which coefficient use for LogFC 
* Output of resultsNames() is printed in the log file for the analysis

