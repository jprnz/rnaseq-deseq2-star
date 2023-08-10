# RNA-Seq Workflow
<div align="center">
    <img src="resources/dag.svg" width="40%" height="40%">
</div>

# Differential Expression using STAR / DESeq2
This workflow is designed to perform simple case-control analysis of one or more variables using DESeq.
The workflow leverages ties together the following pieces of software:
* [fastp](https://github.com/OpenGene/fastp)
* [STAR](https://github.com/alexdobin/STAR)
* [RSeQC](http://rseqc.sourceforge.net/)
* [MultiQC](https://multiqc.info/)
* [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)

All software dependencies are automatically resolved using [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/concepts/index.html).
Local execution and running via [SLURM](https://slurm.schedmd.com/) are supported.  
When running from an HPC all work should be done on an interactive session.

# Setup
## Code
From GitLab, fork this project filling out the project name, description, and ensure under namespace "SGT-Projects" is selected.  
Then, clone this repository into your project directory
```
git clone <url>
```

Checkout new analysis branch and push this to your project repository
```
git checkout -b analysis
git push -u origin analysis
```

Change your default branch to `analysis` in GitLab via setting -> repository -> default branch.

## Conda
Install conda and activate the snakemake environment
```
./setup
source conda/bin/activate snakemake
```

To run the example data using SLURM, from the `.test/` directory, run
```
cd .test
snakemake deliver_analysis -s ../workflow/Snakefile --slurm
```

# Running
See Snakemake documentation for more [command line options](https://snakemake.readthedocs.io/en/stable/executing/cli.html#all-options).
Default command line options are set using `--workflow-profiles` and can be found [here](workflow/profiles/config.yaml).
You can override these pre-sets by using their command-line equivalents or clear these using `--workflow-profiles none`

## Analysis
Setup genome and ensure the settings in `config/config.yaml` are correct, then prepare the workflow configuration files `fastqs.tsv`, `samples.tsv`, `analysis.yaml`.    
See documentation [here](config/README.md) for more information.

To run the analysis workflow on SLURM (from a interactive compute node): 
```
snakemake --slurm
```

Then, to stage an analysis for delivery, run:
```
snakemake deliver_analysis
```

To run only QC:
```
snakemake deliver_qc
```

Or, just to deiver a counts table and BAM files
```
snakemake deliver_counts
```

# Wrapping up

## Pushing your changes
Once an analysis is complete, ensure all your changes are commited and pushed to your project's repository.

Check which files need updating or added to the repository
```
git status
```

Stage existing files to be commited, add any addtional files that are needed, then commit those changes including a little description
```
git add .
git commit -a -m "Analysis as of commit date"
```

Update the project's repository
```
git push -u origin analysis
```

## Finalizing things
Once an analysis is complete, ensure that
- [ ] `Methods.docx` in the delivery folder is accurate and up-to-date
- [ ] The delivery folder is uploaded to Box
- [ ] All your changes are commited pushed to your project's repository
