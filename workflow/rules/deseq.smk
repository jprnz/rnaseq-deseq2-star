# Get all analyses
analyses = {v['analysis']: dict(v) for v in config["analyses"]}

# Determine strandedness
strandedness_file = "config/strandedness.tsv"
if os.path.exists(strandedness_file):
    strand = pd.read_csv(strandedness_file, index_col="sample", sep="\t", dtype='object')
    try:
        strandedness = [strand.loc[v]['strandedness'] for v in samples]
    except:
        raise ValueError(
            f"Check to make sure \'{strandedness_file}\' is formatted "
            "correctly and all samples are acounted for")
else:
    strandedness = [config['strandedness'] for v in samples]


rule count_matrix:
    input:
        expand(stardir + "/{sample}/{sample}_ReadsPerGene.out.tab", sample=samples)
    output:
        deseqdir + "/counts.csv",
    log:
        deseqdir + "/logs/count-matrix.log"
    params:
        samples=samples,
        strand=strandedness
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/count-matrix.py"


rule deseq:
    input:
        samples = lambda wc: analyses[wc.analysis]['samples'],
        config = "config/analysis.yaml",
        counts = deseqdir + "/counts.csv",
        genes = genome_gtf + ".genes"
    output:
        dat = deseqdir + "/{analysis}/analysis.RData",
        xls = deseqdir + "/{analysis}/analysis.xlsx",
        cnt = deseqdir + "/{analysis}/counts.csv",
        norm = deseqdir + "/{analysis}/normalized-counts.csv"
    log:
        deseqdir + "/logs/{analysis}.log"
    params:
        analysis = lambda wc: analyses[wc.analysis]
    conda:
        "../envs/deseq2.yaml"
    resources:
        mem_mb = 8000
    script:
        "../scripts/R/deseq.R"


rule run_deseq:
    input: expand(deseqdir + "/{analysis}/analysis.xlsx", analysis=list(analyses.keys()))

