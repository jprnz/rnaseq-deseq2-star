analyses = {v['analysis']: dict(v) for v in config["analyses"]}
strandedness = [config['strandedness'] for v in samples]


rule count_matrix:
    input:
        expand(stardir + "/{sample}/{sample}_ReadsPerGene.out.tab", sample=samples)
    output:
        deseqdir + "/counts.tsv",
    log:
        deseqdir + "/logs/count-matrix.log"
    params:
        samples=samples,
        strand=strandedness
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/count-matrix.py"


rule deseq_dds:
    input:
        counts = deseqdir + "/counts.tsv"
    output:
        dds = deseqdir + "/dds.rds",
        norm = deseqdir + "/normalized-counts.tsv"
    log:
        deseqdir + "/logs/init.log"
    params:
        samples = config["samples"]
    conda:
        "../envs/deseq2.yaml"
    threads: 1
    script:
        "../scripts/deseq-init.R"


rule deseq_analysis:
    input:
        config = "config/analysis.yaml",
        dds = deseqdir + "/dds.rds",
        genes = genome_gtf + ".genes"
    output:
        dat = deseqdir + "/{analysis}/analysis.RData",
        xls = deseqdir + "/{analysis}/analysis.xlsx"
    log:
        deseqdir + "/logs/{analysis}.log"
    params:
        analysis = lambda wc: wc.analysis
    conda:
        "../envs/deseq2.yaml"
    threads: 1
    script:
        "../scripts/deseq-analysis.R"


rule run_deseq:
    input:
        expand(deseqdir + "/{analysis}/analysis.xlsx", analysis=analyses)

