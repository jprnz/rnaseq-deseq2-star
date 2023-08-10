# Get all analyses
analyses = {v['analysis']: dict(v) for v in config["analyses"]}


rule deseq:
    input:
        samples = lambda wc: analyses[wc.analysis]['samples'],
        config = ancient("config/analysis.yaml"),
        counts = stardir + "/counts.csv",
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
    input: 
        expand(deseqdir + "/{analysis}/analysis.xlsx", analysis=analyses)

