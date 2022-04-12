strandedness = [config['strandedness'] for v in samples]

analyses = config['analyses']

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
        deseqdir + "/dds.rds",
        deseqdir + "/normalized-counts.tsv"
    log:
        deseqdir + "/logs/deseq-init.log"
    params:
        samples = config["samples"]
    conda:
        "../envs/deseq2.yaml"
    threads: 1
    script:
        "../scripts/deseq-init.R"

rule deseq_results:
    input:
        dds = deseqdir + "/dds.rds"
    output:
        dds = deseqdir + "/{analysis}/dds.rds",
        res = deseqdir + "/{analysis}/results.rds",
        xlsx = deseqdir + "/{analysis}/diffexp.xlsx"
    log:
        deseqdir + "/logs/{analysis}-deseq.log"
    params:
        model = lambda wc: analyses[wc.analysis]['model'],
        contrasts = lambda wc: analyses[wc.analysis]['contrast']
    conda:
        "../envs/deseq2.yaml"
    threads: 1
    script:
        "../scripts/deseq2-analysis.R"

rule run_deseq:
    input:
        expand(deseqdir + "/{analysis}/diffexp.xlsx", analysis=analyses)

