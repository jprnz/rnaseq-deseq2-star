rule rseqc_make_bed:
    input:
        genome_gtf
    output:
        bed = rseqcdir + "/annotation.bed",
        db = temp(rseqcdir + "/annotation.db")
    log:
        rseqcdir + "/logs/rseqc_gtf2bed.log"
    conda:
        "../envs/gffutils.yaml"
    script:
        "../scripts/rseqc-bed.py"


rule rseqc_junction_annotation:
    input:
        bam = stardir + "/{sample}.bam",
        bed = rseqcdir + "/annotation.bed"
    output:
        rseqcdir + "/{sample}/{sample}.junction.bed",
    log:
        rseqcdir + "/logs/junction_annotation/{sample}.log",
    params:
        extra = "-q 255",  # STAR uses 255 as a score for unique mappers
        prefix = rseqcdir + "/{sample}/{sample}"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "(set -x; junction_annotation.py "
        "  -i {input.bam}  "
        "  -r {input.bed} "
        "  -o {params.prefix} "
        "  {params.extra} "
        ") &> {log}"


rule rseqc_junction_saturation:
    input:
        bam = stardir + "/{sample}.bam",
        bed = rseqcdir + "/annotation.bed"
    output:
        rseqcdir + "/{sample}/{sample}.junctionSaturation_plot.pdf",
    log:
        rseqcdir + "/logs/junction_saturation/{sample}.log",
    params:
        extra=r"-q 255",
        prefix = rseqcdir + "/{sample}/{sample}"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "(set -x; junction_saturation.py "
        "-i {input.bam} "
        "-r {input.bed} "
        "-o {params.prefix} "
        "{params.extra} "
        ") &> {log}"


rule rseqc_bam_stat:
    input:
        bam = stardir + "/{sample}.bam",
        bed = rseqcdir + "/annotation.bed"
    output:
        rseqcdir + "/{sample}/{sample}.stats.txt",
    log:
        rseqcdir + "/logs/bam_stats/{sample}.log",
    conda:
        "../envs/rseqc.yaml"
    shell:
        "(set -x; bam_stat.py -i {input} > {output}) &> {log}"


rule rseqc_infer_experiment:
    input:
        bam = stardir + "/{sample}.bam",
        bed = rseqcdir + "/annotation.bed"
    output:
        rseqcdir + "/{sample}/{sample}.infer_experiment.txt",
    log:
        rseqcdir + "/logs/infer_experiment/{sample}.log",
    conda:
        "../envs/rseqc.yaml"
    shell:
        "(set -x; infer_experiment.py -r {input.bed} -i {input.bam} > {output}) &> {log}"


rule rseqc_inner_distance:
    input:
        bam = stardir + "/{sample}.bam",
        bed = rseqcdir + "/annotation.bed"
    output:
        rseqcdir + "/{sample}/{sample}.inner_distance.txt"
    log:
        rseqcdir + "/logs/inner_distance/{sample}.log",
    params:
        extra=r"-q 255",
        prefix = rseqcdir + "/{sample}/{sample}"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "(set -x; inner_distance.py "
        "-i {input.bam} "
        "-r {input.bed} "
        "-o {params.prefix} "
        ") &> {log}"


rule rseqc_read_distribution:
    input:
        bam = stardir + "/{sample}.bam",
        bed = rseqcdir + "/annotation.bed"
    output:
        rseqcdir + "/{sample}/{sample}.readdistribution.txt"
    log:
        rseqcdir + "/logs/read_distribution/{sample}.log",
    params:
        extra=r"-q 255",
        prefix = rseqcdir + "/{sample}/{sample}"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "(set -x; read_distribution.py "
        "-i {input.bam} "
        "-r {input.bed} "
        " 1> {output}) 2> {log}"


rule rseqc_read_duplication:
    input:
        bam = stardir + "/{sample}.bam",
    output:
        rseqcdir + "/{sample}/{sample}.pos.DupRate.xls"
    log:
        rseqcdir + "/logs/read_duplication/{sample}.log",
    params:
        extra=r"-q 255",
        prefix = rseqcdir + "/{sample}/{sample}"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "(set -x; read_duplication.py "
        "-i {input} "
        "-o {params.prefix} "
        ") &> {log}"


rule rseqc_readgc:
    input:
        bam = stardir + "/{sample}.bam",
    output:
        rseqcdir + "/{sample}/{sample}.GC.xls"
    log:
        rseqcdir + "/logs/read_gc/{sample}.log",
    params:
        extra=r"-q 255",
        prefix = rseqcdir + "/{sample}/{sample}"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "(set -x; read_GC.py -i {input} -o {params.prefix}) &> {log}"


rule rseqc_gene_body_coverage:
    input:
        bam = stardir + "/{sample}.bam",
        bed = rseqcdir + "/annotation.bed"
    output:
        rseqcdir + "/{sample}/{sample}.geneBodyCoverage.txt"
    log:
        rseqcdir + "/logs/genebody_coverage/{sample}.log",
    params:
        extra=r"-q 255",
        prefix = rseqcdir + "/{sample}/{sample}"
    conda:
        "../envs/rseqc.yaml"
    resources:
        mem_mb = 16000
    shell:
        "(set -x; geneBody_coverage.py "
        "-i {input.bam} "
        "-r {input.bed} "
        "-o {params.prefix} "
        ") &> {log}"


rule run_rseqc:
    input:
        expand(rseqcdir + "/{sample}/{sample}.junctionSaturation_plot.pdf", sample=samples),
        expand(rseqcdir + "/{sample}/{sample}.stats.txt", sample=samples),
        expand(rseqcdir + "/{sample}/{sample}.readdistribution.txt", sample=samples),
        expand(rseqcdir + "/{sample}/{sample}.pos.DupRate.xls", sample=samples),
        expand(rseqcdir + "/{sample}/{sample}.geneBodyCoverage.txt", sample=samples)

