rule fastp:
    input:
        r1 = ancient(fastqdir + "/{sample}_R1.fastq.gz"),
        r2 = ancient(fastqdir + "/{sample}_R2.fastq.gz")
    output:
        r1 = temp(fastpdir + "/{sample}_R1.fastq.gz"),
        r2 = temp(fastpdir + "/{sample}_R2.fastq.gz"),
        json_report = fastpdir + "/json_reports/{sample}.json",
        html_report = fastpdir + "/html_reports/{sample}.html"
    log:
        fastpdir + "/logs/{sample}.log"
    conda:
        "../envs/fastp.yaml"
    resources:
        mem_mb = 16000
    threads: 16
    shell:
        "(set -x; fastp "
        "  -i {input.r1} "
        "  -I {input.r2} "
        "  -o {output.r1} "
        "  -O {output.r2} "
        "  --json {output.json_report} "
        "  --html {output.html_report} "
        "  --thread {threads} "
        "  --detect_adapter_for_pe) &> {log}"

rule run_fastp:
    input:
        expand(fastpdir + "/json_reports/{sample}.json", sample=samples)


