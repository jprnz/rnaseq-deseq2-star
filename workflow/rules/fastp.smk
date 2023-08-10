def get_fastqs(wc):
    paths = files.loc[(wc.sample, wc.pair)].path
    try:
        ret = sorted(paths.tolist())
    except:
        ret = [paths]
    return ret

rule fastq_combine:
    input:
        ancient(get_fastqs)
    output:
        temp(fastpdir + "/fastqs/{sample}_{pair}.fastq.gz")
    log:
        fastpdir + "/logs/combine-{sample}_{pair}.log"
    resources:
        mem_mb = 1000
    group: "fastp"
    run:
        if len(input) > 1:
            cmd = "cat {input} > {output}"
        else:
            cmd = "ln -v {input} {output}"
        shell("echo -e \"Running:\n" + cmd + "\" > {log}")
        shell("(set -x; " + cmd + ") &>> {log}")

rule fastp:
    input:
        r1 = ancient(fastpdir + "/fastqs/{sample}_R1.fastq.gz"),
        r2 = ancient(fastpdir + "/fastqs/{sample}_R2.fastq.gz")
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
    group: "fastp"
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

rule run_fastq_combine:
  input:
      expand(fastpdir + "/fastqs/{sample}_{pair}.fastq.gz", sample=samples, pair=pairs)

rule run_fastp:
    input:
        expand(fastpdir + "/json_reports/{sample}.json", sample=samples)


