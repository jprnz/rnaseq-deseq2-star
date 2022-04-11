def get_fastqs(wc):
    paths = files.loc[(wc.sample, wc.pair)].path
    try:
        ret = sorted(paths.tolist())
    except:
        ret = [paths]
    return ret

rule fastq_combine:
    input:
        get_fastqs
    output:
        temp(fastqdir + "/{sample}_{pair}.fastq.gz")
    log:
        fastqdir + "/logs/{sample}_{pair}.log"
    resources:
        mem_mb = 1000
    run:
        if len(input) > 1:
            cmd = "cat {input} > {output}"
        else:
            cmd = "ln -v {input} {output}"
        shell("echo -e \"Running:\n" + cmd + "\" > {log}")
        shell("(set -x; " + cmd + ") &>> {log}")

rule run_fastq_combine:
  input:
      expand(fastqdir + "/{sample}_{pair}.fastq.gz", sample=samples, pair=pairs)
