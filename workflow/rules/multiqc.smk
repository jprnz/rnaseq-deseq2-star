rule multiqc:
    input:
        rules.run_fastp.input,
        rules.run_star.input,
        rules.run_rseqc.input
    output:
        multiqcdir + "/QC.html"
    log:
        multiqcdir + "/logs/multiqc.log"
    conda:
        "../envs/multiqc.yaml"
    params:
        modules = " ".join([f"-m {v}" for v in ["fastp", "star", "rseqc"]]),
        input_path = [fastpdir, stardir, rseqcdir],
        output_path = multiqcdir,
        output_name = "QC.html",
        config = "config/multiqc.yaml"
    shell:
        "(set -x; multiqc -f "
        "-o {params.output_path} "
        "-n {params.output_name} "
        "-c {params.config} "
        "{params.modules} "
        "{params.input_path} "
        ") &> {log}"

rule run_multiqc:
    input: rules.multiqc.output
