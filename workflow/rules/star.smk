rule star:
    input:
        r1 = fastpdir + "/{sample}_R1.fastq.gz",
        r2 = fastpdir + "/{sample}_R2.fastq.gz",
        index = star_index,
        gtf = genome_gtf
    output:
        stardir + "/{sample}/{sample}_Aligned.sortedByCoord.out.bam",
        stardir + "/{sample}/{sample}_ReadsPerGene.out.tab"
    log:
        stardir + "/logs/{sample}.log"
    params:
        prefix = stardir + "/{sample}/{sample}_",
    conda:
        "../envs/star.yaml"
    resources:
        mem_mb = 48000
    group:
        "star"
    threads: 16
    shell:
        "(set -x; "
        "  STAR "
        "   --genomeDir {input.index}"
        "   --sjdbGTFfile {input.gtf}"
        "   --readFilesIn {input.r1} {input.r2}"
        "   --outFileNamePrefix {params.prefix}"
        "   --runThreadN {threads}"
        "   --outSAMtype BAM SortedByCoordinate "
        "   --quantMode GeneCounts "
        "   --alignSJoverhangMin 999 "
        "   --outFilterMultimapNmax 1 "
        "   --outSAMstrandField intronMotif "
        "   --outSAMattributes Standard "
        "   --readFilesCommand zcat "
        "   --outStd Log "
        ") &> {log}"


rule star_link:
    input:
        stardir + "/{sample}/{sample}_Aligned.sortedByCoord.out.bam"
    output:
        stardir + "/{sample}.bam"
    params:
        link = "{sample}/{sample}_Aligned.sortedByCoord.out.bam",
        name = "{sample}.bam",
        path = stardir
    log:
        stardir + "/logs/{sample}_link.log"
    group:
        "star"
    shell:
        "(set -x; cd {params.path}; ln -vs {params.link} {params.name}) &> {log}"


rule star_index_bam:
    input:
        stardir + "/{sample}.bam"
    output:
        stardir + "/{sample}.bam.bai"
    log:
        stardir + "/logs/{sample}_index.log"
    conda:
        "../envs/star.yaml"
    group:
        "star"
    threads: 5
    shell:
        "(set -x; samtools index -\@ 5 {input}) &> {log}"


rule run_star:
    input:
        expand(stardir + "/{sample}/{sample}_ReadsPerGene.out.tab", sample=samples) +
        expand(stardir + "/{sample}.bam.bai", sample=samples)

