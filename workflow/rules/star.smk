# Determine strandedness
strand_val = config['strandedness']
if strand_val in ['none', 'yes', 'reverse']:
    strandedness = [strand_val for v in samples]
elif os.path.exists(strandedness_file):
    strand = pd.read_csv(strand_val, index_col="sample", sep="\t", dtype='object')
    try:
        strandedness = [strand.loc[v]['strandedness'] for v in samples]
    except:
        raise ValueError(
            f"Check to make sure \'{strandedness_file}\' is formatted "
            "correctly and all samples are acounted for")
else:
    raise ValueError(
        "Check to make sure \'strandedness\' is formatted correctly in config.yaml")

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

rule count_matrix:
    input:
        expand(stardir + "/{sample}/{sample}_ReadsPerGene.out.tab", sample=samples)
    output:
        stardir + "/counts.csv",
    log:
        deseqdir + "/logs/count-matrix.log"
    params:
        samples=samples,
        strand=strandedness
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/count-matrix.py"

rule run_counts:
    input: stardir + "/counts.csv"

rule run_star:
    input:
        expand(stardir + "/{sample}.bam.bai", sample=samples),
        expand(stardir + "/{sample}.bam", sample=samples)

