def get_ensembl_fasta(download=False, path=reference_path, log=None, **args):
    from snakemake import shell
    from operator import itemgetter

    if not log:
        log = "/dev/null"

    datatype, species, build, release = itemgetter('datatype', 'species', 'build', 'release')(args)

    if datatype == "dna":
        suffix = "dna.primary_assembly.fa"
    elif datatype == "cdna":
        suffix = "cdna.all.fa"
    elif datatype == "cds":
        suffix = "cds.all.fa"
    elif datatype == "ncrna":
        suffix = "ncrna.fa"
    elif datatype == "pep":
        suffix = "pep.all.fa"
    else:
        raise ValueError("invalid datatype, must be one of dna, cdna, cds, ncrna, pep")

    branch = "grch37/" if release >= 81 and build == "GRCh37" else ""
    spec = f"{build}" if int(release) > 75 else f"{build}.{release}"
    species_cap = species.capitalize()
    local_path = f"{path}/{species}/{release}-{build}/{species_cap}.{spec}.{suffix}"

    if download:
        url = f"ftp://ftp.ensembl.org/pub/{branch}release-{release}/fasta/{species}/{datatype}/{species_cap}.{spec}.{suffix}.gz"
        try:
            shell("set -x; (curl -S -L {url} | zcat > {local_path}) &> {log}")
        except Exception as e:
            raise ValueError(
                "Unable to download requested sequence data from Ensembl. "
                "Did you check that this combination of species, build, "
                "and release is actually provided?")
    else:
        return local_path


def get_ensembl_gtf(download=False, path=reference_path, log=None, **args):
    from snakemake import shell
    from operator import itemgetter

    if not log:

        log = "/dev/null"

    species, build, release = itemgetter('species', 'build', 'release')(args)

    branch = "grch37/" if release >= 81 and build == "GRCh37" else ""
    spec = f"{build}" if int(release) > 75 else f"{build}.{release}"
    species_cap = species.capitalize()
    local_path = f"{path}/{species}/{release}-{build}/{species_cap}.{spec}.gtf"

    if download:
        url = f"ftp://ftp.ensembl.org/pub/{branch}release-{release}/gtf/{species}/{species_cap}.{build}.{release}.gtf.gz"
        try:
            shell("set -x; (curl -L {url} | zcat > {local_path}) &> {log}")
        except Exception as e:
            raise ValueError(
                "Unable to download requested sequence data from Ensembl. "
                "Did you check that this combination of species, build, "
                "and release is actually provided?")
    else:
        return local_path

genome_fasta = get_ensembl_fasta(datatype="dna", **reference_args)
genome_gtf = get_ensembl_gtf(**reference_args)

rule download_genome:
    output:
        genome_fasta
    log:
        reference_path + "/logs/{species}-{release}-{build}-dna-fasta.log".format(**reference_args)
    cache: True
    run:
        get_ensembl_fasta(datatype='dna', download=True, log=log[0], **reference_args)


rule download_annotation:
    output:
        genome_gtf
    log:
        reference_path + "/logs/{species}-{release}-{build}-gtf.log".format(**reference_args)
    cache: True
    run:
        get_ensembl_gtf(download=True, log=log[0], **reference_args)


rule genome_faidx:
    input:
        genome_fasta
    output:
        f"{genome_fasta}.fai"
    log:
        reference_path + "/logs/{species}-{release}-{build}-dna-fai.log".format(**reference_args)
    conda:
        "../envs/samtools.yaml"
    cache: True
    shell:
        "(set -x; samtools faidx {input}) &> {log}"


rule star_index_genome:
    input:
        fasta = genome_fasta,
        gtf = genome_gtf
    output:
        directory(star_index)
    log:
        star_index + ".log"
    conda:
        "../envs/star.yaml"
    resources:
        mem_mb = 32000
    cache: True
    threads: 16
    shell:
        "(set -x; STAR"
        "  --genomeDir {output} "
        "  --runMode genomeGenerate "
        "  --genomeFastaFiles {input.fasta} "
        "  --sjdbGTFfile {input.gtf} "
        "  --sjdbOverhang 100 "
        "  --runThreadN {threads} "
        ") &> {log}"


