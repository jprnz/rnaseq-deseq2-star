def get_ensembl_fasta(download=False, path=reference_path, log=None, **args):
    from snakemake import shell
    from operator import itemgetter

    if not log:
        log = "/dev/null"

    datatype, species, build, release = itemgetter('datatype', 'species', 'build', 'release')(args)

    # Derived params
    branch = "grch37/" if release >= 81 and build == "GRCh37" else ""
    spec = f"{build}" if int(release) > 75 else f"{build}.{release}"
    species_cap = species.capitalize()

    # Define url
    url_base = f"ftp://ftp.ensembl.org/pub/{branch}release-{release}/fasta/{species}/{datatype}/{species_cap}.{spec}.{{suffix}}.gz"

    if datatype == "dna":
        # Check if primary_assembly is available
        suffix = "dna.primary_assembly.fa"
        url_test = url_base.format(suffix=suffix)
        try:
            shell("set -x; (curl --head -Ss -L {url_test}) &> {log}")
        except:
            suffix = "dna_sm.toplevel.fa"
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


    # Define paths
    local_path = f"{path}/{species}/{release}-{build}/{species_cap}.{spec}.{suffix}"
    url = url_base.format(suffix=suffix)

    if download:
        try:
            shell("set -x; (curl -S -L {url} | zcat > {local_path}) &> {log}")
        except:
            raise ValueError(
                    "Unable to download requested sequence data from Ensembl using urls:\n"
                    "{url}\nDid you check that this combination of species, build, "
                    "and release is actually provided?".format(url=url))
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
                "Unable to download requested sequence data from Ensembl using url:\n"
                "{url}\n Did you check that this combination of species, build, "
                "and release is actually provided?".format(url=url))
    else:
        return local_path


genome_fasta = get_ensembl_fasta(datatype="dna", **reference_args)
genome_gtf = get_ensembl_gtf(**reference_args)

# Regular expression to parse gtf to genes file (gene_id, symbol, biotype)
gtf_genes_sed = 's/.*gene_id "([^"]+).* gene_name "([^"]+).* gene_biotype "([^"]+).*/\1,\2,\3/'


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


rule make_genes_file:
    input:
        genome_gtf
    output:
        f"{genome_gtf}.genes"
    log:
        reference_path + "/logs/{species}-{release}-{build}-gtf-genes.log".format(**reference_args)
    cache: True
    script:
        "../scripts/gtf_symbols.py"


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


