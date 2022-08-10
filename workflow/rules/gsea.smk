import glob

genesets = ["GeneOntologies", "Pathways"]
gsea_set_names = config['gsea_set_names']
gsea_set_paths = config['gsea_set_paths']
gsea_set_descriptions = config['gsea_set_descriptions']


# Wrap this in a checkpoint
def get_gsea_rnkfiles(wc):
    return glob.glob(deseqdir + f"/{wc.analysis}/rnk-files/*.rnk")


def get_gsea_gmtfile(wc):
    gene_set = gsea_set_names[wc.geneset]
    return gsea_set_paths[gene_set]


checkpoint rule gsea:
    input:
        gmt = get_gsea_gmtfile,
        desc = gsea_set_descriptions,
        rnks = get_gsea_rnkfiles,
        xls = deseqdir + "/{analysis}/analysis.xlsx"
    output:
        xls = gseadir + "/{analysis}/GSEA-{geneset}.xlsx"
    log:
        gseadir + "/logs/{analysis}/{geneset}.log"
    conda:
        "../envs/deseq2.yaml"
    resources:
        mem_mb = 8000
    threads: 10
    script:
        "../scripts/R/gsea.R"

rule run_gsea:
    input:
        expand(gseadir + "/{analysis}/GSEA-{geneset}.xlsx", analysis=analyses, geneset=genesets)


