import glob

gsea_set_names = config['gsea_set_names']
gsea_set_paths = config['gsea_set_paths']
gsea_set_descriptions = config['gsea_set_descriptions']

pathways = ["GeneOntologies", "Pathways"]

# Wrap this in a checkpoint
def get_gsea_rnkfiles(wc):
    return glob.glob(deseqdir + f"/{wc.analysis}/rnk-files/*.rnk")

def get_gsea_gmtfile(wc):
    gene_set = gsea_set_names[wc.pathway]
    return gsea_set_paths[gene_set]

checkpoint rule gsea:
    input:
        gmt = get_gsea_gmtfile,
        desc = gsea_set_descriptions,
        rnks = get_gsea_rnkfiles,
    output:
        xls = gseadir + "/{analysis}/{analysis}-{pathway}.xlsx"
    log:
        gseadir + "/logs/{analysis}/{analysis}-{pathway}.log"
    conda:
        "../envs/deseq2.yaml"
    threads: 1
    script:
        "../scripts/gsea.R"

rule run_gsea:
    input:
        expand(gseadir + "/{analysis}/{analysis}-{pathway}.xlsx", analysis=analyses, pathway=pathways)


