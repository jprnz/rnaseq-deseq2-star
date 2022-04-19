#log <- file(snakemake@log[[1]], open="wt")
#sink(log)
#sink(log, type="message")


threads <- snakemake@threads
genesfile <- snakemake@input[['genes']]
ddsfile <- snakemake@input[['dds']]
configfile <- snakemake@input[['config']]
analysis <- snakemake@params[['analysis']]
outdat <- snakemake@output[['dat']]
outxls <- snakemake@output[['xls']]


library("DESeq2")
library("WriteXLS")
library("data.table")
library("pheatmap")
library("yaml")

save.image()

# Parse contrasts
get_contrasts = function(dat, dds) {
    ret = list()
    coldat = as.data.frame(colData(dds))
    for (ctr in dat) {
        if (length(ctr) == 3) {
            n = paste0(ctr[2], "-vs-", ctr[3])
            ret[[n]] = ctr
        } else if (ctr[2] == "all") {
            all = combn(levels(coldat[, ctr[1]]), 2)
            all = apply(all, 2, rev)
            colnames(all) = apply(all, 2, paste0, collapse="-vs-")
            r = apply(all, 2, function(x) c(ctr[1], x[1], x[2]), simplify=FALSE)
            ret = c(ret, r)
        } else {
            n = paste(ctr, collapse="-")
            ret[[n]] = ctr
        }
    }
    return(ret)
}

get_sample_counts = function(cntr, dat, count) {
    ret_list = list()
    col = cntr[1]
    for (lvl in cntr[2:3]) {
        ids = rownames(dat)[dat[, col] == lvl]
        ret_list[[lvl]] = rowMeans(count[, ids])
    }
    names(ret_list) = paste0(names(ret_list), "_mean")
    ret = do.call(data.frame, ret_list)
    ret$gene_id = rownames(ret)
    return(ret)
}

format_tab = function(res, cntr, dds, genes) {
    dat = as.data.frame(colData(dds))

    ret = data.table("gene_id"=rownames(res))
    ret = merge(ret, genes, all.x=TRUE)
    ret = cbind(ret, as.data.table(res))
    ret[, baseMean := NULL]

    if (cntr[1] %in% names(dat)) {
        # Get mean counts
        count = counts(dds, normalized=TRUE)
        sample_counts = get_sample_counts(cntr, dat, count)
        ret = merge(ret, sample_counts)
    }
    return(ret)
}

format_rnk = function(res) {
    if ("stat" %in% names(res)) {
        col = "stat"
    } else {
        col = "log2FoldChange"
    }
    ret = data.table(gene_id=rownames(res), stat=res[, col])
    setorder(ret, stat)
    return(ret)
}


plot_heatmap = function(res, cntr, dds, filename, n=100) {
    # Subset ids to use (top 100, ranked by pvalue)
    dat = as.data.frame(colData(dds))[, cntr[1], drop=FALSE]
    res = as.data.frame(res)
    setorder(res, pvalue)
    ids = rownames(res)[1:n]

    # Z-score rlog counts
    mat = assay(rlog(dds))[ids,]
    mat = t(scale(t(mat)))

    # zlim
    lim = max(abs(range(mat)))
    browser()
    breaks = seq(-lim, lim, length.out=101)

    # Use spectral but replace midpoint with white
    colors = brewer.pal(11, "Spectral")
    colors[6] = "#FFFFFF"
    colors = colorRampPalette(colors)(101)

    # Plot
    pheatmap(
        mat, colors, annotation_col=dat, breaks=breaks, scale="none", border_color=NA,
        show_rownames=FALSE, show_colnames=TRUE, angle_col=45,
        cluster_rows=TRUE, clustering_distance_rows="correlation",
        cluster_cols=TRUE, clustering_distance_cols="correlation",
        treeheight_row=0, treeheight_col=30,
        cellwidth=20, filename=filename)
}

write_xls = function(res, contrasts, filename) {
    ret = list()
    for (cntr in names(contrasts)) {
        ret[[cntr]] = as.data.frame(res[[cntr]])
        if (length(contrasts[[cntr]]) == 3) {
            lfc_name = paste0("logFC (", contrasts[[cntr]][2], " / ", contrasts[[cntr]][3], ")")
            setnames(ret[[cntr]], "log2FoldChange", lfc_name)
        }
    }

    desc_vec = c(
        "gene_id"="Ensembl gene identifier",
        "symbol"="Gene symbol",
        "biotype"="Biological product of gene (see: https://useast.ensembl.org/info/genome/genebuild/biotypes.html)",
        "logFC (group-A / group-B)"="Shrunken estimates of log2 fold change between two groups (effect size, see: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#moreshrink)",
        "lfcSE"="Standard error of logFC column",
        "pvalue"="p-value of test for genes passing independent filtering (see: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#indfilt)",
        "padj"="p-value adjusted for multiple testing using method proposed by Benjamini, Y., and Hochberg, Y. (1995).",
        "<group-A>_mean"="Mean of normalized counts amoung samples associated with group-A",
        "<group-B>_mean"="Mean of normalized counts amoung samples associated with group-B"
    )
    desc = data.table("Columns"=names(desc_vec), "Description"=desc_vec)
    ret = c(list("Column Descriptions"=desc), ret)
    WriteXLS(ret, ExcelFileName=filename, AdjWidth=TRUE, BoldHeaderRow=TRUE, FreezeRow=1)
}

# Read dds and apply model
dds <- readRDS(ddsfile)

# Load config
config = read_yaml(configfile)[['analyses']]
for (cfg in config) {
    if (cfg$analysis == analysis) {
        break()
    }
}

# Load genes file
genes = fread(genesfile)

# Get model and contrast
model = as.formula(cfg$model)
contrasts = get_contrasts(cfg$contrasts, dds)

# Add design to dds
design(dds) = model
dds = DESeq(dds)

# Create places for supplimentary outputs
prefix = dirname(outdat)
rnkpath = file.path(prefix, "rnk-files")
tabpath = file.path(prefix, "tables")
dir.create(rnkpath, showWarnings=FALSE)
dir.create(tabpath, showWarnings=FALSE)

# Run test and do not shrink LFC
res <- lapply(contrasts, function(x) {
    r = results(dds, contrast=x)
    r = lfcShrink(dds, contrast=x, res=r, type="ashr")
    return(r)
})

# Make output tables and rnk files
tab = list()
rnk = list()
for (cntr in names(contrasts)) {
    tab[[cntr]] = format_tab(res[[cntr]], contrasts[[cntr]], dds, genes)
    rnk[[cntr]] = format_rnk(res[[cntr]])
}

# Write results
for (cntr in names(contrasts)) {
    tab_filename = paste0(tabpath, "/", cntr, ".csv")
    fwrite(tab[[cntr]], file=tab_filename)

    rnk_filename = paste0(rnkpath, "/", cntr, ".rnk")
    fwrite(rnk[[cntr]], file=rnk_filename, sep="\t", col.names=FALSE)

    plt_filename = paste0(prefix, "/Heatmap-", cntr, ".pdf")
    plot_heatmap(res[[cntr]], contrasts[[cntr]], dds, filename)
}

# Write Excel
write_xls(res, contrasts, filename=outxls)

# Write workspace
save.image(outdat)


