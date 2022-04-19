log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")


threads <- snakemake@threads
genesfile <- snakemake@input[['genes']]
ddsfile <- snakemake@input[['dds']]
configfile <- snakemake@input[['config']]
analysis <- snakemake@params[['analysis']]
outdat <- snakemake@output[['dat']]
outxls <- snakemake@output[['xls']]

library("DESeq2")
library("RColorBrewer")
library("WriteXLS")
library("data.table")
library("pheatmap")
library("yaml")

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

get_sample_counts = function(ctr, dat, count) {
    ret_list = list()
    col = ctr[1]
    for (lvl in ctr[2:3]) {
        ids = rownames(dat)[dat[, col] == lvl]
        ret_list[[lvl]] = rowMeans(count[, ids])
    }
    ret = do.call(data.frame, ret_list)
    names(ret) = paste0(names(ret_list), "_mean")
    ret$gene_id = rownames(ret)
    return(ret)
}

format_tab = function(res, ctr, dds, genes) {
    dat = as.data.frame(colData(dds))
    res = as.data.frame(res)
    res$gene_id = rownames(res)

    ret = data.table("gene_id"=rownames(res))
    ret = merge(ret, genes, all.x=TRUE)
    ret = merge(ret, res)
    ret[, baseMean := NULL]

    if (ctr[1] %in% names(dat)) {
        # Get mean counts
        count = counts(dds, normalized=TRUE)
        sample_counts = get_sample_counts(ctr, dat, count)
        ret = merge(ret, sample_counts)
    }
    setorder(ret, pvalue, na.last=TRUE)
    return(ret)
}

format_rnk = function(res) {
    if ("stat" %in% names(res)) {
        col = "stat"
    } else {
        col = "log2FoldChange"
    }

    keep = which(!is.na(res$pvalue))
    ret = data.table(gene_id=rownames(res), stat=res[, col])
    ret = ret[keep,]
    setorder(ret, stat)
    return(ret)
}


plot_heatmap = function(res, ctr, dds, filename, n=100) {
    # Subset ids to use (top 100, ranked by pvalue)
    dat = as.data.frame(colData(dds))[, ctr[1], drop=FALSE]
    res = as.data.frame(res)
    setorder(res, pvalue)
    ids = rownames(res)[1:n]

    # Z-score rlog counts
    mat = assay(rlog(dds))[ids,]
    mat = t(scale(t(mat)))

    # zlim
    lim = max(abs(range(mat)))
    breaks = seq(-lim, lim, length.out=101)

    # Use spectral but replace midpoint with white
    colors = brewer.pal(11, "Spectral")
    colors[6] = "#FFFFFF"
    colors = colorRampPalette(colors)(101)

    # Plot
    pheatmap(
        mat, colors, annotation_col=dat, breaks=breaks, scale="none", border_color=NA,
        show_rownames=FALSE, show_colnames=TRUE, angle_col=90,
        cluster_rows=TRUE, clustering_distance_rows="correlation",
        cluster_cols=TRUE, clustering_distance_cols="correlation",
        treeheight_row=0, treeheight_col=30,
        cellwidth=20, filename=filename, legend_col="")
}

write_xls = function(tab, contrs, filename) {
    ret = list()
    for (ctr in names(contrs)) {
        ret[[ctr]] = copy(tab[[ctr]])
        if (length(contrs[[ctr]]) == 3) {
            lfc_name = paste0("logFC (", contrs[[ctr]][2], " / ", contrs[[ctr]][3], ")")
            setnames(ret[[ctr]], "log2FoldChange", lfc_name)
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
        "<group-A>_mean"="Mean of normalized counts among samples associated with group-A",
        "<group-B>_mean"="Mean of normalized counts among samples associated with group-B"
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

save.image()

# Get model and contrast
model = as.formula(cfg$model)
contrs = get_contrasts(cfg$contrasts, dds)

# Add design to dds
design(dds) = model
dds = DESeq(dds)

# Create places for supplementary outputs
prefix = dirname(outdat)
rnkpath = file.path(prefix, "rnk-files")
tabpath = file.path(prefix, "tables")
dir.create(rnkpath, recursive=TRUE, showWarnings=FALSE)
dir.create(tabpath, recursive=TRUE, showWarnings=FALSE)

# Run test and do not shrink LFC
res <- lapply(contrs, function(x) {
    r = results(dds, contrast=x)
    r = lfcShrink(dds, contrast=x, res=r, type="ashr")
    return(r)
})

# Make output tables and rnk files
tab = list()
rnk = list()
for (ctr in names(contrs)) {
    tab[[ctr]] = format_tab(res[[ctr]], contrs[[ctr]], dds, genes)
    rnk[[ctr]] = format_rnk(res[[ctr]])
}

# Write results
for (ctr in names(contrs)) {
    tab_filename = paste0(tabpath, "/", ctr, ".csv")
    fwrite(tab[[ctr]], file=tab_filename)

    rnk_filename = paste0(rnkpath, "/", ctr, ".rnk")
    fwrite(rnk[[ctr]], file=rnk_filename, sep="\t", col.names=FALSE)

    plt_filename = paste0(prefix, "/Heatmap-", ctr, ".pdf")
    plot_heatmap(res[[ctr]], contrs[[ctr]], dds, plt_filename)
}

# Write Excel
write_xls(tab, contrs, filename=outxls)

# Write workspace
save.image(outdat)


