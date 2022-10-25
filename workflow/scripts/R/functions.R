read_config <- function(filename, analysis_name) {
    dat = read_yaml(configfile)[['analyses']]
    found = FALSE
    for (ret in dat) {
        if (ret$analysis == analysis_name) {
            found = TRUE
            break()
        }
    }

    if (found == FALSE) {
        stop('Could not find configuration for analysis "', analysis_name, '"')
    }

    # Consistancy check
    if (is.null(ret$results$names) & is.null(ret$results$contrasts)) {
        stop("At least one value for 'names' or 'contrasts' must be specified under 'results'\n")
    }
    return(ret)
}

read_samplefile = function(filename) {
    ret = fread(filename, check.names=FALSE, data.table=FALSE)

    # Rename samples if requested
    if ("alias" %in% colnames(ret)) {
        rownames(ret) = ret$alias
    } else {
        rownames(ret) = ret$sample
    }

    # Convert strings to factors
    for (col in colnames(ret)) {
        ret[, col] = factor(ret[, col], levels=unique(ret[, col]))
    }

    return(ret)
}

read_counts = function(filename, coldat) {
    ret = fread(
        filename,
        select=c("gene_id", as.character(coldat$sample)),
        stringsAsFactors=FALSE,
        check.names=FALSE,
        data.table=FALSE)

    # Add rownames
    rownames(ret) = ret$gene_id
    ret$gene_id = NULL

    # Rename samples if requested
    if ("alias" %in% colnames(coldat)) {
        alias = coldat$alias
        names(alias) = coldat$sample
        colnames(ret) = alias[colnames(ret)]
    }

    return(ret)
}

get_all_vs_all = function(val, coldat) {
    all = combn(levels(coldat[, val[1]]), 2)
    all = apply(all, 2, rev)
    colnames(all) = apply(all, 2, paste0, collapse="-vs-")
    ret = apply(all, 2, function(x) c(val[1], x[1], x[2]), simplify=FALSE)
    return(ret)
}

parse_contrasts = function(dat, dds) {
    coldat = as.data.frame(colData(dds))
    vars = all.vars(design(dds))

    ret = list()
    for (val in dat) {
        # Test for first form of 'contrasts'
        if (val[1] %in% vars) {
            if (all(val[2:3] %in% levels(coldat[, val[1]]))) {
                # Both conditions specified
                name = paste0(val[2], "-vs-", val[3])
                ret[[name]] = val
            } else if (val[2] == "all") {
                # Determine all-vs-all comparisons
                ret = c(ret, get_all_vs_all(val, coldat))
            } else {
                # Nothing to do
                stop("Could not determine what comparison to make given:", val, "\n")
            }
        } else {
            # We have c("name", "coef", "coef")
            ret[[val[1]]] = val[2:3]
        }
    }
    return(ret)
}

parse_names = function(dat, dds) {
    rn = resultsNames(dds)
    ret = list()
    for (val in dat) {
        if (val[2] %in% rn) {
            ret[[val[1]]] = val[2]
        } else {
            stop("Result name \"", val[2], "\" not in resultsNames()\n")
        }
    }
    return(ret)
}

get_comparisons = function(res, dds) {
    contrasts = parse_contrasts(res$contrasts, dds)
    names = parse_names(res$names, dds)
    return(list("contrast"=contrasts, "name"=names))
}

get_results = function(dds, comparisons, config) {
    ret = list()
    cfg = config[setdiff(names(config), c("contrasts", "names"))]
    for (type in names(comparisons)) {
        for (cmp in names(comparisons[[type]])) {
            val = list(comparisons[[type]][[cmp]])
            names(val) = type
            res = do.call(results, c(list(dds), val, cfg))
            ret[[cmp]] = lfcShrink(dds, res=res, type="ashr")
        }
    }
    return(ret)
}

get_group_counts = function(dds) {
    cnt = counts(dds, normalized=TRUE)
    dat = as.data.frame(colData(dds))
    ret_list = list()
    for (var in all.vars(design(dds))) {
        for (lvl in levels(dat[, var])) {
            ids = which(dat[, var] == lvl)
            if (length(ids) > 1) {
                ret_list[[lvl]] = rowMeans(cnt[, ids])
            } else {
                ret_list[[lvl]] = count[, ids]
            }
        }
    }
    ret = do.call(data.frame, ret_list)
    names(ret) = paste0(names(ret_list), "_mean")
    ret$gene_id = rownames(ret)
    return(ret)
}

format_tab = function(res, dds, genes, cnt) {
    dat = as.data.frame(colData(dds))
    res = as.data.frame(res)
    res$gene_id = rownames(res)

    ret = data.table("gene_id"=rownames(res))
    ret = merge(ret, genes, all.x=TRUE)
    ret = merge(ret, res)
    ret = merge(ret, cnt)
    ret[, baseMean := NULL]
    ret = ret[order(abs(ret$log2FoldChange), decreasing=TRUE, na.last=TRUE),]
    return(ret)
}

format_rnk = function(tab) {
    if ("stat" %in% names(tab)) {
        col = "stat"
    } else {
        col = "log2FoldChange"
    }
    keep = which(!is.na(tab$pvalue))
    ret = tab[keep, c("symbol", col), with=FALSE]
    ret[, symbol := toupper(symbol)]
    setorderv(ret, col)
    return(ret)
}

write_xls = function(tab, cmps, filename) {
    ret = list()
    names = cmps$name
    contrasts = cmps$contrast

    # Look for contrasts / coefs
    for (cmp in names(tab)) {
        ret[[cmp]] = copy(tab[[cmp]])

        if (!is.null(contrasts[[cmp]])) {
            ctr = contrasts[[cmp]]
            lfc_name = paste0("logFC (", ctr[2], " / ", ctr[3], ")")
        } else if (!is.null(names[[cmp]])) {
            lfc_name = paste0("logFC (", names[[cmp]], ")")
        } else {
            lfc_name = "LogFC"
        }
        setnames(ret[[cmp]], "log2FoldChange", lfc_name)
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

write_counts = function(mat, filename) {
    ret = data.table("gene_id"=rownames(mat), mat)
    fwrite(ret, file=filename)
}

# Plotting
plot_pca = function(
        mat, dat, pcs=c("PC1", "PC2"), color_by=NULL, shape_by=NULL,
        samples=TRUE, size=NULL, alpha=NULL, color_legend=NULL, shape_legend=NULL,
        filename=NULL) {

    require(ggplot2)
    require(ggrepel)
    require(data.table)

    # Set defaults
    size <- ifelse(is.null(size), 4, size)
    alpha <- ifelse(is.null(alpha), 1, alpha)

    # rlog values
    mat = assay(rlog(dds))

    # Compute pca and add to colData
    pca <- prcomp(t(mat))
    pca_x <- data.table(pca$x, keep.rownames=TRUE)
    dat <- data.table(dat, keep.rownames=TRUE)
    plt <- merge(pca_x, dat)
    browser()

    # Get percent variance
    pctvar <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)))
    pctvar <- pctvar[which(colnames(pca$x) %in% pcs)]

    # Make base plot
    ret <- ggplot(plt, aes_string(x=pcs[1], y=pcs[2], color=color_by, shape=shape_by))
    ret <- ret + xlab(paste0(pcs[1], " (", pctvar[1], "% variance)"))
    ret <- ret + ylab(paste0(pcs[2], " (", pctvar[2], "% variance)"))
    ret <- ret + geom_point(size=size, alpha=alpha)
    ret <- ret + theme_light()

    if (samples) {
        ret <- ret + geom_text_repel(
            aes(label=.data$rn),
            box.padding=0.5,
            label.padding=0.5,
            min.segment.length=0,
            max.overlaps=20,
            show.legend=FALSE)
    }

    # Save if requested
    if (!is.null(filename)) {
        ggsave(filename=filename, plot=ret, width=7.5, height=6)
    }

    # Return plot if requested
    invisible(ret)
}

plot_heatmap = function(res, dds, filename, n=100) {
    # Subset ids to use (top 100, ranked by pvalue)
    dat = as.data.frame(colData(dds))[, all.vars(design(dds)), drop=FALSE]
    res = na.omit(as.data.frame(res))
    res = res[order(abs(res$log2FoldChange), decreasing=TRUE, na.last=TRUE),]
    ids = na.omit(rownames(res)[1:100])

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

