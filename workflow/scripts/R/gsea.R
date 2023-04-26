log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

gmtfile = snakemake@input[['gmt']]
descfile = snakemake@input[['desc']]
rnkfiles = snakemake@input[['rnks']]
xls_filename = snakemake@output[['xls']]
geneset = snakemake@wildcards[['geneset']]


library("fgsea")
library("data.table")
library("ggplot2")
library("WriteXLS")
library("RColorBrewer")
library("parallel")


# Header for xlsx
header = c(
  "gene_set"="MSigDB gene-set (for more information on gene-sets, visit: https://www.gsea-msigdb.org/gsea/msigdb/search.jsp)",
  "description"="Brief description of gene-set",
  "pval"="Enrichment p-value",
  "padj"="Adjusted p-value (Benjamini and Hochberg, 1995)",
  "log2err"="Expected error for the standard deviation of the p-value logarithm",
  "ES"="Enrichment score (sign indicating direction of change)",
  "NES"="Normalized ES (normalized to the mean ES of random samples of the same size)",
  "size"="Number of genes in the gene-set",
  "leadingEdge"="List of genes driving the enrichment of the gene-set")

run_gsea = function(fn, gmt) {
  rnk <- fread(fn, header=FALSE, col.names=c("symbol", "stat"))

  # Format rnk
  rnk[, symbol := toupper(symbol)]
  stats = rnk$stat
  names(stats) = rnk$symbol

  # Run GSEA
  set.seed(1234)
  ret = fgsea(pathways=gmt, stats=stats, eps=1e-100, minSize=10, maxSize=1000, nproc=1)
  return(ret)
}

format_table <- function(gsea, desc) {
  # Assemble final table
  ret = lapply(gsea, function(dat) {
    r = data.table(pathway=dat$pathway)
    r = merge(r, desc)
    r = merge(r, dat)
    r = na.omit(r, 'pval')
    r[, leadingEdge := sapply(leadingEdge, function(x) paste0(unlist(x), collapse=", "))]
    setnames(r, "pathway", "gene_set")
    setorder(r, pval)

    return(r)
  })
  return(ret)
}

write_table = function(tab, geneset, tab_path) {
  for (n in names(tab)) {
    filename = file.path(tab_path, paste0(n, "-", geneset, ".csv.gz"))
    fwrite(tab[[n]], file=filename)
  }
}

write_excel <- function(tab, header, filename) {
  # Add header
  xls = list("Column Descriptions"=data.frame("Columns"=names(header), "Description"=header))
  xls = c(xls, tab)

  # Write xlsx
  WriteXLS(xls, ExcelFileName=filename, BoldHeaderRow=TRUE, FreezeRow=1)
}

truncate_name = function(sets, maxlen=30) {
  ind = which(nchar(sets) > maxlen)
  ret = substr(sets, 1, maxlen)
  ret[ind] = paste0(ret[ind], "...")

  suf = sub(".*_(UP|DN)$", "\\1", sets[ind]) 
  suf[which(suf == sets[ind])] = ""
  ret[ind] = paste0(ret[ind], suf)

  while(any(duplicated(ret))) {
    ret = truncate_name(sets, maxlen + 5)
  }

  return(ret)
}

plot_nes = function(gsea, filename, ntop=10) {
  gsea = gsea[order(pval),]
  dat = gsea[NES > 0, ][1:ntop,]
  dat = rbind(dat, gsea[NES < 0,][1:ntop,])
  dat = na.omit(dat, "pathway")
  setnames(dat, "pathway", "gene_set")
  setorder(dat, NES)

  # Truncate long names and turn to factor
  dat[, gene_set := truncate_name(gene_set)]
  dat[, gene_set := factor(gene_set, levels=unique(gene_set))]

  # Use spectral but replace midpoint with white
  colors = brewer.pal(11, "Spectral")
  colors[6] = "#FFFFFF"
  colors = colorRampPalette(colors)(101)

  # Plot
  plt <- ggplot(dat, aes(x=NES, y=gene_set, fill=NES)) +
    geom_col(show.legend=FALSE) +
    scale_fill_gradientn(colors=colors) +
    ylab("") +
    theme_minimal() +
    theme(legend.position="none") +
    theme(axis.text.y=element_text(size=7))
  ggsave(filename, plot=plt, width=5, height=4.5)
}

# Make place for output
gsea_path = dirname(xls_filename)
tab_path = file.path(gsea_path, "tables")
dir.create(tab_path, recursive=TRUE, showWarnings=FALSE)

# Read input
desc <- fread(descfile)
gmt <- gmtPathways(gmtfile)

# Do gsea for each contrast
names(rnkfiles) <- sub(".rnk", "", basename(rnkfiles))
gsea <- mclapply(rnkfiles, run_gsea, gmt, mc.cores=10)

# Write plot
for (n in names(gsea)) {
  plt_filename = file.path(gsea_path, paste0(n, "-", geneset, ".pdf"))
  plot_nes(gsea[[n]], plt_filename)
}

# Write results
tab = format_table(gsea, desc)
write_table(tab, geneset, tab_path)
write_excel(tab, header, xls_filename)

