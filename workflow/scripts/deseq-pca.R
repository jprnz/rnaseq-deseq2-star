log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

dds = snakemake@input[[1]]
out = snakemake@output[[1]]

# Use only named params
params = snakemake@params[which(names(snakemake@params) != "")]

library("DESeq2")

plot_pca <- function(mat, dat, pcs=c("PC1", "PC2"), color_by=NULL, shape_by=NULL,
    samples=TRUE, size=NULL, alpha=NULL, color_legend=NULL, shape_legend=NULL,
    filename=NULL) {

  require(ggplot2)
  require(ggrepel)
  require(data.table)

  # Set defaults
  size <- ifelse(is.null(size), 4, size)
  alpha <- ifelse(is.null(alpha), 1, alpha)

  # Compute pca and add to colData
  pca <- prcomp(t(mat))
  pca_x <- data.table(pca$x, keep.rownames=TRUE)
  dat <- data.table(dat, keep.rownames=TRUE)
  plt <- merge(pca_x, dat)
  setnames(plt, "rn", "sample")

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
      aes(label=.data$sample),
      box.padding=0.5,
      min.segment.length=0,
      show.legend=FALSE)
  }

  # Save if requested
  if (!is.null(filename)) {
    ggsave(filename=filename, plot=ret, width=8, height=6.5)
  }

  # Return plot if requested
  invisible(ret)
}


# Load deseq2 data
dds <- readRDS(dds)
save.image()

# obtain normalized counts
mat <- assay(rlog(dds, blind=TRUE))
dat <- as.data.frame(colData(dds))
params[["filename"]] = out

# plot and write
do.call(plot_pca, c(list(mat, dat), params))
