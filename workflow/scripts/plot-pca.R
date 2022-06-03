log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

dds = snakemake@input[[1]]
out = snakemake@output[[1]]
params = snakemake@params

library("DESeq2")

plot_pca <- function(mat, dat, pcs=NULL, color_by=NULL, shape_by=NULL,
    samples=TRUE, size=NULL, alpha=NULL, color_legend=NULL, shape_legend=NULL,
    filename=NULL, ...) {

  require(ggplot2)
  require(ggrepel)

  # Set defaults
  pcs <- ifelse(is.null(pcs), c("PC1", "PC2"), pcs)
  size <- ifelse(is.null(size), 4, size)
  alpha <- ifelse(is.null(alpha), 1, alpha)

  # Compute pca and add to colData
  pca <- prcomp(t(mat))

  # Get percent variance
  pctvar <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)))
  pctvar <- pctvar[which(colnames(pca$x) %in% pcs)]

  # Merge data, set sample_name, and remove "Row.names"
  plt <- merge(pca$x, dat, by="row.names")
  plt$samples <- plt[,'Row.names']
  rownames(plt) <- plt$samples
  plt$Row.name <- NULL

  # Try to preserve ordering
  plt <- plt[intersect(rownames(dat), rownames(plt)),]

  # Make base plot
  ret <- ggplot(plt, aes_string(x=pcs[1], y=pcs[2], color=color_by, shape=shape_by))
  ret <- ret + xlab(paste0(pcs[1], ": ", pctvar[1], "% variance"))
  ret <- ret + ylab(paste0(pcs[2], ": ", pctvar[2], "% variance"))
  ret <- ret + geom_point(size=size, alpha=alpha)

  if (samples) {
    ret <- ret + geom_label_repel(aes(colour=color_by))
  }

  ret <- ret + theme_light()


  # Save if requested
  if (!is.null(filename)) {
    ggsave(filename=filename, plot=ret)
    dev_off()
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
