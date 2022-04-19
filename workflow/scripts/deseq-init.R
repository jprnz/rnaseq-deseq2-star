# Setup log
log <- file(snakemake@log[[1]], open="wt")
sink(log, type="output")
sink(log, type="message")


# Ensure we have all parameters
countsfile <- snakemake@input[["counts"]]
samplefile <- snakemake@params[["samples"]]
ddsfile <- snakemake@output[["dds"]]
normfile <- snakemake@output[["norm"]]
threads <- snakemake@threads

library("DESeq2")
library("data.table")

parallel <- FALSE
if (threads > 1) {
    library("BiocParallel")
    register(MulticoreParam(threads))
    parallel <- TRUE
}


# Load sample data
coldata <- read.table(
    samplefile,
    header=TRUE,
    row.names="sample",
    check.names=FALSE)


# Use given row-order to determine factor levels
for (column in colnames(coldata)) {
    if (is.character(coldata[, column])) {
        vals = coldata[, column]
        coldata[, column] <- factor(vals, levels=unique(vals))
    }
}
str(coldata)

# Counts columns will be reordered to match sample data
cts <- read.table(
    countsfile,
    header=TRUE,
    row.names="gene_id",
    check.names=FALSE)
cts <- cts[, rownames(coldata)]

# Create DDS object
dds <- DESeqDataSetFromMatrix(
    countData=cts,
    colData=coldata,
    design= ~ 1)


# Remove uninformative genes
dds <- dds[ rowSums(counts(dds)) >= 10, ]


# Normalization and preprocessing
dds <- DESeq(dds, parallel=parallel)


# Write dds object as RDS
saveRDS(dds, file=ddsfile)


# Write normalized counts
norm_counts = counts(dds, normalized=TRUE)
fwrite(data.frame("gene_id"=rownames(norm_counts), norm_counts), file=normfile)

