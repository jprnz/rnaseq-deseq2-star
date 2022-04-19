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


get_coldata = function(filename) {
    ret <- read.table(
        samplefile,
        header=TRUE,
        row.names="sample",
        check.names=FALSE)

    # Use given row-order to determine factor levels
    for (col in colnames(ret)) {
        if (is.character(ret[, col])) {
            vals = ret[, col]
            ret[, col] <- factor(vals, levels=unique(vals))
        }
    }
    return(ret)
}


# Counts columns will be reordered to match sample data
coldat = get_coldata(samplefile)

# Get count data
cts <- read.table(
    countsfile,
    header=TRUE,
    row.names="gene_id",
    check.names=FALSE)
cts <- cts[, rownames(coldat)]

# Create DDS object
dds <- DESeqDataSetFromMatrix(
    countData=cts,
    colData=coldat,
    design= ~ 1)


# Remove uninformative genes
dds <- dds[rowSums(counts(dds)) >= 10, ]


# Normalization and preprocessing
dds <- DESeq(dds)


# Write dds object as RDS
saveRDS(dds, file=ddsfile)


# Write normalized counts
norm_counts = counts(dds, normalized=TRUE)
fwrite(data.frame("gene_id"=rownames(norm_counts), norm_counts), file=normfile)

