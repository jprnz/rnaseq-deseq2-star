# logging
sys.stderr = open(snakemake.log[0], "w")

import sys

import pandas as pd

from collections import namedtuple

Sample = namedtuple("Sample", ['name', 'strand'])

def get_column(strand):
    if pd.isnull(strand) or strand == "none":
        return 1  # non stranded protocol
    elif strand == "yes":
        return 2  # 3rd column
    elif strand == "reverse":
        return 3  # 4th column, usually for Illumina truseq
    else:
        raise ValueError(
            (
                "'strandedness' column should be empty or have the "
                "value 'none', 'yes' or 'reverse', instead has the "
                "value {}"
            ).format(repr(strand))
        )

def read_file(file, sample):
    column = get_column(sample.strand)
    return pd.read_table(
        file,
        skiprows=4,
        index_col=0,
        header=None,
        names=["gene_id", sample.name],
        usecols=[0, column],
        dtype=object
    )


# Get all per-sample info
samples = [
    Sample(sample, strand)
    for sample, strand in zip(snakemake.params.samples, snakemake.params.strand)
]

# Read files
counts = [
    read_file(f, sample)
    for f, sample in zip(snakemake.input, samples)
]

# Combine
matrix = pd.concat(counts, axis=1)
matrix.index.name = "gene_id"

# Write
matrix.to_csv(snakemake.output[0], sep="\t")
