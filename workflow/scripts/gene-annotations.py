# logging
sys.stderr = open(snakemake.log[0], "w")

import sys
import re

ID_RE = r'gene_id "([^"]+)"'
NAME_RE = r'gene_name "([^"]+)"'
BIOTYPE_RE = r'gene_biotype "([^"]+)"'
HEADER = "gene_id\tsymbol\tbiotype\n"

def get_value(exp, dat):
    res = re.findall(exp, dat)
    if len(res) == 1:
        ret = res.pop()
    elif len(res) == 0:
        ret = ""
    else:
        raise ValueError(f"Unexpected GTF format: {dat}")
    return ret


infile = open(snakemake.input[0])
outfile = open(snakemake.output[0], 'w')
outfile.write(HEADER)

for line_str in infile:
    if line_str.startswith("#"):
        continue
    line = line_str.split("\t")
    if line[2] == "gene":
        ret = [
            get_value(ID_RE, line[8]),
            get_value(NAME_RE, line[8]),
            get_value(BIOTYPE_RE, line[8])]
        outfile.write("\t".join(ret) + "\n")
    else:
        continue

infile.close()
outfile.close()
