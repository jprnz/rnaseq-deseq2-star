#!/bin/bash

path=$(readlink -f $(dirname $0))
tmpdir="$path/.temp"

[ -e $tmpdir ] || mkdir -p $tmpdir
export TMPDIR=$tmpdir

snakemake \
  --latency-wait 240 \
  --conda-frontend mamba \
  --profile $path/cluster/slurm \
  "$@"
