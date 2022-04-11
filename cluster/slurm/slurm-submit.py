#!/usr/bin/env python3

import os
import re
import sys
import subprocess

from pathlib import Path
from snakemake.utils import read_job_properties

# Default log path
DEFAULT_LOG_PATH = "cluster/logs"

# SLURM Options, except for:
#  cpus-per-task
#  mem
opts = [
    "account", "array", "begin", "constraint", "dependency",
    "error", "job-name", "mail-type", "mail-user", "nodes",
    "ntasks", "output", "partition", "qos", "quiet", "time", "workdir"]

def get_val(val, default, dat, cluster):
    return dat.get(val, cluster.get(val, default))

def get_job_type(props):
    # NOTE: if group is present 'rule' is not valid
    if props["type"] == "group":
        ret = props["groupid"]

    else:
        ret = props["rule"]
    return ret

def get_job_token(props):
    wc = props.get("wildcards", dict())
    if wc:
        ret = "_".join(["{}={}".format(k, v) for k, v in wc.items()])
    else:
        ret = props.get("jobid", None)
    return ret

def get_job_name(props):
    return {"job-name": get_job_type(props)}

def get_resources(props):
    ret = dict()
    resources = props.get("resources", dict())

    # Grab mem
    for opt in ['mem', 'mem_mb', 'mem_gb']:
        val = get_val(opt, None, resources, cluster)
        if val:
            if opt == 'mem_gb':
                ret.update({"mem": str(val) + "G"})
            else:
                ret.update({"mem": str(val)})

    # Look for all over valid options
    for opt in opts:
        val = resources.get(opt, cluster.get(opt, None))
        if val:
            ret.update({opt: val})
    return ret

def get_threads(props):
    val = props.get("threads", cluster.get("cpus-per-task", None))
    if val:
        ret = {"cpus-per-task": val}
    else:
        ret = dict()
    return ret

def get_log_file(props):
    jobtype = str(get_job_type(props))
    jobtoken = str(get_job_token(props))
    if jobtoken:
        ret = jobtype + "-" + jobtoken + ".log"
    else:
        ret = jobtype + ".log"
    return ret

def get_logs(props):
    default = os.path.join(DEFAULT_LOG_PATH, get_log_file(props))
    out = cluster.get("output", "")
    err = cluster.get("error", "")
    ret = dict()

    if err != "":
        os.makedirs(os.path.dirname(err), exist_ok=True)
        ret.update({"error": err})
    if out == "":
        out = default

    os.makedirs(os.path.dirname(out), exist_ok=True)
    ret.update({'output': out})
    return ret

# Get job_properties and cluster info
jobscript = sys.argv[-1]
props = read_job_properties(jobscript)
cluster = props.get("cluster", dict())

# Process job_properties
args = get_resources(props)
args.update(get_threads(props))
args.update(get_job_name(props))
args.update(get_logs(props))

# Create command
opts = " ".join("--{} {}".format(k, v) for k, v in args.items())
cmd = "stdbuf -oL -eL sbatch --parsable {} {}".format(opts, jobscript)

# Run and print results
try:
    res = subprocess.run(cmd, check=True, shell=True, stdout=subprocess.PIPE)
except subprocess.CalledProcessError as e:
    raise e

try:
    print(res.stdout.decode())
except Exception as e:
    raise e

