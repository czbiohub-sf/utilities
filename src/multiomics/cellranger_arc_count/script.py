#!/usr/bin/env python

import os
import subprocess

### VIASH START
par = {
    "peaks": "path/to/bedfile",
    "libraries_path": "path/to/lib",
    "run_id": "1",
    "reference_genome": "path/to/reference_genome",
    "output": "path/to/results",
}
### VIASH END

CELLRANGER = "cellranger-arc" 

run_id = par["run_id"]

os.chdir(par["output"])
command = [
    CELLRANGER,
    "count",
    f"--id={run_id}",
    f"--reference={par['reference_genome']}",
    f"--libraries={par['libraries_path']}",
    "--localmem=256",
    "--localcores=64",
]

if par["peaks"]:
    peaks = par["peaks"]
    command.append(f"--peaks={str(peaks)}")

# run cellranger-arc count
with subprocess.Popen(
    command,
    stdout=subprocess.PIPE,
    stderr=subprocess.STDOUT,
) as p:
    for line in p.stdout:
        print(line.decode(), end="")

if p.returncode > 0:
    raise RuntimeError(f"cellranger-arc count failed with exit code {p.returncode}")
