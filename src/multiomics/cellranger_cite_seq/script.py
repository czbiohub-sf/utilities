#!/usr/bin/env python

import os
import subprocess

### VIASH START
par = {
    "feature_ref_path": "path/to/feature",
    "expect_cells": 3000,
    "run_id": "1",
    "libraries_path": "path/to/libraries",
    "reference_genome": "path/to/reference_genome",
    "output": "path/to/results",
}
### VIASH END


CELLRANGER = "cellranger"

run_id = par["run_id"]

# Run cellranger
os.chdir(str(par["output"]))
command = [
    CELLRANGER,
    "count",
    f"--id={par['run_id']}",
    f"--transcriptome={par['reference_genome']}",
    f"--libraries={par['libraries_path']}",
    f"--feature-ref={par['feature_ref_path']}",
    f"--expect-cells={par['expect_cells']}",
    "--localmem=256",
    "--localcores=64",
]

# run cellranger-arc count
with subprocess.Popen(
    command,
    stdout=subprocess.PIPE,
    stderr=subprocess.STDOUT,
) as p:
    for line in p.stdout:
        print(line.decode(), end="")

if p.returncode > 0:
    raise RuntimeError(f"cellranger count failed with exit code {p.returncode}")