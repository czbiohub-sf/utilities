
import sys
import glob
import os
import re
import subprocess

## VIASH START
par = {
    "sample_sheet" : None,
    "input": "../resources_test/bs_195891710/bcl_data/",
    "output": "../resources_test/bs_195891710/fastqs/",
    "reports": "../resources_test/bs_195891710/reports/",
}
## VIASH END

if not os.path.isdir(par["output"]):
    print(f"creating output directory {par['output']}")
    os.makedirs(par["output"])

# construct command args
print("Help")
command = [
    "cellranger mkfastq",
    "--run=" + par["input"],
    "--output-dir=" + par["output"],
]
if par["sample_sheet"] is not None:
    command = command + [ "--sample-sheet="+par["sample_sheet"] ]

# run 10x_mkfastq
with subprocess.Popen(
    command,
    stdout=subprocess.PIPE,
    stderr=subprocess.STDOUT,
    shell=True,
) as p:
    for line in p.stdout:
        print(line.decode(), end="")

if p.returncode > 0:
    raise RuntimeError(f"10x_mkfastq failed with exit code {p.returncode}")

print("demux/10x_mkfastq run completed.")
sys.stdout.flush()