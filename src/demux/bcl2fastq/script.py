#!/usr/bin/env python

# Dependencies:
# pip: argparse, glob, os, re, subprocess, sys

import sys

import glob
import os
import re
import subprocess
import sys

## VIASH START

par = {
    "sample_sheet" : "resources_test/bcl/sample_sheet.csv",
    "input": "resources_test/bcl/input/",
    "output": "resources_test/bcl/output/",
    "reports": "resources_test/bcl/reports/",
    "skip_undetermined": True,
    "star_structure": True,
}

## VIASH END

# Run memory monitor
command = (
    "while true;"
    ' do echo "memory usage" `cat /sys/fs/cgroup/memory/memory.usage_in_bytes`;'
    ' echo "disk usage" `df -h | grep "/mnt"`;'
    " sleep 300;"
    " done"
)
p = subprocess.Popen([command], shell=True)

# Run bcl2 fastq
command = [
    "bcl2fastq",
    "--sample-sheet", par["sample_sheet"],
    "--runfolder-dir", par["input"],
    "--output-dir", par["output"],
]
pipe = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
if pipe.communicate()[1]:
    p.kill()
    raise RuntimeError("bcl2fastq failed: " + str(pipe[1]))

# Fix directory structure of the files
fastqgz_files = glob.glob(os.path.join(par['output'], "*fastq.gz"))
print("all fastq.gz files\n{}\n\n".format("\n".join(fastqgz_files)))

for fastq_file in fastqgz_files:
    if par["skip_undetermined"] and os.path.basename(fastq_file).startswith(
        "Undetermined"
    ):
        print("removing {}".format(os.path.basename(fastq_file)))
        os.remove(fastq_file)
    elif par["star_structure"]:
        m = re.match("(.+)(_R[12]_001.fastq.gz)", os.path.basename(fastq_file))
        if m:
            sample = m.group(1)
            if not os.path.exists(os.path.join(par['output'], sample)):
                print(
                    "creating {}".format(os.path.join(par['output'], sample))
                )
                os.mkdir(os.path.join(par['output'], sample))
            print("moving {}".format(fastq_file))
            os.rename(
                fastq_file,
                os.path.join(par['output'], sample, os.path.basename(fastq_file)),
            )
        else:
            print("Warning: regex didn't match {}".format(fastq_file))

sys.stdout.flush()
print("demux/bcl2fastq run completed.")