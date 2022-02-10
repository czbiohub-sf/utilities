#!/usr/bin/env python

# Dependencies:
# pip: argparse, glob, os, re, subprocess, sys
# conda: bcl2fastq
# apt: awscli

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
    "output_reports": ""
}

meta = { 'resources_dir': '../../helpers' }

## VIASH END

sys.path.append(meta['resources_dir'])
from helpers.log_util import log_command

# Run bcl2 fastq
command = [
    "bcl2fastq",
    "--sample-sheet", par["sample_sheet"],
    "--runfolder-dir", par["input"],
    "--output-dir", par["output"],
]
failed = log_command(
    logger, command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True
)
if failed:
    p.kill()
    raise RuntimeError("bcl2fastq failed, see above for error")

# fix directory structure of the files *before* sync!
fastqgz_files = glob.glob(os.path.join(output_path, "*fastq.gz"))
logger.debug("all fastq.gz files\n{}\n\n".format("\n".join(fastqgz_files)))

for fastq_file in fastqgz_files:
    if par["skip_undetermined"] and os.path.basename(fastq_file).startswith(
        "Undetermined"
    ):
        logger.info("removing {}".format(os.path.basename(fastq_file)))
        os.remove(fastq_file)
    elif par["star_structure"]:
        m = re.match("(.+)(_R[12]_001.fastq.gz)", os.path.basename(fastq_file))
        if m:
            sample = m.group(1)
            if not os.path.exists(os.path.join(output_path, sample)):
                logger.debug(
                    "creating {}".format(os.path.join(output_path, sample))
                )
                os.mkdir(os.path.join(output_path, sample))
            logger.debug("moving {}".format(fastq_file))
            os.rename(
                fastq_file,
                os.path.join(output_path, sample, os.path.basename(fastq_file)),
            )
        else:
            logger.warning("Warning: regex didn't match {}".format(fastq_file))

sys.stdout.flush()

# Move reports data back to S3
# reports_path = subprocess.check_output(
#     "ls -d {}".format(
#         os.path.join(output_path, "Reports", "html", "*", "all", "all", "all")
#     ),
#     shell=True,
# ).rstrip()
# command = [
#     "aws",
#     "s3",
#     "cp",
#     "--quiet",
#     reports_path.decode(),
#     os.path.join(par["report_path"], par["exp_id"]),
#     "--recursive",
# ]
# for i in range(S3_RETRY):
#     if not log_command(logger, command, shell=True):
#         break
#     logger.info("retrying cp reports")
# else:
#     raise RuntimeError("couldn't cp reports")
