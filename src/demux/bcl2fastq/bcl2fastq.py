#!/usr/bin/env python

# Dependencies:
# pip: argparse, glob, os, re, subprocess, sys
# conda: bcl2fastq
# apt: awscli

import sys

import argparse
import glob
import os
import re
import subprocess
import sys
from unittest import result

## VIASH START

par = {

}

meta = { 'resources_dir': '../../utilities' }

## VIASH END

sys.path.append(meta['resources_dir'])
from utilities.log_util import get_logger, log_command

BCL2FASTQ = "bcl2fastq"

def main(logger):

    S3_RETRY = 5

    root_dir = "./"

    print(par)

    if par["sample_sheet_name"] is None:
        par["sample_sheet_name"] = "{}.csv".format(par["exp_id"])

    # local directories
    result_path = os.path.join(root_dir, "data", "hca", par["exp_id"])
    bcl_path = os.path.join(result_path, "bcl")
    output_path = os.path.join(result_path, "fastqs")

    # download sample sheet
    os.makedirs(result_path)
    os.mkdir(bcl_path)

    # if sample sheet is on s3
    if par["input_path"][0:2] == "s3":
        command = [
            "aws",
            "s3",
            "cp",
            "--quiet",
            os.path.join(par["sample_sheet"], par["sample_sheet_name"]),
            result_path,
        ]
        for i in range(S3_RETRY):
            if not log_command(logger, command, shell=True):
                break
            logger.info("retrying s3 copy")
        else:
            raise RuntimeError(
                "couldn't download sample sheet {}".format(
                    os.path.join(par["sample_sheet"], par["sample_sheet_name"])
                )
            )

        # download the bcl files
        command = [
            "aws",
            "s3",
            "sync",
            "--quiet",
            "--force-glacier-transfer" if par["force_glacier"] else "",
            os.path.join(par["input_path"], par["exp_id"]),
            bcl_path,
        ]
        for i in range(S3_RETRY):
            if not log_command(logger, command, shell=True):
                break
            logger.info("retrying s3 sync bcl")
        else:
            raise RuntimeError(
                "couldn't sync {}".format(os.path.join(par["input_path"], par["exp_id"]))
            )
    else:
        # copying data locally into correct directories
        command = [
            "cp",
            os.path.join(par["sample_sheet"], par["sample_sheet_name"]),
            result_path,
        ]
        log_command(logger, command, shell=True)

        command = [
            "cp",
            os.path.join(par["input_path"], par["exp_id"]),
            bcl_path,
        ]
        log_command(logger, command, shell=True)
        

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
        BCL2FASTQ,
        " ".join(par["bcl2fastq_options"]),
        "--sample-sheet",
        os.path.join(result_path, par["sample_sheet_name"]),
        "-R",
        bcl_path,
        "-o",
        output_path,
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

    if par["s3"]:
        # upload fastq files to destination folder
        command = [
            "aws",
            "s3",
            "sync",
            "--quiet",
            output_path,
            os.path.join(par["output_path"], par["exp_id"]),
            "--exclude",
            '"*"',
            "--include",
            '"*fastq.gz"',
        ]
        for i in range(S3_RETRY):
            if not log_command(logger, command, shell=True):
                break
            logger.info("retrying sync fastq")
        else:
            raise RuntimeError("couldn't sync fastqs")

        # Move reports data back to S3
        reports_path = subprocess.check_output(
            "ls -d {}".format(
                os.path.join(output_path, "Reports", "html", "*", "all", "all", "all")
            ),
            shell=True,
        ).rstrip()
        command = [
            "aws",
            "s3",
            "cp",
            "--quiet",
            reports_path.decode(),
            os.path.join(par["report_path"], par["exp_id"]),
            "--recursive",
        ]
        for i in range(S3_RETRY):
            if not log_command(logger, command, shell=True):
                break
            logger.info("retrying cp reports")
        else:
            raise RuntimeError("couldn't cp reports")

    p.kill()


if __name__ == "__main__":
    mainlogger, log_file, file_handler = get_logger(__name__, debug=True)
    main(mainlogger)