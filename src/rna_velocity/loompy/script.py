#!/usr/bin/env python
import argparse
import datetime
import os
import pathlib
import re
import subprocess
import tarfile
import time
import csv

from collections import defaultdict

### VIASH START
par = {
    "taxon": "homo",
    "metadata": "path/to/metadata",
    "input": "path/to/fastqs",
    "output": "path/to/output",
    "reference_genome": "path/to/ref",
    "num_partitions": 10,
    "partition_id": None,
    "cell_count": 3000
}
### VIASH END

# valid reference genomes
reference_genomes_indexes = {
    "homo": "human_GRCh38_gencode.v31",
}

# other helpful constants
LOOMPY = "loompy"
CURR_MIN_VER = datetime.datetime(2017, 3, 1, tzinfo=datetime.timezone.utc)

run_dir = par["output"] / "data"
run_dir.mkdir(parents=True)

# extract sample name(s) and technology from the metadata tsv file
metadata_name = os.path.basename(par["metadata"])

technology, sample_name = "", ""

with open(par["metadata"]) as fd:
    rd = csv.reader(fd, delimiter="\t", quotechar='"')
    file_content = list()
    for row in rd:
        file_content.append(row)
    file_content = file_content[1:]
    sample_name = file_content[0][0]
    technology = file_content[0][
        1
    ]  # need to fix this later to fit tsv file with multiple samples

# check if the input genome is valid
if par["taxon"] in reference_genomes_indexes:
    genome_name = reference_genomes_indexes[par["taxon"]]
else:
    raise ValueError(f"unknown taxon {par["taxon"]}")


# extract valid fastq files
sample_re_smartseq2 = re.compile("([^/]+)_R\d(?:_\d+)?.fastq.gz$")
sample_re_10x = re.compile("([^/]+)_L\d+_R\d(?:_\d+)?.fastq.gz$")
s3_output_bucket, s3_output_prefix = s3u.s3_bucket_and_key(args.s3_output_path)

logger.info(
    "Running partition {} of {}".format(args.partition_id, args.num_partitions)
)

# fastq files are either stored directly under the s3 input folder, or in sample sub-folders under the s3 input folder
if list(s3u.get_files(s3_input_bucket, s3_input_prefix)):
    fastq_key_and_size = [
        (fn, s)
        for fn, s in s3u.get_size(s3_input_bucket, s3_input_prefix)
        if fn.endswith("fastq.gz")
    ]
else:
    sample_folder_paths = s3u.get_folders(s3_input_bucket, s3_input_prefix + "/")
    fastq_key_and_size = []
    for sample in sample_folder_paths:
        files = [
            (fn, s)
            for fn, s in s3u.get_size(s3_input_bucket, s3_input_prefix)
            if fn.endswith("fastq.gz")
        ]
        fastq_key_and_size += files

sample_name_to_fastq_keys = defaultdict(list)
fastq_sizes = defaultdict(list)
fastqs_key_to_name = dict()

for fn, s in fastq_key_and_size:
    matched = False
    if "10x" in technology:
        matched = sample_re_10x.search(os.path.basename(fn))
    elif "smartseq2" in technology:
        matched = sample_re_smartseq2.search(os.path.basename(fn))
    if matched:
        sample_name_to_fastq_keys[matched.group(1)].append(fn)
        fastq_sizes[matched.group(1)].append(s)
        fastqs_key_to_name[fn] = os.path.basename(fn)

logger.info(f"number of samples: {len(sample_name_to_fastq_keys)}")

# download input fastqs from S3 to an EC2 instance
fastq_dir = run_dir / "fastqs"
fastq_dir.mkdir(parents=True)

for key in fastqs_key_to_name.keys():
    s3c.download_file(
        Bucket=s3_input_bucket,
        Key=key,
        Filename=str(fastq_dir / fastqs_key_to_name[key]),
    )

# run kallisto alignment and RNA velocity analysis on the valid fastq files
for sample in sorted(sample_name_to_fastq_keys)[
    args.partition_id :: args.num_partitions
]:
    result_path = run_dir / "results"
    result_path.mkdir(parents=True)

    command = [
        "loompy",
        "fromfq",
        str(result_path / f"{sample_name}.loom"),
        sample_name,
        str(par["reference_genome"]),
        str(par["metadata"]),
    ]
    fastq_names = [
        fastqs_key_to_name[key] for key in sample_name_to_fastq_keys[sample]
    ]
    fastq_dirs = [str(fastq_dir / fastq) for fastq in fastq_names]
    command += fastq_dirs
    print(command)  # for testing purpose

    failed = ut_log.log_command(
        logger,
        command,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        universal_newlines=True,
    )

    t_config = TransferConfig(use_threads=False)
    if failed:
        raise RuntimeError("loompy failed")
    else:
        logger.info(f"Uploading {sample_name}.loom")
        s3c.upload_file(
            Filename=str(result_path / f"{sample_name}.loom"),
            Bucket=s3_output_bucket,
            Key=os.path.join(s3_output_prefix, f"{sample_name}.loom"),
            Config=t_config,
        )

    command = ["rm", "-rf", str(result_path)]
    ut_log.log_command(logger, command, shell=True)

    time.sleep(30)

logger.info("Job completed")
