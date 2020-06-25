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

import utilities.log_util as ut_log
import utilities.s3_util as s3u

import boto3
from boto3.s3.transfer import TransferConfig


# reference genome bucket name for different regions
S3_REFERENCE = {"east": "czbiohub-reference-east", "west": "czbiohub-reference"}

# valid reference genomes
reference_genomes_indexes = {
    "homo": "human_GRCh38_gencode.v31",
}

# other helpful constants
LOOMPY = "loompy"
CURR_MIN_VER = datetime.datetime(2017, 3, 1, tzinfo=datetime.timezone.utc)


def get_default_requirements():
    return argparse.Namespace(vcpus=16, memory=64000, storage=500, ecr_image="velocyto")


def get_parser():
    """ Construct and return the ArgumentParser object that parses input command.
    """

    parser = argparse.ArgumentParser(
        prog="run_loompy.py",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Run kallisto alignment and RNA velocity analysis with loompy",
    )

    # required arguments
    requiredNamed = parser.add_argument_group("required arguments")

    requiredNamed.add_argument(
        "--taxon",
        required=True,
        choices=list(reference_genomes_indexes.keys()),
        help="Reference genome index for the alignment and RNA velocity run",
    )

    requiredNamed.add_argument(
        "--metadata",
        required=True,
        help="Path to the metadata file with sample name, technology, and target cell numbers",
    )

    requiredNamed.add_argument(
        "--s3_input_path",
        required=True,
        help="The folder with fastq.gz files to align and perform RNA velocity on. Can either have plain fastq files from multiple samples or include sample subfolders.",
    )

    requiredNamed.add_argument(
        "--s3_output_path",
        required=True,
        help="The folder to store the resulting loom file",
    )

    requiredNamed.add_argument(
        "--num_partitions",
        type=int,
        required=True,
        default=10,
        help="Number of groups to divide samples " "into for the loompy run",
    )

    requiredNamed.add_argument(
        "--partition_id", type=int, required=True, help="Index of the sample group",
    )

    # optional arguments
    parser.add_argument("--cell_count", type=int, default=3000)
    parser.add_argument("--glacier", action="store_true")
    parser.add_argument("--root_dir", default="/mnt")

    return parser


def main(logger):
    """ Download reference genome from the S3 bucket to an EC2 instance, run alignment jobs with kallisto, then generate loom files based on the alignment results, and upload the loom files to S3 bucket.

        logger - Logger object that exposes the interface the code directly uses
    """

    parser = get_parser()

    args = parser.parse_args()

    root_dir = pathlib.Path(args.root_dir)

    if os.environ.get("AWS_BATCH_JOB_ID"):
        root_dir = root_dir / os.environ["AWS_BATCH_JOB_ID"]

    # local directories on the EC2 instance
    if args.s3_input_path.endswith("/"):
        args.s3_input_path = args.s3_input_path[:-1]

    run_dir = root_dir / "data"
    run_dir.mkdir(parents=True)

    # extract sample name(s) and technology from the metadata tsv file
    metadata_name = os.path.basename(args.metadata)
    metadata_dir = run_dir / "metadata"
    metadata_dir.mkdir(parents=True)
    metadata_dir = metadata_dir / metadata_name
    s3_metadata_bucket, s3_metadata_prefix = s3u.s3_bucket_and_key(args.metadata)

    s3c.download_file(
        Bucket=s3_metadata_bucket,  # just always download this from us-west-2...
        Key=s3_metadata_prefix,
        Filename=str(metadata_dir),
    )

    technology, sample_name = "", ""
    with open(metadata_dir) as fd:
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
    if args.taxon in reference_genomes_indexes:
        genome_name = reference_genomes_indexes[args.taxon]
    else:
        raise ValueError(f"unknown taxon {args.taxon}")

    if "10x" in technology:
        genome_dir = root_dir / "genome" / "10X" / genome_name
    elif (
        "smartseq2" in technology
    ):  # may need to update these after confirming what technology name looks like for smartseq2 data
        genome_dir = (
            root_dir / "genome" / "smartseq2" / genome_name
        )  # necessary to separate the reference genome location path for 10x and smartseq2?
    genome_dir.mkdir(parents=True)

    s3_input_bucket, s3_input_prefix = s3u.s3_bucket_and_key(args.s3_input_path)

    logger.info(
        f"""Run Info: partition {args.partition_id} out of {args.num_partitions}
                   genome_dir:\t{genome_dir}
                        taxon:\t{args.taxon}
                s3_input_path:\t{args.s3_input_path}"""
    )

    s3 = boto3.resource("s3")

    # download the reference genome index data
    logger.info("Downloading reference genome index files of {}".format(genome_name))

    if "10x" in technology:
        s3_genome_index = f"s3://{S3_REFERENCE['west']}/loompy/10X/{genome_name}"
    elif "smartseq2" in technology:
        s3_genome_index = f"s3://{S3_REFERENCE['west']}/loompy/smartseq2/{genome_name}"

    s3_genome_index_bucket, s3_genome_index_prefix = s3u.s3_bucket_and_key(
        s3_genome_index
    )
    s3_genome_files_prefix = list(
        s3u.get_files(s3_genome_index_bucket, s3_genome_index_prefix)
    )
    s3_genome_files_prefix = s3_genome_files_prefix[1:]
    file_names = list(
        os.path.basename(file_path) for file_path in s3_genome_files_prefix
    )
    genome_name_to_prefix = dict(zip(file_names, s3_genome_files_prefix))
    
    for file in genome_name_to_prefix.keys():
        s3c.download_file(
            Bucket=s3_genome_index_bucket,
            Key=genome_name_to_prefix[file],
            Filename=str(genome_dir / file),
        )

    # extract valid fastq files
    sample_re_smartseq2 = re.compile("([^/]+)_R\d(?:_\d+)?.fastq.gz$")
    sample_re_10x = re.compile("([^/]+)_L\d+_R\d(?:_\d+)?.fastq.gz$")
    s3_output_bucket, s3_output_prefix = s3u.s3_bucket_and_key(args.s3_output_path)

    logger.info(
        "Running partition {} of {}".format(args.partition_id, args.num_partitions)
    )

    # fastq files are stored all under the s3 input folder, or in sample sub-folders
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

    sample_name_to_fastq_key_and_size = defaultdict(list)
    fastq_sizes = defaultdict(list)
    fastq_names = list()

    for fn, s in fastq_key_and_size:
        fastq_names.append(os.path.basename(fn))
        matched = False
        if "10x" in technology:
            matched = sample_re_10x.search(os.path.basename(fn))
        elif "smartseq2" in technology:
            matched = sample_re_smartseq2.search(os.path.basename(fn))
        if matched:
            sample_name_to_fastq_key_and_size[matched.group(1)].append(fn)
            fastq_sizes[matched.group(1)].append(s)

    logger.info(f"number of samples: {len(sample_name_to_fastq_key_and_size)}")
    
    # download input fastqs from S3 to an EC2 instance
    fastq_dir = run_dir / "fastqs"
    fastq_dir.mkdir(parents=True)
    assert len(fastq_names) == len(fastq_key_and_size), 'fastq_names and fastq_key_and_size should have same length'
    fastqs_key_to_name = dict(zip(fastq_key_and_size, fastq_names))
    
    for key_size in fastq_key_and_size:
        s3c.download_file(
            Bucket=s3_input_bucket,
            Key=key_size[0],
            Filename=str(fastq_dir / fastqs_key_to_name[key_size]),
        )

    # run kallisto alignment and RNA velocity analysis on the valid fastq files
    for sample in sorted(sample_name_to_fastq_key_and_size)[args.partition_id :: args.num_partitions]:
        result_path = run_dir / "results"
        result_path.mkdir(parents=True)
        os.chdir(result_path)

        command = [
            "loompy",
            "fromfq",
            f"{sample_name}.loom",
            sample_name,
            str(genome_dir),
            str(metadata_dir),
        ]
        fastqs = [str(fastq_dir / fastq) for fastq in fastqs_key_to_name[sample_name_to_fastq_key_and_size[sample][0]]]
        command += fastqs

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

        command = ["rm", "-rf", result_path]
        ut_log.log_command(logger, command, shell=True)

        time.sleep(30)

    logger.info("Job completed")


if __name__ == "__main__":
    mainlogger, log_file, file_handler = ut_log.get_logger(__name__)
    s3c = boto3.client("s3")
    main(mainlogger)
