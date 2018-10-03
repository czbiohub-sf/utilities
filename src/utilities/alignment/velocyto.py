#!/usr/bin/env python
import argparse
import os
import re
import subprocess

import utilities.log_util as ut_log
import utilities.s3_util as s3u

import boto3
from boto3.s3.transfer import TransferConfig


def get_default_requirements():
    return argparse.Namespace(vcpus=2, memory=8000, storage=500, ecr_image="velocyto")


def get_parser():
    parser = argparse.ArgumentParser(
        prog="velocyto.py", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument("--taxon", choices=("homo",), default="homo")

    parser.add_argument(
        "--s3_input_path",
        default="s3://czbiohub-seqbot/fastqs",
        help="Location of input folders",
    )
    parser.add_argument("--s3_output_path", required=True, help="Location for output")

    parser.add_argument("--num_partitions", type=int, required=True)
    parser.add_argument("--partition_id", type=int, required=True)
    parser.add_argument(
        "--input_dirs",
        nargs="+",
        required=True,
        help="List of input folders to process",
    )
    parser.add_argument(
        "--plates", nargs="+", required=True, help="List of plates to run"
    )

    return parser


def run_sample(
    sample_key,
    mask_path,
    gtf_path,
    s3_input_bucket,
    s3_output_bucket,
    s3_output_prefix,
    run_dir,
    logger,
):

    s3c = boto3.client("s3")
    t_config = TransferConfig(use_threads=False, num_download_attempts=25)

    sample_name = os.path.basename(sample_key)
    sample_id = sample_name.split(".")[0]  # this is brittle!
    local_sample = os.path.join(run_dir, "input", sample_name)

    s3c.download_file(
        Bucket=s3_input_bucket, Key=sample_key, Filename=local_sample, Config=t_config
    )

    veloctyo_command = [
        "velocyto",
        "run-smartseq2",
        "-o",
        run_dir,
        "-m",
        mask_path,
        "-e",
        sample_id,
        local_sample,
        gtf_path,
    ]

    # try:
    output = subprocess.check_output(" ".join(veloctyo_command), shell=True)
    # except subprocess.CalledProcessError:
    #     output = "Command failed!"

    output_file = os.path.join(run_dir, "velocyto", f"{sample_id}.loom")

    logger.info("Uploading {}".format(output_file))
    s3c.upload_file(
        Filename=output_file,
        Bucket=s3_output_bucket,
        Key=s3_output_prefix,
        Config=t_config,
    )

    os.remove(local_sample)
    os.remove(output_file)


def main(logger):
    parser = get_parser()

    args = parser.parse_args()

    if os.environ.get("AWS_BATCH_JOB_ID"):
        root_dir = os.path.join("/mnt", os.environ["AWS_BATCH_JOB_ID"])
    else:
        root_dir = "/mnt"

    run_dir = os.path.join(root_dir, "data")
    os.makedirs(run_dir)

    os.mkdir(os.path.join(run_dir, "reference"))
    os.mkdir(os.path.join(run_dir, "input"))

    if args.taxon == "homo":
        gtf_file = "hg38-plus.gtf"
        mask_file = "hg38_rmsk.gtf"
    elif args.taxon == "mus":
        gtf_file = os.path.join(root_dir, "mm10-plus.gtf")
        mask_file = os.path.join(root_dir, "mm10_rmsk.gtf")
    else:
        raise ValueError("Invalid taxon {}".format(args.taxon))

    s3_input_bucket, s3_input_prefix = s3u.s3_bucket_and_key(args.s3_input_path)

    logger.info(
        """Run Info: partition {} out of {}
                    gtf_file:\t{}
                   mask_file:\t{}
               s3_input_path:\t{}
                  input_dirs:\t{}""".format(
            args.partition_id,
            args.num_partitions,
            gtf_file,
            mask_file,
            args.taxon,
            args.s3_input_path,
            ", ".join(args.input_dirs),
        )
    )

    gtf_path = os.path.join(run_dir, "reference", gtf_file)
    mask_path = os.path.join(run_dir, "reference", mask_file)
    s3u.download_files(
        [f"velocyto/{gtf_file}", f"velocyto/{mask_file}"],
        [gtf_path, mask_path],
        b="czbiohub-reference",
        n_proc=2,
    )

    sample_re = re.compile(f"([^/]+).{args.taxon}.Aligned.out.sorted.bam$")
    plate_set = set(args.plates)

    s3_output_bucket, s3_output_prefix = s3u.s3_bucket_and_key(args.s3_output_path)

    for input_dir in args.input_dirs:
        logger.info(
            "Running partition {} of {} for {}".format(
                args.partition_id, args.num_partitions, input_dir
            )
        )

        sample_files = [
            fn
            for fn in s3u.get_files(
                s3_input_bucket, os.path.join(s3_input_prefix, input_dir)
            )
            if fn.endswith(f"{args.taxon}.Aligned.out.sorted.bam")
        ]

        plate_samples = []

        for fn in sample_files:
            matched = sample_re.search(os.path.basename(fn))
            if matched.group(1).split("_")[1] in plate_set:
                plate_samples.append(fn)

        logger.info(f"number of bam files: {len(plate_samples)}")

        for sample_name in sorted(plate_samples)[
            args.partition_id :: args.num_partitions
        ]:
            run_sample(
                sample_name,
                mask_path,
                gtf_path,
                s3_input_bucket,
                s3_output_bucket,
                os.path.join(s3_output_prefix, input_dir),
                run_dir,
                logger,
            )

    logger.info("Job completed")


if __name__ == "__main__":
    mainlogger, log_file, file_handler = ut_log.get_logger(__name__)

    main(mainlogger)
