#!/usr/bin/env python3

import argparse
import warnings
import posixpath

from utilities.alignment.run_10x_count import reference_genomes, deprecated
import utilities.s3_util as s3u


def main():
    parser = argparse.ArgumentParser(
        description="Create a shell to run alignment jobs"
        " with 10x for multiple samples all together"
    )

    # required arguments
    requiredNamed = parser.add_argument_group("required arguments")

    requiredNamed.add_argument(
        "--taxon",
        choices=list(reference_genomes.keys()),
        required=True,
        help="Reference genome for the alignment run, "
        "selected from the reference_genomes dictionary keys from "
        "alignment.run_10x_count.py",
    )

    requiredNamed.add_argument(
        "--s3_input_path",
        required=True,
        help="The folder containing sample folders, "
        "each of which have fastq.gz files to align",
    )

    requiredNamed.add_argument(
        "--s3_output_path",
        required=True,
        help="The folder to store the alignment results",
    )

    # optional arguments
    parser.add_argument(
        "--branch", default="master", help="Branch of utilities repo to use"
    )

    parser.add_argument(
        "script_args",
        nargs=argparse.REMAINDER,
        help="Extra arguments are passed to run_10x_count",
    )

    args = parser.parse_args()

    # check if the input genome is valid
    if args.taxon in reference_genomes:
        if args.taxon in deprecated:
            warnings.warn(
                f"The name '{args.taxon}' will be removed in the future,"
                f" start using '{deprecated[args.taxon]}'"
            )
    else:
        raise ValueError(f"unknown taxon {args.taxon}")

    # get the list of sample folder paths under the input folder
    s3_input_bucket, s3_input_prefix = s3u.s3_bucket_and_key(args.s3_input_path)
    s3_input_prefix += "/"
    sample_folder_paths = [
        folder_path for folder_path in s3u.get_folders(s3_input_bucket, s3_input_prefix)
    ]
    complete_input_paths = [
        "s3://" + s3_input_bucket + "/" + path for path in sample_folder_paths
    ]

    # print input arguments of running alignment.run_10x_count for each sample folder
    num_partitions = len(complete_input_paths)
    for i in range(num_partitions):
        s3_input_path = complete_input_paths[i]
        print(
            " ".join(
                (
                    "evros",
                    f"--branch {args.branch}",
                    "alignment.run_10x_count",
                    f"--taxon {args.taxon}",
                    f"--num_partitions {num_partitions}",
                    f"--partition_id {i}",
                    f"--s3_input_path {s3_input_path}",
                    f"--s3_output_path {args.s3_output_path}",
                    " ".join(args.script_args),
                )
            )
        )
        print("sleep 10")
