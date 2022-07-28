#!/usr/bin/env python3

import argparse
import warnings

from utilities.alignment.run_star_and_htseq import reference_genomes, deprecated


def main():
    parser = argparse.ArgumentParser(
        description="Create a shell to run alignment jobs"
        " with STAR and htseq for multiple samples all together"
    )

    # required arguments
    requiredNamed = parser.add_argument_group("required arguments")

    requiredNamed.add_argument(
        "--taxon",
        choices=list(reference_genomes.keys()),
        required=True,
        help="Reference genome for the alignment run, "
        "selected from the reference_genomes dictionary keys from "
        "alignment.run_star_and_htseq.py",
    )
    requiredNamed.add_argument(
        "--num_partitions",
        type=int,
        required=True,
        help="Number of groups to divide samples into for the alignment run",
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
        help="Extra arguments are passed to run_star_and_htseq",
    )

    args = parser.parse_args()

    # check if the input genome is valid
    if args.taxon not in reference_genomes:
        raise ValueError(f"{args.taxon} is currently unavailable for velocyto run")
        
    # if args.taxon in reference_genomes:
    #     if args.taxon in deprecated:
    #         warnings.warn(
    #             f"The name '{args.taxon}' will be removed in the future,"
    #             f" start using '{deprecated[args.taxon]}'"
    #         )
    # else:
    #     raise ValueError(f"unknown taxon {args.taxon}")

    # print input arguments for running alignment.run_star_and_htseq for each group of sample
    for i in range(args.num_partitions):
        print(
            " ".join(
                (
                    "evros",
                    f"--branch {args.branch}",
                    "alignment.run_star_and_htseq",
                    f"--taxon {args.taxon}",
                    f"--num_partitions {args.num_partitions}",
                    f"--partition_id {i}",
                    f"--s3_input_path {args.s3_input_path}",
                    f"--s3_output_path {args.s3_output_path}",
                    " ".join(args.script_args),
                )
            )
        )
        print("sleep 10")
