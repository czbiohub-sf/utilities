#!/usr/bin/env python3

import argparse

from utilities.velocyto.run_velocyto_star import reference_genomes
import utilities.s3_util as s3u


def main():
    parser = argparse.ArgumentParser(
        description="Create a shell script locally to run velocyto on smartseq2 data aligned with STAR"
    )

    # required arguments
    requiredNamed = parser.add_argument_group("required arguments")

    requiredNamed.add_argument(
        "--taxon",
        choices=("hg38-plus", "mm10-plus"),
        required=True,
        help="Reference genome for the velocyto run on smartseq2 data aligned with STAR",
    )

    requiredNamed.add_argument(
        "--s3_input_path",
        required=True,
        help="Location of multiple input sample folders with STAR alignment results",
    )

    requiredNamed.add_argument(
        "--s3_output_path", required=True, help="Location for output",
    )

    requiredNamed.add_argument("--num_partitions", type=int)
    requiredNamed.add_argument(
        "--input_dirs",
        nargs="+",
        help="The folder(s) containing bam files to velocytize",
    )

    # optional arguments
    parser.add_argument(
        "--branch", default="master", help="branch of utilities repo to use"
    )

    parser.add_argument(
        "script_args",
        nargs=argparse.REMAINDER,
        help="Extra arguments are passed to velocyto",
    )

    args = parser.parse_args()

    # check if the input genome is valid
    if args.taxon not in reference_genomes:
        raise ValueError(f"{args.taxon} is currently unavailable for velocyto run")

    for i in range(args.num_partitions):
        print(
            " ".join(
                (
                    "evros",
                    f"--branch {args.branch}",
                    "rna_velocity.run_velocyto_star",
                    f"--taxon {args.taxon}",
                    f"--num_partitions {args.num_partitions}",
                    f"--partition_id {i}",
                    f"--s3_input_path {args.s3_input_path}",
                    f"--s3_output_path {args.s3_output_path}",
                    f"--input_dirs {' '.join(args.input_dirs)}",
                    " ".join(args.script_args),
                )
            )
        )
        print("sleep 10")
