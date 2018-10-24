#!/usr/bin/env python3

import argparse


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--branch", default="master", help="branch of utilities repo to use"
    )

    parser.add_argument("taxon", choices=("homo",))
    parser.add_argument("num_partitions", type=int)
    parser.add_argument(
        "input_dirs", nargs="+", help="The folder(s) containing bam files to velocytize"
    )

    parser.add_argument(
        "script_args",
        nargs=argparse.REMAINDER,
        help="Extra arguments are passed to velocyto",
    )

    args = parser.parse_args()

    for i in range(args.num_partitions):
        print(
            " ".join(
                (
                    "evros",
                    f"--branch {args.branch}",
                    "alignment.velocyto",
                    f"--taxon {args.taxon}",
                    f"--num_partitions {args.num_partitions}",
                    f"--partition_id {i}",
                    f"--input_dirs {' '.join(args.input_dirs)}",
                    " ".join(args.script_args),
                )
            )
        )
        print("sleep 10")
