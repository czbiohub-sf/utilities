#!/usr/bin/env python3

import argparse

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("taxon", choices=("mus", "homo"))
    parser.add_argument("num_partitions", type=int)
    parser.add_argument(
        "input_dirs", nargs="+", help="The folder(s) containing fastq.gz files to align"
    )
    parser.add_argument(
        "script_args",
        nargs=argparse.REMAINDER,
        help="Extra arguments are passed to run_star_and_htseq",
    )

    args = parser.parse_args()

    for i in range(args.num_partitions):
        print(
            " ".join(
                (
                    "evros",
                    "alignment.run_star_and_htseq",
                    "--taxon {}".format(args.taxon),
                    "--num_partitions {}".format(args.num_partitions),
                    "--partition_id {}".format(i),
                    "--input_dirs {}".format(" ".join(args.input_dirs)),
                    " ".join(args.script_args),
                )
            )
        )
        print("sleep 10")
