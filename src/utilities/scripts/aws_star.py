#!/usr/bin/env python3

import argparse

reference_genomes = {
    "homo": "HG38-PLUS",
    "hg38-plus": "HG38-PLUS",
    "mus": "MM10-PLUS",
    "mm10-plus": "MM10-PLUS",
    "microcebus": "MicMur3-PLUS",
    "gencode.vM19": "gencode.vM19",
}

deprecated = {"homo", "mus", "mus-premrna"}


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--branch", default="master", help="branch of utilities repo to use"
    )

    parser.add_argument("taxon", choices=list(reference_genomes.keys()))
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

    if args.taxon in reference_genomes:
        if args.taxon in deprecated:
            print(
                f"The name '{args.taxon}' will be removed in the future,"
                f" start using '{reference_genomes[args.taxon]}'"
            )
    else:
        raise ValueError(f"unknown taxon {args.taxon}")

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
                    f"--input_dirs {' '.join(args.input_dirs)}",
                    " ".join(args.script_args),
                )
            )
        )
        print("sleep 10")
