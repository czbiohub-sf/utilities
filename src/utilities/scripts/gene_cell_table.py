#!/usr/bin/env python

import argparse
import csv
import io
import os

import boto3

from utilities.log_util import get_logger


def get_htseq_counts(client, bucket, htseq_file):
    fb = io.BytesIO()

    client.download_fileobj(Bucket=bucket, Key=htseq_file, Fileobj=fb)

    fb.seek(0)
    gene_list, counts = list(
        zip(*[map(str.strip, line.decode().split("\t")) for line in fb])
    )

    return gene_list, counts


def get_log_file(client, bucket, log_file):
    fb = io.BytesIO()

    client.download_fileobj(Bucket=bucket, Key=log_file, Fileobj=fb)

    fb.seek(0)
    metric_names, values = list(
        zip(*[map(str.strip, line.decode().split("|")) for line in fb if b"|" in line])
    )

    return metric_names, values


def gene_cell_table(args, logger, dryrun):
    logger.info("Starting")

    if args.output_file.endswith(".txt"):
        sep = "\t"
    elif args.output_file.endswith(".csv"):
        sep = ","
    else:
        raise ValueError(
            "Unfamiliar file format {}".format(os.path.splitext(args.output_file)[1])
        )

    logger.info("Starting S3 client")
    client = boto3.client("s3")
    paginator = client.get_paginator("list_objects")

    htseq_files = []
    log_files = []

    logger.info("Getting htseq file list")
    response_iterator = paginator.paginate(Bucket=args.s3_bucket, Prefix=args.s3_path)
    for result in response_iterator:
        if "Contents" in result:
            htseq_files.extend(
                r["Key"]
                for r in result["Contents"]
                if r["Key"].endswith("htseq-count.txt")
            )
            log_files.extend(
                r["Key"]
                for r in result["Contents"]
                if r["Key"].endswith("log.final.out")
            )
    logger.info("{} htseq files found".format(len(htseq_files)))

    sample_names = tuple(os.path.basename(fn)[:-16] for fn in htseq_files)

    gene_lists = set()
    gene_counts = list()

    for htseq_file in htseq_files:
        logger.debug("Downloading {}".format(htseq_file))
        if not dryrun:
            gene_list, gene_count = get_htseq_counts(client, args.s3_bucket, htseq_file)
            gene_lists.add(gene_list)
            gene_counts.append(gene_count)

    logger.info("Downloaded {} files".format(len(htseq_files)))
    if not dryrun:
        assert len(gene_lists) == 1
        gene_list = gene_lists.pop()

        logger.info("Writing to {}".format(args.output_file))

        with open(args.output_file, "w") as OUT:
            wtr = csv.writer(OUT, delimiter=sep)
            wtr.writerow(("gene",) + sample_names)
            for i, g in enumerate(gene_list):
                wtr.writerow((g,) + tuple(gc[i] for gc in gene_counts))

    if args.reflow:
        logger.info("Done!")
        return

    log_metrics = set()
    log_values = list()

    for log_file in log_files:
        logger.debug("Downloading {}".format(log_file))
        if not dryrun:
            metric_names, values = get_log_file(client, args.s3_bucket, log_file)
            log_metrics.add(metric_names)
            log_values.append(values)

    logger.info("Downloaded {} files".format(len(log_files)))
    if not dryrun:
        assert len(log_metrics) == 1
        log_metrics = log_metrics.pop()

    log_file = ".log".join(os.path.splitext(args.output_file))

    logger.info("Writing to {}".format(log_file))
    if not dryrun:
        with open(log_file, "w") as OUT:
            wtr = csv.writer(OUT, delimiter=sep)
            wtr.writerow(("metric",) + sample_names)
            for i, m in enumerate(log_metrics):
                wtr.writerow((m,) + tuple(mv[i] for mv in log_values))

    logger.info("Done!")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Construct the gene-cell table for an experiment\n"
            "e.g. gene_cell_table --s3_bucket czbiohub-maca"
            " data/170823_A00111_etc path/to/output.csv"
        ),
        epilog="See https://github.com/czbiohub/utilities for more examples",
        add_help=False,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # basic usage
    basic_group = parser.add_argument_group("basic arguments")
    basic_group.add_argument(
        "s3_path",
        help="Path to experiment. e.g. fastqs/171101_NB501961_0026_AHL33MBGX3",
    )
    basic_group.add_argument(
        "output_file", help="File to save the output, e.g. my_gc_table.csv"
    )

    # other arguments
    other_group = parser.add_argument_group("other options")
    other_group.add_argument(
        "--s3_bucket", help="S3 bucket. e.g. czb-seqbot", default="czb-seqbot"
    )
    other_group.add_argument(
        "--reflow", action="store_true", help="Output is from reflow, don't download logs"
    )
    other_group.add_argument(
        "--dryrun", action="store_true", help="Don't actually download any files"
    )
    other_group.add_argument(
        "--debug", action="store_true", help="Set logging to debug level"
    )
    other_group.add_argument(
        "-h", "--help", action="help", help="Show this help message and exit"
    )

    args = parser.parse_args()

    main_logger, _lf, _fh = get_logger(__name__, args.debug, args.dryrun)

    gene_cell_table(args, main_logger, args.dryrun)
