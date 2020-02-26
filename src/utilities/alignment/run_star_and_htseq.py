#!/usr/bin/env python
import argparse
import datetime
import os
import re
import subprocess
import tarfile
import time

from collections import defaultdict

import utilities.log_util as ut_log
import utilities.s3_util as s3u

import boto3
from boto3.s3.transfer import TransferConfig


S3_REFERENCE = {"east": "czi-hca", "west": "czbiohub-reference"}

reference_genomes = {
    "homo": "HG38-PLUS",
    "hg38-plus": "HG38-PLUS",
    "homo.gencode.v30.ERCC.chrM": "homo.gencode.v30.annotation.ERCC92",
    "mus": "MM10-PLUS",
    "mm10-plus": "MM10-PLUS",
    "microcebus": "MicMur3-PLUS",
    "gencode.vM19": "gencode.vM19",
    "gencode.vM19.ERCC": "gencode.vM19.ERCC.SP1",
    "zebrafish-plus": "danio_rerio_plus_STAR2.6.1d",
    "homo.gencode.v30-plus-HAV18": "gencode.v30.annotation.ERCC92.HAV_18f_KP879216"
}

deprecated = {"homo": "hg38-plus", "mus": "mm10-plus"}

STAR = "STAR"
HTSEQ = "htseq-count"
SAMTOOLS = "samtools"

COMMON_PARS = [
    STAR,
    "--outFilterType",
    "BySJout",
    "--outFilterMultimapNmax",
    "20",
    "--alignSJoverhangMin",
    "8",
    "--alignSJDBoverhangMin",
    "1",
    "--outFilterMismatchNmax",
    "999",
    "--outFilterMismatchNoverLmax",
    "0.04",
    "--alignIntronMin",
    "20",
    "--alignIntronMax",
    "1000000",
    "--alignMatesGapMax",
    "1000000",
    "--outSAMstrandField",
    "intronMotif",
    "--outSAMtype",
    "BAM",
    "Unsorted",
    "--outSAMattributes",
    "NH",
    "HI",
    "NM",
    "MD",
    "--genomeLoad",
    "LoadAndKeep",
    "--outReadsUnmapped",
    "Fastx",
    "--readFilesCommand",
    "zcat",
]

CURR_MIN_VER = datetime.datetime(2017, 3, 1, tzinfo=datetime.timezone.utc)


def get_default_requirements():
    return argparse.Namespace(vcpus=16, memory=64000, storage=500, ecr_image="aligner")


def get_parser():
    parser = argparse.ArgumentParser(
        prog="run_star_and_htseq.py",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--taxon", required=True, choices=list(reference_genomes.keys())
    )

    parser.add_argument(
        "--s3_input_path", required=True, help="Location of input folder"
    )
    parser.add_argument(
        "--s3_output_path", required=True, help="Location for output files"
    )

    parser.add_argument(
        "--region",
        default="west",
        choices=("east", "west"),
        help=(
            "Region you're running jobs in."
            " Should match the location of"
            " the fastq.gz files"
        ),
    )

    parser.add_argument("--num_partitions", type=int, required=True)
    parser.add_argument("--partition_id", type=int, required=True)

    parser.add_argument(
        "--star_proc",
        type=int,
        default=16,
        help="Number of processes to give to each STAR run",
    )

    parser.add_argument(
        "--force_realign",
        action="store_true",
        help="Align files even when results already exist",
    )
    parser.add_argument(
        "--min_size",
        type=int,
        default=50000,
        help="Minimum file size (in bytes) for a sample to be aligned.",
    )

    return parser


def run_sample(
    s3_input_bucket, sample_name, sample_fns, genome_dir, run_dir, star_proc, logger
):
    t_config = TransferConfig(use_threads=False, num_download_attempts=25)

    dest_dir = os.path.join(run_dir, sample_name)

    if not os.path.exists(dest_dir):
        os.makedirs(dest_dir)
        os.mkdir(os.path.join(dest_dir, "rawdata"))
        os.mkdir(os.path.join(dest_dir, "results"))
        os.mkdir(os.path.join(dest_dir, "results", "Pass1"))

    for sample_fn in sample_fns:
        s3c.download_file(
            Bucket=s3_input_bucket,
            Key=sample_fn,
            Filename=os.path.join(dest_dir, os.path.basename(sample_fn)),
            Config=t_config,
        )

    # start running STAR
    # getting input files first

    reads = sorted(
        os.path.join(dest_dir, os.path.basename(sample_fn)) for sample_fn in sample_fns
    )

    command = COMMON_PARS[:]
    command.extend(
        (
            "--runThreadN",
            str(star_proc),
            "--genomeDir",
            genome_dir,
            "--readFilesIn",
            " ".join(reads),
        )
    )
    failed = ut_log.log_command(
        logger,
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        shell=True,
        cwd=os.path.join(dest_dir, "results", "Pass1"),
    )

    # running sam tools
    command = [
        SAMTOOLS,
        "sort",
        "-m",
        "6000000000",
        "-o",
        "./Pass1/Aligned.out.sorted.bam",
        "./Pass1/Aligned.out.bam",
    ]
    failed = failed or ut_log.log_command(
        logger,
        command,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        cwd=os.path.join(dest_dir, "results"),
    )

    # running samtools index -b
    command = [SAMTOOLS, "index", "-b", "Aligned.out.sorted.bam"]
    failed = failed or ut_log.log_command(
        logger,
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        shell=True,
        cwd=os.path.join(dest_dir, "results", "Pass1"),
    )

    # generating files for htseq-count
    command = [
        SAMTOOLS,
        "sort",
        "-m",
        "6000000000",
        "-n",
        "-o",
        "./Pass1/Aligned.out.sorted-byname.bam",
        "./Pass1/Aligned.out.sorted.bam",
    ]
    failed = failed or ut_log.log_command(
        logger,
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        shell=True,
        cwd=os.path.join(dest_dir, "results"),
    )

    return failed, dest_dir


def run_htseq(dest_dir, sjdb_gtf, id_attr, logger):
    command = [
        HTSEQ,
        "-r",
        "name",
        "-s",
        "no",
        "-f",
        "bam",
        f"--idattr={id_attr}",
        "-m",
        "intersection-nonempty",
        os.path.join(dest_dir, "results", "Pass1", "Aligned.out.sorted-byname.bam"),
        sjdb_gtf,
        ">",
        "htseq-count.txt",
    ]
    failed = ut_log.log_command(
        logger,
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        shell=True,
        cwd=os.path.join(dest_dir, "results"),
    )

    return failed


def upload_results(sample_name, taxon, dest_dir, s3_output_path, logger):
    t_config = TransferConfig(use_threads=False)

    s3_output_bucket, s3_output_prefix = s3u.s3_bucket_and_key(s3_output_path)

    src_files = [
        os.path.join(dest_dir, "results", "htseq-count.txt"),
        os.path.join(dest_dir, "results", "Pass1", "Log.final.out"),
        os.path.join(dest_dir, "results", "Pass1", "SJ.out.tab"),
        os.path.join(dest_dir, "results", "Pass1", "Aligned.out.sorted.bam"),
        os.path.join(dest_dir, "results", "Pass1", "Aligned.out.sorted.bam.bai"),
    ]

    dest_names = [
        "{}.{}.htseq-count.txt".format(sample_name, taxon),
        "{}.{}.log.final.out".format(sample_name, taxon),
        "{}.{}.SJ.out.tab".format(sample_name, taxon),
        "{}.{}.Aligned.out.sorted.bam".format(sample_name, taxon),
        "{}.{}.Aligned.out.sorted.bam.bai".format(sample_name, taxon),
    ]

    for src_file, dest_name in zip(src_files, dest_names):
        logger.info("Uploading {}".format(dest_name))
        s3c.upload_file(
            Filename=src_file,
            Bucket=s3_output_bucket,
            Key=os.path.join(s3_output_prefix, dest_name),
            Config=t_config,
        )


def main(logger):
    parser = get_parser()

    args = parser.parse_args()

    if os.environ.get("AWS_BATCH_JOB_ID"):
        root_dir = os.path.join("/mnt", os.environ["AWS_BATCH_JOB_ID"])
    else:
        root_dir = "/mnt"

    run_dir = os.path.join(root_dir, "data")
    os.makedirs(run_dir)

    if args.taxon in reference_genomes:
        if args.taxon in deprecated:
            logger.warn(
                f"The name '{args.taxon}' will be removed in the future,"
                f" start using '{deprecated[args.taxon]}'"
            )

        genome_name = reference_genomes[args.taxon]
    else:
        raise ValueError(f"unknown taxon {args.taxon}")

    if args.taxon == "gencode.vM19" or args.taxon == "gencode.vM19.ERCC":
        id_attr = "gene_name"
    else:
        id_attr = "gene_id"

    if args.region != "west" and genome_name not in ("HG38-PLUS", "MM10-PLUS"):
        raise ValueError(f"you must use --region west for {genome_name}")

    genome_dir = os.path.join(root_dir, "genome", "STAR", genome_name)
    ref_genome_star_file = f"STAR/{genome_name}.tgz"
    sjdb_gtf = os.path.join(root_dir, f"{genome_name}.gtf")

    if args.region == "east":
        ref_genome_star_file = os.path.join("ref-genome", ref_genome_star_file)

    s3_input_bucket, s3_input_prefix = s3u.s3_bucket_and_key(args.s3_input_path)

    logger.info(
        f"""Run Info: partition {args.partition_id} out of {args.num_partitions}
                   genome_dir:\t{genome_dir}
         ref_genome_star_file:\t{ref_genome_star_file}
                     sjdb_gtf:\t{sjdb_gtf}
                      id_attr:\t{id_attr}
                        taxon:\t{args.taxon}
                s3_input_path:\t{args.s3_input_path}"""
    )

    s3 = boto3.resource("s3")

    # download the gtf file
    os.mkdir(os.path.join(root_dir, "genome"))
    logger.info("Downloading and extracting gtf data {}".format(sjdb_gtf))

    s3c.download_file(
        Bucket=S3_REFERENCE["west"],  # just always download this from us-west-2...
        Key=f"velocyto/{genome_name}.gtf",
        Filename=sjdb_gtf,
    )

    # download STAR stuff
    os.mkdir(os.path.join(root_dir, "genome", "STAR"))
    logger.info("Downloading and extracting STAR data {}".format(ref_genome_star_file))

    s3_object = s3.Object(S3_REFERENCE[args.region], ref_genome_star_file)

    with tarfile.open(fileobj=s3_object.get()["Body"], mode="r|gz") as tf:
        tf.extractall(path=os.path.join(root_dir, "genome", "STAR"))

    # Load Genome Into Memory
    command = [STAR, "--genomeDir", genome_dir, "--genomeLoad", "LoadAndExit"]
    if ut_log.log_command(
        logger, command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True
    ):
        raise RuntimeError("Failed to load genome into memory")

    sample_re = re.compile("([^/]+)_R\d(?:_\d+)?.fastq.gz$")
    s3_output_bucket, s3_output_prefix = s3u.s3_bucket_and_key(args.s3_output_path)

    logger.info(
        "Running partition {} of {}".format(args.partition_id, args.num_partitions)
    )

    # Check the input folder for existing runs
    if not args.force_realign:
        output = s3u.prefix_gen(
            s3_output_bucket, s3_output_prefix, lambda r: (r["LastModified"], r["Key"])
        )
    else:
        output = []

    output_files = {
        tuple(os.path.basename(fn).rsplit(".", 2)[0].split(".", 1)[:2])
        for dt, fn in output
        if fn.endswith(".htseq-count.txt") and dt > CURR_MIN_VER
    }

    logger.info("Skipping {} existing results".format(len(output_files)))

    sample_files = [
        (fn, s)
        for fn, s in s3u.get_size(s3_input_bucket, s3_input_prefix)
        if fn.endswith("fastq.gz")
    ]

    sample_lists = defaultdict(list)
    sample_sizes = defaultdict(list)

    for fn, s in sample_files:
        matched = sample_re.search(os.path.basename(fn))
        if matched:
            sample_lists[matched.group(1)].append(fn)
            sample_sizes[matched.group(1)].append(s)

    logger.info(f"number of samples: {len(sample_lists)}")

    for sample_name in sorted(sample_lists)[args.partition_id :: args.num_partitions]:
        if (sample_name, args.taxon) in output_files:
            logger.debug(f"{sample_name} already exists, skipping")
            continue

        if sum(sample_sizes[sample_name]) < args.min_size:
            logger.info(f"{sample_name} is below min_size, skipping")
            continue

        failed, dest_dir = run_sample(
            s3_input_bucket,
            sample_name,
            sorted(sample_lists[sample_name]),
            genome_dir,
            run_dir,
            args.star_proc,
            logger,
        )

        failed = failed or run_htseq(dest_dir, sjdb_gtf, id_attr, logger)

        if not failed:
            upload_results(
                sample_name, args.taxon, dest_dir, args.s3_output_path, logger
            )

        command = ["rm", "-rf", dest_dir]
        ut_log.log_command(logger, command, shell=True)

        time.sleep(30)

    logger.info("Job completed")


if __name__ == "__main__":
    mainlogger, log_file, file_handler = ut_log.get_logger(__name__)
    s3c = boto3.client("s3")
    main(mainlogger)
