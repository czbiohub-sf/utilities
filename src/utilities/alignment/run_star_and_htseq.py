#!/usr/bin/env python
import argparse
import datetime
import logging
import os
import re
import subprocess
import tarfile

import multiprocessing as mp

from collections import defaultdict

import utilities.log_util as ut_log
import utilities.s3_util as s3u

import boto3
from boto3.s3.transfer import TransferConfig


S3_LOG_DIR = "s3://jamestwebber-logs/star_logs/"
S3_REFERENCE = {"east": "czi-hca", "west": "czbiohub-reference"}

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

    parser.add_argument("--taxon", choices=("homo", "mus"))

    parser.add_argument(
        "--s3_input_path",
        default="s3://czb-seqbot/fastqs",
        help="Location of input folders",
    )
    parser.add_argument(
        "--s3_output_path",
        default=None,
        help="Location for output, default [input_dir]/results",
    )
    parser.add_argument(
        "--region",
        default="east",
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
        "--input_dirs",
        nargs="+",
        required=True,
        help="List of input folders to process",
    )

    parser.add_argument(
        "--star_proc",
        type=int,
        default=8,
        help="Number of processes to give to each STAR run",
    )
    parser.add_argument(
        "--htseq_proc", type=int, default=4, help="Number of htseq processes to run"
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
        help=("Minimum file size (in bytes) for" " a sample to be aligned."),
    )

    return parser


def run_sample(
    star_queue, htseq_queue, log_queue, s3_input_bucket, genome_dir, run_dir, n_proc
):

    s3c = boto3.client("s3")
    t_config = TransferConfig(use_threads=False, num_download_attempts=25)

    for input_dir, sample_name, sample_fns in iter(star_queue.get, "STOP"):
        log_queue.put(("{} - {}".format(input_dir, sample_name), logging.INFO))
        dest_dir = os.path.join(run_dir, input_dir, sample_name)
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
            os.path.join(dest_dir, os.path.basename(sample_fn))
            for sample_fn in sample_fns
        )

        command = COMMON_PARS[:]
        command.extend(
            (
                "--runThreadN",
                str(n_proc),
                "--genomeDir",
                genome_dir,
                "--readFilesIn",
                " ".join(reads),
            )
        )
        failed = ut_log.log_command_to_queue(
            log_queue,
            command,
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
        failed = failed or ut_log.log_command_to_queue(
            log_queue, command, shell=True, cwd=os.path.join(dest_dir, "results")
        )

        # running samtools index -b
        command = [SAMTOOLS, "index", "-b", "Aligned.out.sorted.bam"]
        failed = failed or ut_log.log_command_to_queue(
            log_queue,
            command,
            shell=True,
            cwd=os.path.join(dest_dir, "results", "Pass1"),
        )

        # remove unsorted bam files
        if not failed:
            os.remove(os.path.join(dest_dir, "results", "Pass1", "Aligned.out.bam"))

        # remove fastq files
        for fastq_file in reads:
            os.remove(fastq_file)

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
        failed = failed or ut_log.log_command_to_queue(
            log_queue, command, shell=True, cwd=os.path.join(dest_dir, "results")
        )

        # ready to be htseq-ed and cleaned up
        if not failed:
            htseq_queue.put((input_dir, sample_name, dest_dir))


def run_htseq(htseq_queue, log_queue, s3_input_path, s3_output_path, taxon, sjdb_gtf):
    s3c = boto3.client("s3")
    t_config = TransferConfig(use_threads=False)

    for input_dir, sample_name, dest_dir in iter(htseq_queue.get, "STOP"):
        # running htseq
        command = [
            HTSEQ,
            "-r",
            "name",
            "-s",
            "no",
            "-f",
            "bam",
            "-m",
            "intersection-nonempty",
            os.path.join(dest_dir, "results", "Pass1", "Aligned.out.sorted-byname.bam"),
            sjdb_gtf,
            ">",
            "htseq-count.txt",
        ]
        failed = ut_log.log_command_to_queue(
            log_queue, command, shell=True, cwd=os.path.join(dest_dir, "results")
        )
        if failed:
            command = ["rm", "-rf", dest_dir]
            ut_log.log_command_to_queue(log_queue, command, shell=True)
            continue

        os.remove(
            os.path.join(dest_dir, "results", "Pass1", "Aligned.out.sorted-byname.bam")
        )

        # compress the results dir and move it to s3
        command = ["tar", "-cvzf", "{}.{}.tgz".format(sample_name, taxon), "results"]
        ut_log.log_command_to_queue(log_queue, command, shell=True, cwd=dest_dir)

        if s3_output_path is None:
            s3_output_path = os.path.join(s3_input_path, input_dir, "results")

        s3_output_bucket, s3_output_prefix = s3u.s3_bucket_and_key(s3_output_path)

        src_files = [
            os.path.join(dest_dir, "{}.{}.tgz".format(sample_name, taxon)),
            os.path.join(dest_dir, "results", "htseq-count.txt"),
            os.path.join(dest_dir, "results", "Pass1", "Log.final.out"),
            os.path.join(dest_dir, "results", "Pass1", "SJ.out.tab"),
            os.path.join(dest_dir, "results", "Pass1", "Aligned.out.sorted.bam"),
            os.path.join(dest_dir, "results", "Pass1", "Aligned.out.sorted.bam.bai"),
        ]

        dest_names = [
            "{}.{}.tgz".format(sample_name, taxon),
            "{}.{}.htseq-count.txt".format(sample_name, taxon),
            "{}.{}.log.final.out".format(sample_name, taxon),
            "{}.{}.SJ.out.tab".format(sample_name, taxon),
            "{}.{}.Aligned.out.sorted.bam".format(sample_name, taxon),
            "{}.{}.Aligned.out.sorted.bam.bai".format(sample_name, taxon),
        ]

        for src_file, dest_name in zip(src_files, dest_names):
            log_queue.put(("Uploading {}".format(dest_name), logging.INFO))
            s3c.upload_file(
                Filename=src_file,
                Bucket=s3_output_bucket,
                Key=os.path.join(s3_output_prefix, dest_name),
                Config=t_config,
            )

        # rm all the files
        command = ["rm", "-rf", dest_dir]
        ut_log.log_command_to_queue(log_queue, command, shell=True)


def main(logger):
    parser = get_parser()

    args = parser.parse_args()

    if os.environ.get("AWS_BATCH_JOB_ID"):
        root_dir = os.path.join("/mnt", os.environ["AWS_BATCH_JOB_ID"])
    else:
        root_dir = "/mnt"

    run_dir = os.path.join(root_dir, "data", "hca")
    os.makedirs(run_dir)

    if args.taxon == "homo":
        genome_dir = os.path.join(root_dir, "genome/STAR/HG38-PLUS/")
        ref_genome_file = "hg38-plus.tgz"
        ref_genome_star_file = "STAR/HG38-PLUS.tgz"
        sjdb_gtf = os.path.join(root_dir, "genome", "hg38-plus", "hg38-plus.gtf")
    elif args.taxon == "mus":
        genome_dir = os.path.join(root_dir, "genome/STAR/MM10-PLUS/")
        ref_genome_file = "mm10-plus.tgz"
        ref_genome_star_file = "STAR/MM10-PLUS.tgz"
        sjdb_gtf = os.path.join(root_dir, "genome", "mm10-plus", "mm10-plus.gtf")
    else:
        raise ValueError("Invalid taxon {}".format(args.taxon))

    if args.region == "east":
        ref_genome_file = os.path.join("ref-genome", ref_genome_file)
        ref_genome_star_file = os.path.join("ref-genome", ref_genome_star_file)

    if args.star_proc > mp.cpu_count():
        raise ValueError(
            "Not enough CPUs to give {} processes to STAR".format(args.star_proc)
        )

    s3_input_bucket, s3_input_prefix = s3u.s3_bucket_and_key(args.s3_input_path)

    logger.info(
        """Run Info: partition {} out of {}
                   star_proc:\t{}
                  htseq_proc:\t{}
                  genome_dir:\t{}
             ref_genome_file:\t{}
        ref_genome_star_file:\t{}
                    sjdb_gtf:\t{}
                       taxon:\t{}
               s3_input_path:\t{}
                  input_dirs:\t{}""".format(
            args.partition_id,
            args.num_partitions,
            args.star_proc,
            args.htseq_proc,
            genome_dir,
            ref_genome_file,
            ref_genome_star_file,
            sjdb_gtf,
            args.taxon,
            args.s3_input_path,
            ", ".join(args.input_dirs),
        )
    )

    s3 = boto3.resource("s3")

    # download the genome data
    os.mkdir(os.path.join(root_dir, "genome"))
    logger.info("Downloading and extracting genome data {}".format(ref_genome_file))

    s3_object = s3.Object(S3_REFERENCE[args.region], ref_genome_file)

    with tarfile.open(fileobj=s3_object.get()["Body"], mode="r|gz") as tf:
        tf.extractall(path=os.path.join(root_dir, "genome"))

    # download STAR stuff
    os.mkdir(os.path.join(root_dir, "genome", "STAR"))
    logger.info("Downloading and extracting STAR data {}".format(ref_genome_star_file))

    s3_object = s3.Object(S3_REFERENCE[args.region], ref_genome_star_file)

    with tarfile.open(fileobj=s3_object.get()["Body"], mode="r|gz") as tf:
        tf.extractall(path=os.path.join(root_dir, "genome", "STAR"))

    # Load Genome Into Memory
    command = [STAR, "--genomeDir", genome_dir, "--genomeLoad", "LoadAndExit"]
    ut_log.log_command(logger, command, shell=True)

    log_queue, log_thread = ut_log.get_thread_logger(logger)

    star_queue = mp.Queue()
    htseq_queue = mp.Queue()

    n_star_procs = mp.cpu_count() // args.star_proc

    star_args = (
        star_queue,
        htseq_queue,
        log_queue,
        s3_input_bucket,
        genome_dir,
        run_dir,
        args.star_proc,
    )
    star_procs = [
        mp.Process(target=run_sample, args=star_args) for i in range(n_star_procs)
    ]

    for p in star_procs:
        p.start()

    htseq_args = (
        htseq_queue,
        log_queue,
        args.s3_input_path,
        args.s3_output_path,
        args.taxon,
        sjdb_gtf,
    )
    htseq_procs = [
        mp.Process(target=run_htseq, args=htseq_args) for i in range(args.htseq_proc)
    ]

    for p in htseq_procs:
        p.start()

    sample_re = re.compile("([^/]+)_R\d_\d+.fastq.gz$")

    for input_dir in args.input_dirs:
        if args.s3_output_path is None:
            s3_output_path = os.path.join(args.s3_input_path, input_dir, "results")
        else:
            s3_output_path = args.s3_output_path

        s3_output_bucket, s3_output_prefix = s3u.s3_bucket_and_key(s3_output_path)

        # Check the input_dir folder for existing runs
        if not args.force_realign:
            output = s3u.prefix_gen(
                s3_output_bucket,
                s3_output_prefix,
                lambda r: (r["LastModified"], r["Key"]),
            )
        else:
            output = []

        output_files = {
            tuple(os.path.basename(fn).split(".")[:2])
            for dt, fn in output
            if fn.endswith("htseq-count.txt") and dt > CURR_MIN_VER
        }

        logger.info("Skipping {} existing results".format(len(output_files)))

        logger.info(
            "Running partition {} of {} for {}".format(
                args.partition_id, args.num_partitions, input_dir
            )
        )

        output = [
            (fn, s)
            for fn, s in s3u.get_size(
                s3_input_bucket, os.path.join(s3_input_prefix, input_dir)
            )
            if fn.endswith("fastq.gz")
        ]

        logger.info(f"number of fastq.gz files: {len(output)}")

        sample_lists = defaultdict(list)
        sample_sizes = defaultdict(list)

        for fn, s in output:
            matched = sample_re.search(os.path.basename(fn))
            if matched:
                sample_lists[matched.group(1)].append(fn)
                sample_sizes[matched.group(1)].append(s)

        for sample_name in sorted(sample_lists)[
            args.partition_id :: args.num_partitions
        ]:
            if (sample_name, args.taxon) in output_files:
                logger.info(f"{sample_name} already exists, skipping")
                continue

            if sum(sample_sizes[sample_name]) < args.min_size:
                logger.info(f"{sample_name} is below min_size, skipping")
                continue

            logger.info(f"Adding sample {sample_name} to queue")
            star_queue.put((input_dir, sample_name, sorted(sample_lists[sample_name])))

    for i in range(n_star_procs):
        star_queue.put("STOP")

    for p in star_procs:
        p.join()

    for i in range(args.htseq_proc):
        htseq_queue.put("STOP")

    for p in htseq_procs:
        p.join()

    log_queue.put("STOP")
    log_thread.join()

    # Remove Genome from Memory
    command = [STAR, "--genomeDir", genome_dir, "--genomeLoad", "Remove"]
    ut_log.log_command(logger, command, shell=True)

    logger.info("Job completed")


if __name__ == "__main__":
    mainlogger, log_file, file_handler = ut_log.get_logger(__name__)

    try:
        main(mainlogger)
    except:
        mainlogger.info("An exception occurred", exc_info=True)
        raise
    finally:
        if log_file:
            log_cmd = "aws s3 cp --quiet {} {}".format(log_file, S3_LOG_DIR)
            mainlogger.info(log_cmd)

            file_handler.close()
            subprocess.check_output(log_cmd, shell=True)
