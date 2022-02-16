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

# reference genome bucket name for different regions
S3_REFERENCE = {"east": "czbiohub-reference-east", "west": "czbiohub-reference"}

# valid and deprecated reference genomes
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
    "homo.gencode.v30-plus-HAV18": "gencode.v30.annotation.ERCC92.HAV_18f_KP879216",
    "gencode_mouse_MTB": "gencode_mouse_MTB"
}
deprecated = {"homo": "hg38-plus", "mus": "mm10-plus"}

# other helpful constants
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

## VIASH START
par = {
    "taxon" : "homo",
    "input": "../resources_test/fastqs",
    "output": "../resources_test/output/",
    "reports": "../resources_test/bs_195891710/reports/",
}
## VIASH END

def run_sample(
    s3_input_bucket, sample_name, sample_fns, genome_dir, run_dir, star_proc, logger
):
    """ Run alignment jobs with STAR.

        s3_input_bucket - Name of the bucket with input fastq files to align
        sample_name - Sequenced sample name (joined by "_")
        sample_fns - Sample file names. Each file name is concatenated by sample_name,
                     "_R1_" or"_R2_", a number, and ".fastq.gz"
        genome_dir - Path to reference genome
        run_dir - Path local to the machine on EC2 under which alignment results
                  are stored before uploaded to S3
        star_proc - Number of processes to give to each STAR run
        logger - Logger object that exposes the interface the code directly uses

        Return two values. FAILED is a boolean value of whether the alignment run
        fails. DEST_DIR is the path under which STAR alignment results are stored.
        This path is local to the machine on EC2 running the alignment, and that's
        where we copy the alignment results to upload to S3 later.
    """

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

    input_command = COMMON_PARS[:]
    input_command.extend(
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
        input_command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        shell=True,
        cwd=os.path.join(dest_dir, "results", "Pass1"),
    )

    # running sam tools
    sample_command = [
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
        sample_command,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        cwd=os.path.join(dest_dir, "results"),
    )

    # running samtools index -b
    sample_index_command = [SAMTOOLS, "index", "-b", "Aligned.out.sorted.bam"]
    failed = failed or ut_log.log_command(
        logger,
        sample_index_command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        shell=True,
        cwd=os.path.join(dest_dir, "results", "Pass1"),
    )

    # generating files for htseq-count
    output_command = [
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
        output_command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        shell=True,
        cwd=os.path.join(dest_dir, "results"),
    )

    return failed, dest_dir


def run_htseq(dest_dir, sjdb_gtf, id_attr, logger):
    """ Run alignment job with htseq.

        dest_dir - Path local to the machine on EC2 under which alignment results
                   are stored before uploaded to S3. Child path of run_dir/sample_name
        sjdb_gtf - Path of reference genome .gtf files used to detect splice junctions
        id_attr - Determine naming format in the count file for different genomes
        logger - Logger object that exposes the interface the code directly uses

        Return FAILED, a boolean value of whether the alignment run fails
    """

    htseq_command = [
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
        htseq_command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        shell=True,
        cwd=os.path.join(dest_dir, "results"),
    )

    return failed


def upload_results(sample_name, taxon, dest_dir, s3_output_path, logger):
    """ Upload alignment results copied from EC2 machine directory onto S3.

        sample_name - Sequenced sample name (joined by "_")
        taxon - Reference genome name
        dest_dir - Path local to the machine on EC2 under which alignment results
                   are stored before uploaded to S3. Child path of run_dir/sample_name
        s3_output_path - S3 path of where the alignment results are stored
        logger - Logger object that exposes the interface the code directly uses
    """

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

run_dir = par["output"]
os.makedirs(run_dir)

# check if the input genome and region are valid
if par.taxon in reference_genomes:
    if par.taxon in deprecated:
        print(
            f"The name '{par.taxon}' will be removed in the future,"
            f" start using '{deprecated[par.taxon]}'"
        )

    genome_name = reference_genomes[par.taxon]
else:
    raise ValueError(f"unknown taxon {par.taxon}")

if par.taxon == "gencode.vM19" or par.taxon == "gencode.vM19.ERCC":
    id_attr = "gene_name"
else:
    id_attr = "gene_id"

genome_dir = os.path.join(par['output'], "genome", "STAR", genome_name)
ref_genome_star_file = f"STAR/{genome_name}.tgz"
sjdb_gtf = os.path.join(par['output'], f"{genome_name}.gtf")

print(
    f"""Run Info: partition {args.partition_id} out of {args.num_partitions}
                genome_dir:\t{genome_dir}
        ref_genome_star_file:\t{ref_genome_star_file}
                    sjdb_gtf:\t{sjdb_gtf}
                    id_attr:\t{id_attr}
                    taxon:\t{args.taxon}
            s3_input_path:\t{args.s3_input_path}"""
)


# download the reference genome data
os.mkdir(os.path.join(par["output"], "genome"))
print("Downloading and extracting gtf data {}".format(sjdb_gtf))

s3c.download_file(
    Bucket=S3_REFERENCE["west"],  # just always download this from us-west-2...
    Key=f"velocyto/{genome_name}.gtf",
    Filename=sjdb_gtf,
)

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