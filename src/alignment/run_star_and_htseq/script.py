#!/usr/bin/env python

import datetime
import os
import re
import subprocess
import time

from collections import defaultdict

## VIASH START
par = {
    "input": "../resources_test/bs_195891710/fastqs/",
    "output": "../resources_test/bs_195891710/alignedqs/",
    "reference_genome": "../resources_test/reference_genome/",
    "taxon": "homo",
    "star_proc": 16,
    "partition_id": 0,
    "num_partitions": 10,
    "min_size": 50000,
}
## VIASH END

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
    "--outFilterType", "BySJout",
    "--outFilterMultimapNmax", "20",
    "--alignSJoverhangMin", "8",
    "--alignSJDBoverhangMin", "1",
    "--outFilterMismatchNmax", "999",
    "--outFilterMismatchNoverLmax", "0.04",
    "--alignIntronMin", "20",
    "--alignIntronMax", "1000000",
    "--alignMatesGapMax", "1000000",
    "--outSAMstrandField", "intronMotif",
    "--outSAMtype", "BAM", "Unsorted", 
    "--outSAMattributes", "NH", "HI", "NM", "MD",
    "--genomeLoad", "LoadAndKeep",
    "--outReadsUnmapped", "Fastx",
    "--readFilesCommand", "zcat",
]

CURR_MIN_VER = datetime.datetime(2017, 3, 1, tzinfo=datetime.timezone.utc)

def run_sample(
    input, sample_name, sample_fns, genome_dir, output, star_proc
):
    """ Run alignment jobs with STAR.

        s3_input_bucket - Input fastq files to align
        sample_name - Sequenced sample name (joined by "_")
        sample_fns - Sample file names. Each file name is concatenated by sample_name,
                     "_R1_" or"_R2_", a number, and ".fastq.gz"
        genome_dir - Path to reference genome
        run_dir - Path to output dir
        star_proc - Number of processes to give to each STAR run

        Return two values. FAILED is a boolean value of whether the alignment run
        fails. DEST_DIR is the path under which STAR alignment results are stored.
    """
    failed = False

    dest_dir = os.path.join(output, sample_name)

    if not os.path.exists(dest_dir):
        os.makedirs(dest_dir)
        os.mkdir(os.path.join(dest_dir, "rawdata"))
        os.mkdir(os.path.join(dest_dir, "results"))
        os.mkdir(os.path.join(dest_dir, "results", "Pass1"))

    # start running STAR
    # getting input files first
    reads = sorted(
        os.path.join(dest_dir, os.path.basename(sample_fn)) for sample_fn in sample_fns
    )

    input_command = COMMON_PARS[:]
    input_command.extend(
        (
            "--runThreadN", str(star_proc),
            "--genomeDir", genome_dir,
            "--readFilesIn",
            " ".join(reads),
        )
    )
    # run STAR
    with subprocess.Popen(
        input_command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    ) as p:
        for line in p.stdout:
            print(line.decode(), end="")

    if p.returncode > 0:
        print(f"STAR failed with exit code {p.returncode}")
        failed = True

    # run samtools sort
    sample_command = [
        SAMTOOLS,
        "sort",
        "-m", "6000000000",
        "-o",
        "./Pass1/Aligned.out.sorted.bam",
        "./Pass1/Aligned.out.bam",
    ]
    with subprocess.Popen(
        sample_command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        cwd=os.path.join(dest_dir, "results"),
    ) as p:
        for line in p.stdout:
            print(line.decode(), end="")

    if p.returncode > 0:
        print(f"samtools sort failed with exit code {p.returncode}")
        failed = True

    # run samtools index
    sample_index_command = [SAMTOOLS, "index", "-b", "Aligned.out.sorted.bam"]
    with subprocess.Popen(
        sample_index_command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        cwd=os.path.join(dest_dir, "results", "Pass1"),
    ) as p:
        for line in p.stdout:
            print(line.decode(), end="")

    if p.returncode > 0:
        print(f"samtools index failed with exit code {p.returncode}")
        failed = True

    # generating files for htseq-count
    # run samtools output
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
    with subprocess.Popen(
        output_command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        cwd=os.path.join(dest_dir, "results"),
    ) as p:
        for line in p.stdout:
            print(line.decode(), end="")

    if p.returncode > 0:
        print(f"samtools output failed with exit code {p.returncode}")
        failed = True

    return failed, dest_dir


def run_htseq(dest_dir, sjdb_gtf, id_attr):
    """ Run alignment job with htseq.

        dest_dir - Path local to the machine  under which alignment results
                   are stored. Child path of run_dir/sample_name
        sjdb_gtf - Path of reference genome .gtf files used to detect splice junctions
        id_attr - Determine naming format in the count file for different genomes

        Return FAILED, a boolean value of whether the alignment run fails
    """
    failed = False
    htseq_command = [
        HTSEQ,
        "-r", "name",
        "-s", "no",
        "-f", "bam",
        f"--idattr={id_attr}",
        "-m", "intersection-nonempty",
        os.path.join(dest_dir, "results", "Pass1", "Aligned.out.sorted-byname.bam"),
        sjdb_gtf,
        ">",
        "htseq-count.txt",
    ]
    with subprocess.Popen(
        htseq_command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        cwd=os.path.join(dest_dir, "results"),
    ) as p:
        for line in p.stdout:
            print(line.decode(), end="")

    if p.returncode > 0:
        print(f"htseq_command output failed with exit code {p.returncode}")
        failed = True

    return failed


""" Download reference genome, run alignment jobs, and upload results to S3.

    logger - Logger object that exposes the interface the code directly uses
"""
run_dir = os.path.join(par["output"])
os.makedirs(run_dir)

genome_name = reference_genomes[par["taxon"]]
if par["taxon"] == "gencode.vM19" or par["taxon"] == "gencode.vM19.ERCC":
    id_attr = "gene_name"
else:
    id_attr = "gene_id"
sjdb_gtf = par["reference_genome"] + f"/{genome_name}.gtf"

# Load Genome Into Memory
command = [STAR, "--genomeDir", par["reference_genome"], "--genomeLoad", "LoadAndExit"]
with subprocess.Popen(
    command,
    stdout=subprocess.PIPE,
    stderr=subprocess.STDOUT,
) as p:
    for line in p.stdout:
        print(line.decode(), end="")
print("Loaded reference genome into memory")

sample_re = re.compile("([^/]+)_R\d(?:_\d+)?.fastq.gz$")

output = []

output_files = {
    tuple(os.path.basename(fn).rsplit(".", 2)[0].split(".", 1)[:2])
    for dt, fn in output
    if fn.endswith(".htseq-count.txt") and dt > CURR_MIN_VER
}

sample_files = [
    (fn, os.path.getsize(fn))
    for fn in os.scandir(par["input"])
    if fn.endswith("fastq.gz")
]

sample_lists = defaultdict(list)
sample_sizes = defaultdict(list)

for fn, s in sample_files:
    matched = sample_re.search(os.path.basename(fn))
    if matched:
        sample_lists[matched.group(1)].append(fn)
        sample_sizes[matched.group(1)].append(s)

print(f"number of samples: {len(sample_lists)}")

for sample_name in sorted(sample_lists)[par["partition_id"] :: par["num_partitions"]]:
    if (sample_name, par["taxon"]) in output_files:
        continue

    if sum(sample_sizes[sample_name]) < par["min_size"]:
        continue

    failed, dest_dir = run_sample(
        par["input"],
        sample_name,
        sorted(sample_lists[sample_name]),
        par["reference_genome"],
        run_dir,
        par["star_proc"],
    )

    failed = failed or run_htseq(dest_dir, sjdb_gtf, id_attr)

    time.sleep(30)
