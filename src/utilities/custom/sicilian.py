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

STAR_params = [
    "STAR",
    "--runThreadN 4",
    "--twopassMode Basic",
    "--alignIntronMax 1000000",
    "--outSAMtype BAM Unsorted",
    "--outSAMattributes All",
    "--chimOutType WithinBAM SoftClip Junctions",
    "--chimJunctionOverhangMin 10",
    "--chimSegmentReadGapMax 0",
    "--chimOutJunctionFormat 1",
    "--chimSegmentMin 12",
    "--chimScoreJunctionNonGTAG -4",
    "--chimNonchimScoreDropMin 10",
    "--quantMode GeneCounts",
    "--outReadsUnmapped Fastx"
]

def get_default_requirements():
    """
    Define the default hardware requirements for this job
    """
    return argparse.Namespace(vcpus=4, memory=8000, storage=500)


def get_parser():
    """
    Define arguments
    """
    parser = argparse.ArgumentParser(
        prog="sicilian.py",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Run SICILIAN")

    # required arguments
    requiredNamed = parser.add_argument_group("required arguments")

    requiredNamed.add_argument(
        "--s3_input_path",
        required=True,
        help="Location of files to process with SICILIAN")

    requiredNamed.add_argument(
        "--s3_output_path",
        required=True,
        help="Location to store SICILIAN results")

    requiredNamed.add_argument(
        "--name",
        required=True,
        help="Name of sample files")

    requiredNamed.add_argument(
        "--bc_pattern",
        required=True,
        help="Pattern of barcode and UMI")

    requiredNamed.add_argument(
        "--single",
        required=True,
        help="Is single ended")

    requiredNamed.add_argument(
        "--tenX",
        required=True,
        help="Is tenX library")

    requiredNamed.add_argument(
        "--stranded_library",
        required=True,
        help="Is stranded library")

    requiredNamed.add_argument(
        "--R1_end",
        required=True,
        help="R1 file fastq extension")

    requiredNamed.add_argument(
        "--R2_end",
        required=True,
        help="R2 file fastq extension")

    requiredNamed.add_argument(
        "--gtf_file",
        required=True,
        help="GTF file")

    requiredNamed.add_argument(
        "--star_ref_path",
        required=True,
        help="path to STAR regerence genome")

    requiredNamed.add_argument(
        "--domain_file",
        required=True,
        help="Domain file")

    requiredNamed.add_argument(
        "--annotator_file",
        required=True,
        help="annotator_file")

    requiredNamed.add_argument(
        "--exon_pickle_file",
        required=True,
        help="exon_pickle_file")

    requiredNamed.add_argument(
        "--splice_pickle_file",
        required=True,
        help="splice_pickle_file")
    return parser


def whitelist(
    dest_dir, R1_file, bc_pattern, name):
    """
    Whitelist step
    """
    # Make directory to store whitelist results
    whitelist_dir = os.path.join(dest_dir, "whitelist")
    if not os.path.exists(whitelist_dir):
        os.makedirs(whitelist_dir)
    # umi tools input command
    input_command = [
        "umi_tools whitelist",
        "--stdin={}/{}".format(dest_dir, R1_file),
        "--bc-pattern={}".format(bc_pattern),
        "--stdout={}/{}_whitelist.txt".format(whitelist_dir, name),
        "--plot-prefix={}/{}".format(whitelist_dir, name)]
    # Run whitelist command
    failed = ut_log.log_command(
        input_command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        shell=True,
        cwd=whitelist_dir)
    return failed, whitelist_dir


def extract(
    dest_dir, whitelist_dir, R1_file, R2_file, bc_pattern, name, R1_end, R2_end):
    """
    Extract step
    """
    # Make directory to store extract results
    extract_dir = os.path.join(dest_dir, "extract")
    if not os.path.exists(extract_dir):
        os.makedirs(extract_dir)
    input_command = [
        "umi_tools extract ",
        "--error-correct-cell",
        "--bc-pattern {}".format(bc_pattern),
        "--stdin {}/{}".format(dest_dir, R1_file),
        "--stdout {}/{}_extracted{}".format(extract_dir, name, R1_end),
        "--read2-in {}/{}".format(dest_dir, R2_file),
        "--read2-out {}/{}_extracted{}".format(extract_dir, name, R2_end),
        "--whitelist {}/{}_whitelist.txt".format(whitelist_dir, name)]
    # Run whitelist command
    failed = ut_log.log_command(
        input_command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        shell=True,
        cwd=extract_dir)
    return failed, extract_dir


def map(
    dest_dir, star_ref_path, name, read_index, gtf_file, gzip, tenX, extract_dir, r_end):
    """
    Mapping step
    """
    # Make directory to store mapping results
    map_dir = os.path.join(dest_dir, "STAR_map")
    if not os.path.exists(map_dir):
        os.makedirs(map_dir)
    # Create mapping command
    input_command = STAR_params[:]
    input_command.extend(
        (
            "--genomeDir {}".format(star_ref_path),
            "--outFileNamePrefix {}/{}".format(map_dir, read_index+1),
            "--sjdbGTFfile {}".format(gtf_file),))
    if gzip:
        input_command.append("--readFilesCommand zcat")
    if tenX:
        input_command.append("--readFilesIn {}/{}_extracted{}".format(extract_dir, name, r_end))
    else:
        input_command.append("--readFilesIn {}/{}{}".format(dest_dir, name, r_end))
    failed = ut_log.log_command(
        input_command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        shell=True,
        cwd=map_dir)
    return failed, map_dir

def class_input(
    dest_dir, map_path, name, gtf_file, annotator_file, tenX, single, stranded_library):
    """
    Class input step
    """
    input_command = [
        "python3 scripts/light_class_input.py",
        "--outpath {}/".format(dest_dir),
        "--gtf {}".format(gtf_file),
        "--annotator {}".format(annotator_file),
        "--bams"]
    if single == 'True':
        input_command.append("{}/2Aligned.out.bam".format(map_path))
    else:
        input_command.append("{}/1Aligned.out.bam".format(map_path))
        input_command.append("{}/2Aligned.out.bam".format(map_path))
    if tenX == 'True':
        input_command.append("--stranded_library")
    if single != 'True':
        input_command.append("--paired")
    failed = ut_log.log_command(
        input_command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        shell=True,
        cwd=dest_dir)
    return failed


def GLM(
    dest_dir, name, gtf_file, tenX, single, stranded_library, domain_file, exon_pickle_file, splice_pickle_file):
    input_command = [
        "Rscript scripts/GLM_script_light.R",
        "{}/".format(dest_dir),
        gtf_file]
    if single == 'True':
        input_command.append("1")
    else:
        input_command.append("0")
    if tenX == 'True':
        input_command.append("1")
    else:
        input_command.append("0")
    if stranded_library == 'True':
        input_command.append("1")
    else:
        input_command.append("0")
    input_command.append("{} {} {}".format(domain_file, exon_pickle_file, splice_pickle_file))
    failed = ut_log.log_command(
            input_command,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            shell=True,
            cwd=dest_dir)
    return failed


def main():
    """
    sicilian --s3_input_path s3://salzman-lab/sicilian --s3_output_path s3://salzman-lab/sicilian_out --name sample --bc_pattern CCCCCCCCCCCCCCNNNNNNNNNNNN --R1_end _1.fastq --R2_end _2.fastq --single True --tenX True --star_ref_path /scratch/PI/horence/kaitlin/refs/hg38_STAR_2.7.5.a --gtf_file /scratch/PI/horence/kaitlin/covid/SZS_pipeline/util_files/GRCh38_latest_genomic.gtf --stranded_library True --domain_file domain_file.txt --annotator_file hg38_refseq.pkl  --exon_pickle_file hg38_refseq_exon_bounds.pkl --splice_pickle_file hg38_refseq_splices.pkl

    sicilian --s3_input_path s3://salzman-lab/sicilian --s3_output_path s3://salzman-lab/sicilian_out --name sample --bc_pattern CCCCCCCCCCCCCCNNNNNNNNNNNN --R1_end _1.fastq --R2_end _2.fastq --single True --tenX True --star_ref_path s3://salzman-lab/hg38_STAR_ref --gtf_file s3://salzman-lab/GRCh38_latest_genomic.gtf --stranded_library True --domain_file s3://salzman-lab/hg38_STAR_ref/domain_file.txt --annotator_file s3://salzman-lab/hg38_STAR_ref/hg38_refseq.pkl  --exon_pickle_file s3://salzman-lab/hg38_STAR_ref/hg38_refseq_exon_bounds.pkl --splice_pickle_file s3://salzman-lab/hg38_STAR_ref/hg38_refseq_splices.pkl
    """
    parser = get_parser()
    args = parser.parse_args()

    # Get input bucket name
    s3_input_bucket, s3_input_prefix = s3u.s3_bucket_and_key(args.s3_input_path)

    # Make destination directory
    dest_dir = os.path.join(os.getcwd(), args.s3_output_path, args.name)
    if not os.path.exists(dest_dir):
        os.makedirs(dest_dir)

    # boto3 download options
    s3 = boto3.client('s3')
    t_config = TransferConfig(use_threads=False, num_download_attempts=25)

    # Download R1 and R2 fastq files
    r_ends = [args.R1_end, args.R2_end]

    R1_file = "{}{}".format(args.name, args.R1_end)
    R2_file = "{}{}".format(args.name, args.R2_end)
    for sample_fn in [R1_file, R2_file]:
        s3.download_file(
            Bucket=s3_input_bucket, # Bucket = s3 bucket of file to download
            Key="sicilian/{}".format(sample_fn), #Key = file name(include path of s3 location ie 'sicilian/test.txt')
            Filename=os.path.join(dest_dir, os.path.basename(sample_fn)), # where to download to
            Config=t_config)

    if args.single:
        failed, whitelist_dir = whitelist(
            dest_dir,
            R1_file,
            args.bc_pattern,
            args.name)

        failed, extract_dir = extract(
            dest_dir,
            whitelist_dir,
            R1_file,
            R2_file,
            args.bc_pattern,
            args.name,
            args.R1_end,
            args.R2_end)

    if args.single == True:
        read_indices = 1
    else:
        read_indices = 0
    for read_index in range(read_indices, 2):
        r_end = r_ends[read_index]
        if r_end.split(".")[-1] == "gz":
            gzip = True
        else:
            gzip = False
        failed, map_dir = map(
            dest_dir,
            args.star_ref_path,
            args.name,
            read_index,
            args.gtf_file,
            gzip,
            args.tenX,
            extract_dir,
            r_ends[read_index]
        )

    failed = class_input(
        dest_dir,
        map_dir,
        args.name,
        args.gtf_file,
        args.annotator_file,
        args.tenX,
        args.single,
        args.stranded_library
    )

    failed = GLM(
        dest_dir,
        args.name,
        args.gtf_file,
        args.tenX,
        args.single,
        args.stranded_library,
        args.domain_file,
        args.exon_pickle_file,
        args.splice_pickle_file
    )


if __name__ == "__main__":
    main()
