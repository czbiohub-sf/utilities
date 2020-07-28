#!/usr/bin/env python
import argparse
import datetime
import os
import pathlib
import posixpath
import re
import subprocess
import tarfile
import time
import sys
import shutil
import glob

from collections import defaultdict

import utilities.log_util as ut_log
import utilities.s3_util as s3u

import boto3
from boto3.s3.transfer import TransferConfig


def get_default_requirements():
    return argparse.Namespace(
        vcpus=16, memory=64000, storage=500, ecr_image="velocyto"
    )

def display_info(logger):
    """Displays kb, kallisto and bustools version + citation information, along
    with a brief description and examples.
    
    Keyword Argument:
    mainlogger - Logger of main function (type: logging.Logger)
    """
    info_command = ["kb", "info"]
    
#     if ut_log.log_command(
#         logger,
#         info_command,
#         stdout=subprocess.PIPE,
#         stderr=subprocess.STDOUT,
#         shell=True,
#     ):
#         logger.info("Failed to view kb_python package details")
#     sys.exit(1)
    
    proc = subprocess.run(" ".join(info_command), **kwargs)

    if proc.returncode != 0:
        raise RuntimeError("`info` command failed")
        if proc.stdout and isinstance(proc.stdout, str):
            raise RuntimeError(proc.stdout)
        elif isinstance(proc.stdout, bytes):
            raise RuntimeError(proc.stdout.decode())

        return True
    else:
        return False


def display_technologies(logger):
    """Displays a list of supported technologies along with whether kb provides
    a whitelist for that technology and the FASTQ argument order for kb count.
    
    Keyword Argument:
    mainlogger - Logger of main function (type: logging.Logger)
    """
    technology_command = ["kb", "--list"]
    
    if ut_log.log_command(
        logger,
        technology_command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        shell=True,
    ):
        logger.info("Failed to view single-cell technology list compatible with the kb_python package") 
    sys.exit(1)

#     proc = subprocess.run(" ".join(technology_command), **kwargs)

#     if proc.returncode != 0:
#         raise RuntimeError("`--list` command failed")
#         if proc.stdout and isinstance(proc.stdout, str):
#             raise RuntimeError(proc.stdout)
#         elif isinstance(proc.stdout, bytes):
#             raise RuntimeError(proc.stdout.decode())

#         return True
#     else:
#         return False

    
def parse_ref(args, run_dir, logger):
    """Parser for the `ref` command. Build kallisto index files.
    
    Keyword Arguments:
    args -- Command-line arguments dictionary, as parsed by argparse (type: dict)
    run_dir -- Path on the EC2 instance under which jobs are run (type: pathlib.Path)
    logger -- Logger object that exposes the interface the code directly uses (type: logging.Logger)
    """
    # Build paths on the EC2 instance to store inputs and outputs of kallisto index building. Download kallisto index building files from the AWS S3 bucket to the EC2 instance, if not using the in-built kallisto indices
    kallisto_index_dir = run_dir / "kallisto_index"
    kallisto_index_inputs = kallisto_index_dir / "inputs"
    kallisto_index_outputs = kallisto_index_dir / "outputs"
    kallisto_index_inputs.mkdir(parents=True)
    kallisto_index_outputs.mkdir(parents=True)
    kb_ref_paths = dict()
    s3_kb_ref = dict()
    kb_ref_output_to_s3 = dict()

    if "-d" not in sys.argv:
        if "--workflow" in sys.argv and "kite" in sys.argv:
            kb_ref_paths["feature_path"] = kallisto_index_inputs / os.path.basename(args.feature)
            s3_kb_ref["s3_feature_bucket"] = s3u.s3_bucket_and_key(args.feature)[0]
            s3_kb_ref["s3_feature_prefix"] = s3u.s3_bucket_and_key(args.feature)[1]
            s3c.download_file(
                Bucket=s3_kb_ref["s3_feature_bucket"],
                Key=s3_kb_ref["s3_feature_prefix"],
                Filename=str(kb_ref_paths["feature_path"]),
            )
        for arg in ["fasta", "gtf"]:
            print(f"testing purpose - see if args.arg output the values of fasta or gtf this time. change the format of all args.arg if this works: {args.arg[1:-1]}") # testing purpose
            kb_ref_paths[arg + "_path"] = kallisto_index_inputs / os.path.basename(args.arg)
            s3_kb_ref["s3_" + arg + "_bucket"] = s3u.s3_bucket_and_key(args.arg)[0]
            s3_kb_ref["s3_" + arg + "_prefix"] = s3u.s3_bucket_and_key(args.arg)[1]
            s3c.download_file(
                Bucket=s3_kb_ref["s3_" + arg + "_bucket"],
                Key=s3_kb_ref["s3_" + arg + "_prefix"],
                Filename=str(kb_ref_paths[arg + "_path"]),
            )
            
    for arg in ["-i", "-g", "-f1", "-f2", "-c1", "-c2", "--tmp"]:
        if arg in sys.argv:
            arg = arg[2:] if arg == "--tmp" else arg[1:]
            print(f"testing purpose - see if args.arg output all argument names from '-i', '-g', '-f1', '-f2', '-c1', '-c2', '--tmp': {args.arg}") # testing purpose
            if arg == "tmp": 
                kb_ref_paths[arg + "_path"] = kallisto_index_dir / "alter_tmp"
                s3_kb_ref["s3_" + arg + "_bucket"] = s3u.s3_bucket_and_key(args.arg)[0]
                s3_kb_ref["s3_" + arg + "_prefix"] = posixpath.join(s3u.s3_bucket_and_key(args.arg)[1], "alter_tmp")
            else:
                kb_ref_paths[arg + "_path"] = kallisto_index_outputs / os.path.basename(args.arg)
                s3_kb_ref["s3_" + arg + "_bucket"] = s3u.s3_bucket_and_key(args.arg)[0]
                s3_kb_ref["s3_" + arg + "_prefix"] = s3u.s3_bucket_and_key(args.arg)[1]


    # Build the command of running `kb ref` to generate kallisto index files
    ref_input_boolean = ['--lamanno', '--overwrite', '--keep-tmp', '--verbose']
    ref_input_upload_required = ['-i', '-g', '-f1', '-f2', '-c1', '-c2', '--tmp']
    ref_input_left_args = ['-d', '-n', '-k', '--workflow']

    kb_ref_command = ['kb', 'ref']
    for input in ["fasta", "gtf"]:
        if "-d" not in sys.argv:
            print(f"testing purpose: `ref` positional argument: {args.input}") # testing purpose
            if "--workflow" in sys.argv and "kite" in sys.argv:
                print(f"testing purpose: `ref` positional argument: {args.input}") # testing purpose
                kb_ref_command += [str(kb_ref_paths["feature_path"])]
            kb_ref_command += [str(kb_ref_paths[input + "_path"])]
    for input in ref_input_boolean:
        if input in sys.argv:
            kb_ref_command += [input]
    for input in ref_input_upload_required:
        if input in sys.argv:
            kb_ref_command += [input, str(kb_ref_paths[input[2:] + "_path"]) if input == "--tmp" else str(kb_ref_paths[input[1:] + "_path"])]
    for input in ref_input_left_args:
        if input in sys.argv:
            kb_ref_command += [input, args.input]
    print(f"testing purpose - check if kb_ref_command is of correct format: {kb_ref_command}") # testing purpose

    # Run the command to generate kallisto index files, and upload the output files on the EC2 instance back to AWS S3 
    kb_ref_failed = ut_log.log_command(
        logger,
        kb_ref_command,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        universal_newlines=True,
    )
    kb_ref_t_config = TransferConfig(use_threads=False) # testing purpose: comment this line if this runs into error. kallisto indices are not pretty big and don't necessarily need a transfer config
    if kb_ref_failed:
        raise RuntimeError("kallisto index building failed")
    else:
        kb_ref_upload_files = glob.glob(str(kallisto_index_outputs / "*"))
        print(f"testing purpose - check kallisto index files: {kb_ref_upload_files}") # testing purpose
        for file in kb_ref_upload_files:
            logger.info(f"Uploading {os.path.basename(file)} from the EC2 instance to AWS S3")
            s3c.upload_file(
                Filename=file,
                Bucket=kb_ref_output_to_s3[file][0],
                Key=kb_ref_output_to_s3[file][1],
                Config=kb_ref_t_config,
            )
            time.sleep(30)
            print(f"testing purpose - see if kb_ref upload_file function intakes correct bucket and prefix names for kallisto index output files: {kb_ref_output_to_s3[file][0]}, {kb_ref_output_to_s3[file][1]}") # testing purpose
    return


def parse_count(args, run_dir, logger):
    """Parser for the `count` command. Data quantification with kallisto and bustools.
    
    Keyword Arguments:
    args -- Command-line arguments dictionary, as parsed by argparse (type: dict)
    run_dir -- Path on the EC2 instance under which jobs are run (type: pathlib.Path)
    logger -- Logger object that exposes the interface the code directly uses (type: logging.Logger)
    """
    # Build paths on the EC2 instance to store inputs and outputs of kb data quantification results.
    kb_count_dir = run_dir / "kb_count"
    kb_fastqs = kb_count_dir / "fastqs"
    kb_count_inputs = kb_count_dir / "inputs"
    kb_count_outputs = kb_count_dir / "outputs"
    kb_fastqs.mkdir(parents=True)
    kb_count_outputs.mkdir(parents=True)
    kb_count_paths = dict()
    s3_kb_count = dict()

    for arg in ["--tmp", "-o", "-w", "-i", "-g", "-c1", "-c2"]:
        if arg in sys.argv:
            arg = arg[2:] if arg == "--tmp" else arg[1:]
            print(f"testing purpose - see if args.arg returns the correct values for `kb count` inputs: {args.arg}") # testing purpose
            if arg == "tmp": 
                kb_count_paths[arg + "_path"] = kb_count_dir / "alter_tmp"
                s3_kb_count["s3_" + arg + "_prefix"] = posixpath.join(s3u.s3_bucket_and_key(args.arg)[1], "alter_tmp")
            elif arg == "o":
                kb_count_paths[arg + "_path"] = kb_count_outputs
                s3_kb_count["s3_" + arg + "_prefix"] = posixpath.join(s3u.s3_bucket_and_key(args.arg)[1], "outputs")
            else:
                kb_count_paths[arg + "_path"] = kb_count_inputs / os.path.basename(args.arg)
                s3_kb_count["s3_" + arg + "_prefix"] = s3u.s3_bucket_and_key(args.arg)[1]
                s3_kb_count["s3_" + arg + "_bucket"] = s3u.s3_bucket_and_key(args.arg)[0]
                s3c.download_file(
                    Bucket=s3_kb_count["s3_" + arg + "_bucket"],
                    Key=s3_kb_count["s3_" + arg + "_prefix"],
                    Filename=str(kb_count_paths[arg + "_path"]),
                )
            s3_kb_count["s3_" + arg + "_bucket"] = s3u.s3_bucket_and_key(args.arg)[0]
                
                
                

    # Download fastq files from the AWS S3 bucket to the EC2 instance.
    kb_count_paths["fastqs_path"], s3_kb_count["s3_fastqs_bucket"], s3_kb_count["s3_fastqs_prefix"] = dict(), dict(), dict()
    kb_count_fastq_files_paths, s3_kb_count_fastqs_bucket, s3_kb_count_fastqs_prefix = kb_count_paths["fastqs_path"], s3_kb_count["s3_fastqs_bucket"], s3_kb_count["s3_fastqs_prefix"]
    
    s3_fastq_folder_bucket, s3_fastq_folder_prefix = s3u.s3_bucket_and_key(args.fastq_folder)
    s3_fastq_files_prefix = list(s3u.get_files(bucket=s3_fastq_folder_bucket, prefix=s3_fastq_folder_prefix))[1:]
    print("testing purpose - see if all fastq files prefix are extracted: {s3_fastq_files_prefix}") # testing purpose
    fastq_format = re.compile("([^/]+)_R\d(?:_\d+)?.fastq.gz$")

    for fastq_prefix in s3_fastq_files_prefix:
        if not fastq_format.search(os.path.basename(fastq_prefix)):
            continue
        kb_count_fastq_files_paths[os.path.basename(fastq_prefix)] = kb_fastqs / os.path.basename(fastq_prefix)
        s3_kb_count_fastqs_bucket[os.path.basename(fastq_prefix)] = s3_fastq_folder_bucket
        s3_kb_count_fastqs_prefix[os.path.basename(fastq_prefix)] = fastq_prefix
        fastq_t_config = TransferConfig(use_threads=False) # testing purpose: comment this line if this runs into error.
        s3c.download_file(
            Bucket=s3_kb_count_fastqs_bucket[os.path.basename(fastq_prefix)],
            Key=s3_kb_count_fastqs_prefix[os.path.basename(fastq_prefix)],
            Filename=str(kb_count_fastq_files_paths[os.path.basename(fastq_prefix)]),
            Config=fastq_t_config,
        )

    # Build the command of running `kb count` to generate count matrices        
    count_input_boolean = ['--keep-tmp', '--verbose', '--mm', '--tcc', '--overwrite', '--lamanno', '--nucleus', '--loom', '--h5ad']
    count_input_file_transfer_required = ['--tmp', '-o', '-w']
    count_input_kb_indices = ['-i', '-g', '-c1', '-c2']
    count_input_left_args = ['-t', '-m', '--workflow', '--filter', '-x']

    kb_count_command = ['kb', 'count']
    for fastq_path in kb_count_fastq_files_paths.values():
        print(f"testing purpose - view the paths of individual fastqs on the EC2 instance, i.e. values of dictionary `kb_count_fastq_files_paths`: {kb_count_fastq_files_paths.values()}") # testing purpose
        kb_count_command += [str(fastq_path)]
    for input in count_input_boolean:
        if input in sys.argv:
            kb_count_command += [input]
    for input in count_input_file_transfer_required:
        if input in sys.argv:
            kb_count_command += [input, str(kb_count_paths[input[2:] + "_path"]) if input == "--tmp" else str(kb_count_paths[input[1:] + "_path"])]
    for input in count_input_kb_indices:
        if input in sys.argv:
            kb_count_command += [input, str(kb_count_paths[input[1:] + "_path"])]
    for input in count_input_left_args:
        if input in sys.argv:
            kb_count_command += [input, args.input]
    print(f"testing purpose - view kb count command: {kb_count_command}") # testing purpose

    # Run the command to generate count matrices, and upload the output files on the EC2 instance back to AWS S3 
    kb_count_failed = ut_log.log_command(
        logger,
        kb_count_command,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        universal_newlines=True,
    )
    kb_count_t_config = TransferConfig(use_threads=False) # testing purpose: comment this line if this runs into error. count matrices are not pretty big and don't necessarily need a transfer config
    if failed:
        raise RuntimeError("kb data quantification failed")
    else:
        logger.info(f"Uploading kb quantification data from the EC2 instance to AWS S3")
        for root,dirs,files in os.walk(str(kb_count_outputs)):
            for file in files:
                s3c.upload_file(
                    Filename=os.path.join(root,file),
                    Bucket=s3_kb_count["s3_o_bucket"],
                    Key=posixpath.join(s3_kb_count["s3_o_prefix"], file),
                    Config=kb_count_t_config,
                )
                time.sleep(30)
                print(f"testing purpose - see if kb data quantification outputs are uploaded to the correct s3 paths: " + s3_kb_count["s3_o_bucket"] + ", " + posixpath.join(s3_kb_count["s3_o_prefix"], file)) # testing purpose
    return


COMMAND_TO_FUNCTION = {
    'ref': parse_ref,
    'count': parse_count,
}


def setup_info_args(parser, parent):
    """Helper function to set up a subparser for the `info` command.
    
    Keyword Arguments:
    parser -- argparse parser to add the `info` command to (type: argparse.ArgumentParser)
    parent -- argparse parser parent of the newly added subcommand. used to inherit shared commands/flags (type: argparse.ArgumentParser)

    Return:
    the newly added parser (type: argparse.ArgumentParser)
    """
    parser_info = parser.add_parser(
        'info',
        description='Display kb-python package and citation information',
        help='Display kb-python package and citation information',
        parents=[parent],
        add_help=False,
    )
    return parser_info


def setup_ref_args(parser, parent):
    """Helper function to set up a subparser for the `ref` command.
    
    Keyword Arguments:
    parser -- argparse parser to add the `ref` command to (type: argparse.ArgumentParser)
    parent -- argparse parser parent of the newly added subcommand. used to inherit shared commands/flags (type: argparse.ArgumentParser)
    
    Return:
    the newly added parser (type: argparse.ArgumentParser)
    """
    workflow = sys.argv[sys.argv.index('--workflow') +
                        1] if '--workflow' in sys.argv else 'standard'

    parser_ref = parser.add_parser(
        'ref',
        description='Build a kallisto index and transcript-to-gene mapping',
        help='Build a kallisto index and transcript-to-gene mapping',
        parents=[parent],
    )
    parser_ref._actions[0].help = parser_ref._actions[0].help.capitalize()

    # required arguments
    required_ref = parser_ref.add_argument_group('required arguments')
    required_ref.add_argument(
        '-i',
        metavar='INDEX',
        help='Path to the kallisto index to be constructed',
        type=str,
        required=True
    )
    required_ref.add_argument(
        '-g',
        metavar='T2G',
        help='Path to transcript-to-gene mapping to be generated',
        type=str,
        required=True
    )
    required_ref.add_argument(
        '-f1',
        metavar='FASTA',
        help=(
            '[Optional with -d] Path to the cDNA FASTA (lamanno, nucleus) '
            'or mismatch FASTA (kite) to be generated '
        ),
        type=str,
        required='-d' not in sys.argv
    )
    
    # required arguments for lamanno and nucleus workflows
    required_lamanno = parser_ref.add_argument_group(
        'required arguments for `lamanno` and `nucleus` workflows'
    )
    required_lamanno.add_argument(
        '-f2',
        metavar='FASTA',
        help='Path to the intron FASTA to be generated',
        type=str,
        required=workflow in {'lamanno', 'nucleus'}
        or any(arg in sys.argv for arg in {'--lamanno', '--nucleus'})
    )
    required_lamanno.add_argument(
        '-c1',
        metavar='T2C',
        help='Path to generate cDNA transcripts-to-capture',
        type=str,
        required=workflow in {'lamanno', 'nucleus'}
        or any(arg in sys.argv for arg in {'--lamanno', '--nucleus'})
    )
    required_lamanno.add_argument(
        '-c2',
        metavar='T2C',
        help='Path to generate intron transcripts-to-capture',
        type=str,
        required=workflow in {'lamanno', 'nucleus'}
        or any(arg in sys.argv for arg in {'--lamanno', '--nucleus'})
    )

    # other optional and conditionally required arguments
    parser_ref.add_argument(
        '-d',
        help=(
            'Download a pre-built kallisto index (along with all necessary files) '
            'instead of building it locally'
        ),
        type=str,
        choices=['human', 'mouse', 'linnarsson'],
        required=False
    )
    parser_ref.add_argument(
        '--lamanno',
        help='Deprecated. Use `--workflow lamanno` instead.',
        action='store_true'
    )
    parser_ref.add_argument(
        '--overwrite',
        help='Overwrite existing kallisto index',
        action='store_true'
    )
    parser_ref.add_argument(
        'fasta',
        help='Genomic FASTA file',
        type=str,
        nargs=None if '-d' not in sys.argv and workflow != 'kite' else '?'
    )
    parser_ref.add_argument(
        'gtf',
        help='Reference GTF file',
        type=str,
        nargs=None if '-d' not in sys.argv and workflow != 'kite' else '?'
    )
    parser_ref.add_argument(
        'feature',
        help=(
            '[`kite` workflow only] Path to TSV containing barcodes and feature names.'
        ),
        type=str,
        nargs=None if '-d' not in sys.argv and workflow == 'kite' else '?'
    )
    return parser_ref


def setup_count_args(parser, parent):
    """Helper function to set up a subparser for the `count` command.
    
    Keyword Arguments:
    parser -- argparse parser to add the `count` command to (type: argparse.ArgumentParser)
    parent -- argparse parser parent of the newly added subcommand. used to inherit shared commands/flags (type: argparse.ArgumentParser)
    
    Return:
    the newly added parser (type: argparse.ArgumentParser)
    """
    workflow = sys.argv[sys.argv.index('--workflow') +
                        1] if '--workflow' in sys.argv else 'standard'

    # count
    parser_count = parser.add_parser(
        'count',
        description=('Generate count matrices from a set of single-cell FASTQ files. '
                     'Run `--list` to view single-cell technology information.'),  # noqa
        help='Generate count matrices from a set of single-cell FASTQ files',
        parents=[parent],
    )
    parser_count._actions[0].help = parser_count._actions[0].help.capitalize()
    
    # positional argument
    parser_count.add_argument('fastq_folder', help='Path to the folder containing FASTQ files to be processed')
    
    # required arguments
    required_count = parser_count.add_argument_group('required arguments')
    required_count.add_argument(
        '-i',
        metavar='INDEX',
        help='Path to kallisto index',
        type=str,
        required=True
    )
    required_count.add_argument(
        '-g',
        metavar='T2G',
        help='Path to transcript-to-gene mapping',
        type=str,
        required=True
    )
    required_count.add_argument(
        '-x',
        metavar='TECHNOLOGY',
        help='Single-cell technology used (`--list` to view)',
        type=str,
        required=True
    )
    
    # optional arguments
    parser_count.add_argument(
        '-o',
        metavar='OUT',
        help='Path to output directory (default: current directory)',
        type=str,
        default='.',
    )
    parser_count.add_argument(
        '-w',
        metavar='WHITELIST',
        help=(
            'Path to file of whitelisted barcodes to correct to. '
            'If not provided and bustools supports the technology, '
            'a pre-packaged whitelist is used. If not, the bustools '
            'whitelist command is used. (`--list` to view whitelists)'
        ),
        type=str
    )
    parser_count.add_argument(
        '-t',
        metavar='THREADS',
        help='Number of threads to use (default: 8)',
        type=int,
        default=8
    )
    parser_count.add_argument(
        '-m',
        metavar='MEMORY',
        help='Maximum memory used (default: 4G)',
        type=str,
        default='4G'
    )
    parser_count.add_argument(
        '--tcc',
        help='Generate a TCC matrix instead of a gene count matrix.',
        action='store_true'
    )
    parser_count.add_argument(
        '--overwrite',
        help='Overwrite existing output.bus file',
        action='store_true'
    )
    parser_count.add_argument(
        '--filter',
        help='Produce a filtered gene count matrix (default: bustools)',
        type=str,
        const='bustools',
        nargs='?',
        choices=['bustools']
    )
    
    # required arguments for lamanno and nucleus workflows
    required_lamanno = parser_count.add_argument_group(
        'required arguments for `lamanno` and `nucleus` workflows'
    )
    required_lamanno.add_argument(
        '-c1',
        metavar='T2C',
        help='Path to cDNA transcripts-to-capture',
        type=str,
        required=workflow in {'lamanno', 'nucleus'}
        or any(arg in sys.argv for arg in {'--lamanno', '--nucleus'})
    )
    required_lamanno.add_argument(
        '-c2',
        metavar='T2C',
        help='Path to intron transcripts-to-captured',
        type=str,
        required=workflow in {'lamanno', 'nucleus'}
        or any(arg in sys.argv for arg in {'--lamanno', '--nucleus'})
    )

    # only one workflow (i.e. lamanno or nucleus) can be used in one run
    velocity_group = parser_count.add_mutually_exclusive_group()
    velocity_group.add_argument(
        '--lamanno',
        help='Deprecated. Use `--workflow lamanno` instead.',
        action='store_true'
    )
    velocity_group.add_argument(
        '--nucleus',
        help='Deprecated. Use `--workflow nucleus` instead.',
        action='store_true'
    )
    
    # loom and h5ad files can't be generated simultaneously in one run
    conversion_group = parser_count.add_mutually_exclusive_group()
    conversion_group.add_argument(
        '--loom',
        help='Generate loom file from count matrix',
        action='store_true'
    )
    conversion_group.add_argument(
        '--h5ad',
        help='Generate h5ad file from count matrix',
        action='store_true'
    )
    return parser_count


def get_parser():
    """ Construct and return the ArgumentParser object that parses input command.
    """
    
    # Main parser
    parser = argparse.ArgumentParser(
        prog="run_kb_python.py",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='kb_python package'
    )
    parser._actions[0].help = parser._actions[0].help.capitalize()
    parser.add_argument(
        '--list',
        help='Display list of supported single-cell technologies',
        action='store_true'
    )
    parser.add_argument("--root_dir", default="/mnt")
    subparsers = parser.add_subparsers(
        dest='command',
        metavar='<CMD>',
    )

    # Add common options to this parent parser
    parent = argparse.ArgumentParser(add_help=False)
    parent.add_argument(
        '--workflow',
        help=(
            'Type of workflow. Use `lamanno` to calculate '
            'RNA velocity based on La Manno et al. 2018 logic. Use `nucleus` to '
            'calculate RNA velocity on single-nucleus RNA-seq reads (default: standard)'
        ),
        type=str,
        default='standard',
        choices=['standard', 'lamanno', 'nucleus', 'kite']
    )
    parent.add_argument(
        '--keep-tmp',
        help='Do not delete the tmp directory',
        action='store_true'
    )
    parent.add_argument(
        '--tmp_dir',
        help='AWS S3 bucket path to store temporary files generated in the run',
        type=str,
        required = '--keep-tmp' in sys.argv,
    )
    parent.add_argument(
        '--verbose', help='Print debugging information', action='store_true'
    )

    # Command parsers
    setup_info_args(subparsers, argparse.ArgumentParser(add_help=False))
    parser_ref = setup_ref_args(subparsers, parent)
    parser_count = setup_count_args(subparsers, parent)

    command_to_parser = {
        'ref': parser_ref,
        'count': parser_count,
    }   

    # Show help when no arguments are given
    if len(sys.argv) == 0:
        parser.print_help(sys.stderr)
        sys.exit(1)
    if len(sys.argv) == 1:
        if sys.argv[1] in command_to_parser:
            command_to_parser[sys.argv[1]].print_help(sys.stderr)
        else:
            parser.print_help(sys.stderr)
        sys.exit(1)
        
    return parser
        
    
def main(logger):
    """Command-line entrypoint. Download fastqs from the S3 bucket to an EC2 instance, build kallisto index, run kallisto peusoalignment, generate count matrices, and upload the data quantification results back to the S3 bucket.

    Keyword Argument:    
    logger -- Logger object that exposes the interface the code directly uses (type: logging.Logger)
    """
    # Parse input arguments
    parser = get_parser()
    args = parser.parse_args()
    print(f"testing purpose - what is sys.argv: {sys.argv}")
    print(f"testing purpose - what is args: {args}")
    
    # Root direcotry path on the EC2 instance
    root_dir = pathlib.Path(args.root_dir)

    if os.environ.get("AWS_BATCH_JOB_ID"):
        root_dir = root_dir / os.environ["AWS_BATCH_JOB_ID"]

    # Directory on the EC2 instance where relevant files for the run are stored
    run_dir = root_dir / "data"
    run_dir.mkdir(parents=True)
    
    # Run `kb --info` to see kb_python package details, or `kb --list` to see supported technologies list
    if 'info' in sys.argv:
        display_info(logger)
    elif '--list' in sys.argv:
        display_technologies(logger) 
    
    # Use the updated input format for non-standard workflows
    if any(arg in sys.argv for arg in {'--lamanno', '--nucleus'}):
        logger.warning((
            'The `--lamanno` and `--nucleus` flags are deprecated. '
            'These options will be removed in a future release. '
            'Please use `--workflow lamanno` or `--workflow nucleus` instead when the --workflow argument is available.'
        ))
    
    # Create tmp directory on the EC2 instance to store temporary files
    logger.debug('Creating tmp directory')
    tmp_dir = run_dir / "tmp"
    tmp_dir.mkdir(parents=True)
    
    # Run the command, and delete the temporary directory unless otherwise noted in the input
    try:
        logger.debug(args)
        print(f"testing purpose - see if args.command returns `ref` or `count`: {args.command}") # testing purpose
        COMMAND_TO_FUNCTION[args.command](args, run_dir, logger)
    finally:
        # Always clean temp dir
        if not args.keep_tmp:
            logger.debug('Removing tmp directory')
            shutil.rmtree(tmp_dir, ignore_errors=True)
            
    logger.info("Job completed")


if __name__ == "__main__":
    mainlogger, log_file, file_handler = ut_log.get_logger(__name__)
    s3c = boto3.client("s3")
    main(mainlogger)
