#!/usr/bin/env python

import argparse
import logging
import os
import subprocess

CELLRANGER = 'cellranger'

S3_RETRY = 5
S3_LOG_DIR = 's3://jamestwebber-logs/mkfastq_logs/'


def get_default_requirements():
    return argparse.Namespace(vcpus=64, memory=256000, storage=2000,
                              queue='aegea_batch_demux',
                              ulimits=['nofile:1000000'])


def get_parser():
    parser = argparse.ArgumentParser(
            prog='10x_mkfastq.py',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--exp_id', required=True)

    parser.add_argument('--s3_input_dir',
                        default='s3://czbiohub-seqbot/bcl')
    parser.add_argument('--s3_output_dir',
                        default='s3://czbiohub-seqbot/fastqs')
    parser.add_argument('--s3_report_dir',
                        default='s3://czbiohub-seqbot/reports')
    parser.add_argument('--s3_sample_sheet_dir',
                        default='s3://czbiohub-seqbot/sample-sheets')

    parser.add_argument('--sample_sheet_name', default=None,
                        help='Defaults to [exp_id].csv')
    parser.add_argument('--root_dir', default='/mnt')

    return parser


def log_command(logger, command, **kwargs):
    logger.info(' '.join(command))
    output = subprocess.check_output(' '.join(command), **kwargs)
    logger.debug(output)


def main(logger):
    parser = get_parser()

    args = parser.parse_args()

    if os.environ.get('AWS_BATCH_JOB_ID'):
        args.root_dir = os.path.join(args.root_dir,
                                     os.environ['AWS_BATCH_JOB_ID'])


    if args.sample_sheet_name is None:
        args.sample_sheet_name = '{}.csv'.format(args.exp_id)

    # local directories
    result_path = os.path.join(args.root_dir, 'data', 'hca', args.exp_id)
    bcl_path = os.path.join(result_path, 'bcl')
    output_path = os.path.join(result_path, 'fastqs')

    os.makedirs(result_path)
    os.mkdir(bcl_path)

    # download sample sheet
    command = ['aws', 's3', 'cp', '--quiet',
               os.path.join(args.s3_sample_sheet_dir, args.sample_sheet_name),
               result_path]
    for i in range(S3_RETRY):
        try:
            log_command(logger, command, shell=True)
            break
        except subprocess.CalledProcessError:
            logger.info("retrying s3 copy")
    else:
        raise RuntimeError("couldn't download sample sheet {}".format(
                os.path.join(args.s3_sample_sheet_dir, args.sample_sheet_name))
        )


    # download the bcl files
    command = ['aws', 's3', 'sync', '--quiet',
               os.path.join(args.s3_input_dir, args.exp_id), bcl_path]
    for i in range(S3_RETRY):
        try:
            log_command(logger, command, shell=True)
            break
        except subprocess.CalledProcessError:
            logger.info("retrying s3 sync bcl")
    else:
        raise RuntimeError("couldn't sync {}".format(
                os.path.join(args.s3_input_dir, args.exp_id))
        )


    # Run cellranger mkfastq
    command = [CELLRANGER, 'mkfastq', '--localmem=60',
               '--sample-sheet={}'.format(os.path.join(result_path, args.sample_sheet_name)),
               '--run={}'.format(os.path.join(bcl_path)),
               '--output-dir={}'.format(output_path)]
    log_command(logger, command, shell=True)


    # upload fastq files to destination folder
    command = ['aws', 's3', 'sync', '--quiet',
               output_path, args.s3_output_dir]
    for i in range(S3_RETRY):
        try:
            log_command(logger, command, shell=True)
            break
        except subprocess.CalledProcessError:
            logger.info("retrying sync fastq")
    else:
        raise RuntimeError("couldn't sync fastqs")


if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)
    mainlogger = logging.getLogger(__name__)

    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)

    # create a logging format
    formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )

    stream_handler.setFormatter(formatter)

    mainlogger.addHandler(stream_handler)

    if os.environ.get('AWS_BATCH_JOB_ID'):
        log_file = '{}.log'.format(os.environ['AWS_BATCH_JOB_ID'])
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)

        # add the handlers to the logger
        mainlogger.addHandler(file_handler)
    else:
        log_file = None

    try:
        main(mainlogger)
    except:
        mainlogger.info("An exception occurred", exc_info=True)
        raise
    finally:
        if log_file:
            log_cmd = 'aws s3 cp --quiet {} {}'.format(log_file, S3_LOG_DIR)
            mainlogger.info(log_cmd)

            file_handler.close()
            subprocess.check_output(log_cmd, shell=True)
