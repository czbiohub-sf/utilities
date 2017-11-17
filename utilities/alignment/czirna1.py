#!/usr/bin/env python

import os
import click
import logging
import subprocess

# from ..demux.bcl2fastq import log_command


READ_NAMES = 'R1', 'R2'

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

S3_LOG_DIR = 's3://olgabot-transcript-assembly/logs/'

stream_handler = logging.StreamHandler()
stream_handler.setLevel(logging.INFO)

# create a logging format
formatter = logging.Formatter(
    '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)

stream_handler.setFormatter(formatter)

logger.addHandler(stream_handler)

# Set up where to write the logs to
if os.environ.get('AWS_BATCH_JOB_ID'):
    log_file = '{}.log'.format(os.environ['AWS_BATCH_JOB_ID'])
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)

    # add the handlers to the logger
    logger.addHandler(file_handler)
else:
    log_file = None


def log_command(logger, command, **kwargs):
    logger.info(' '.join(command))
    output = subprocess.check_output(' '.join(command), **kwargs)
    logger.debug(output)


def maybe_run_command(logger, command, retry_message, error_message, retry=5):
    for i in range(retry):
        try:
            log_command(logger, command, shell=True)
            break
        except subprocess.CalledProcessError:
            logger.info(retry_message)
    else:
        raise RuntimeError(error_message)


@click.command()
@click.argument('read1')
@click.argument('read2')
@click.option('--output-bucket', default='s3://olgabot-transcript-assembly')
@click.option('--temp-folder', default='original_data')
def cli(read1, read2, output_bucket='s3://olgabot-transcript-assembly',
        temp_folder='original_data'):
    """
    
    Parameters
    ----------
    read1, read2 : str
        Path(s) to fastq.gz files. Can be comma separated for multiple files
    temp_folder : str
        Location to store temporary data
    """
    temp_folder = temp_folder if temp_folder.endswith('/') else temp_folder + '/'

    command = ('while true;'
               ' do echo "memory usage" `cat /sys/fs/cgroup/memory/memory.usage_in_bytes`;'
               ' echo "disk usage" `df -h | grep "/mnt"`;'
               ' sleep 90;'
               ' done')
    memory_checker = subprocess.Popen([command], shell=True)
    read_pair = read1, read2
    basename = None

    for reads in read_pair:

        # Split to separate filenames if applicable
        if ',' in reads:
            reads = reads.split(',')
        else:
            reads = [reads]
        for read in reads:
            # Get a name for the sample using the first read
            if basename is None:
                basename = os.path.basename(read)

            command = ['aws', 's3', 'cp', read, temp_folder]
            maybe_run_command(logger, command, f'Retrying copying {read}',
                              f"Couldn't download {read}")

    # Concatenate gzip files
    data_folder = 'data'
    os.makedirs(data_folder)
    for read_name in READ_NAMES:
        command = ['zcat', f'{temp_folder}/*{read_name}*',
                   f'{data_folder}/{read_name}.fastq']
        maybe_run_command(logger, command,
                          f"Retrying concatenating {read_name} fastq files",
                          f"Couldn't concatenate {read_name} fastq files")

    # Run steps 1, 2, 7 of the command
    flags = map(str, (1, 2, 7))
    for flag in flags:
        command = ['time', 'Spyros-3cells.py', flag]
        maybe_run_command(logger, command,
                          f"Retrying Paolo's code with flag: {flag}",
                          f"Couldn't run Paolo's code with flag: {flag}")

    # Copy assembled transcripts and data back to s3
    output = f'{output_bucket}/{basename}'
    command = ['aws', 's3', 'cp', '--exclude', '*.fastq', 'data/', output]
    maybe_run_command(logger, command,
                      f"Retrying copying outputs to {output}",
                      f"Couldn't copy outputs to {output}")
    memory_checker.kill()


if __name__ == "__main__":

    try:
        cli()
    except:
        # Save to file
        if log_file:
            log_cmd = 'aws s3 cp {} {}'.format(log_file, S3_LOG_DIR)
            logger.info(log_cmd)

            file_handler.close()
            subprocess.check_output(log_cmd, shell=True)
