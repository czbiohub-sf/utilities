import os
import click
import logging
import subprocess

from ..demux.bcl2fastq import log_command


READ_NAMES = 'R1', 'R2'


def maybe_run_command(logger, command, retry_message, error_message, retry=5):
    for i in range(retry):
        try:
            log_command(logger, command, shell=True)
            break
        except subprocess.CalledProcessError:
            logger.info(retry_message)
    else:
        raise RuntimeError(error_message)


def assemble_transcripts(read1, read2, temp_folder='original_data'):
    """
    
    Parameters
    ----------
    read1, read2 : str
        Path(s) to fastq.gz files. Can be comma separated for multiple files
    """
    logging.basicConfig(level=logging.DEBUG)
    logger = logging.getLogger(__name__)


    read_pair = read1, read2

    for reads in read_pair:
        # Split to separate filenames if applicable
        if ',' in reads:
            reads = reads.split(',')
        else:
            reads = [reads]
        for read in reads:
            command = ['aws', 's3', 'cp', read, temp_folder + '/']
            maybe_run_command(logger, command, f'Retrying copying {read}',
                              f"Couldn't download {read}")

    # Concatenate gzip files
    data_folder = 'data'
    os.makedirs(data_folder)
    for read_name in READ_NAMES:
        command = ['zcat', f'{temp_folder}/*{read_name}*', '{data_folder}/{read_name}.fastq']

