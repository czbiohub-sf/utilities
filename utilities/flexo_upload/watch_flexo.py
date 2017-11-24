#!/usr/bin/env python
# script intended for cronjob to scan SEQS folders and upload new runs when
# they are finished.


import glob
import logging
import os
import subprocess
import sys
import time

from logging.handlers import TimedRotatingFileHandler

import boto3


ROOT_DIR = '/mnt/SEQS'
SEQS = ['MiSeq-01', 'NextSeq-01', 'NovaSeq-01']

# I believe these are the correct files for each sequencer
SEQ_FILES = {'MiSeq-01'  : 'RTAComplete.txt',
             'NextSeq-01': 'RunCompletionStatus.xml',
             'NovaSeq-01': 'SequenceComplete.txt'}

S3_BUCKET = 'czbiohub-seqbot'
S3_BCL_DIR = 'bcl'

# time to sleep between uploads
SLEEPY_TIME = 0.001 # 1/1000th of a second between every file...


def maybe_exit_process():
    # get all python pids
    pids = subprocess.check_output("pgrep python", shell=True).split()

    cmds = 0
    for pid in pids:
        with open(os.path.join('/proc', pid, 'cmdline')) as f:
            # get the cmdline and match against this script
            line = f.read().split('\x00')[0]
            if line == sys.argv[0]:
                cmds += 1

    # if there are more than one, exit this one
    if cmds > 1:
        sys.exit(0)
    else:
        # otherwise, take a short nap before starting
        time.sleep(5)


def scan_dir(seq_dir, client, logger):
    run_name = os.path.basename(seq_dir)
    logger.info("Getting file list from {}".format(
            os.path.join(S3_BUCKET, S3_BCL_DIR, run_name))
    )

    paginator = client.get_paginator('list_objects')
    response_iterator = paginator.paginate(
            Bucket=S3_BUCKET, Prefix=os.path.join(S3_BCL_DIR, run_name)
    )

    file_set = {r['Key'] for result in response_iterator
                for r in result.get('Contents', [])}
    logger.info("Found {} objects in {}".format(
            len(file_set), os.path.join(S3_BUCKET, S3_BCL_DIR, run_name))
    )

    return file_set


def main(logger, upload_set):
    logger.debug("Creating S3 client")
    client = boto3.client('s3')

    logger.info("Scanning {}...".format(ROOT_DIR))
    total_uploads = 0

    # for each sequencer, check for newly completed runs
    for seq in SEQS:
        logger.info(seq)
        fns = glob.glob(os.path.join(ROOT_DIR, seq, '[0-9]*', SEQ_FILES[seq]))
        logger.debug('{} ~complete runs in {}'.format(
                len(fns), os.path.join(ROOT_DIR, seq))
        )
        for fn in fns:
            seq_dir = os.path.dirname(fn)
            if seq_dir in upload_set:
                logger.debug('skipping {}, already uploaded'.format(seq_dir))
                continue

            if (seq != 'NovaSeq-01'
                or os.path.exists(os.path.join(seq_dir, 'CopyComplete.txt'))):
                    file_set = scan_dir(seq_dir, client, logger)
                    num_files = 0
                    uploads = 0
                    logger.info('syncing {}'.format(seq_dir))

                    seq_root = os.path.dirname(seq_dir)
                    for root, dirs, files in os.walk(seq_dir, topdown=True):
                        num_files += len(files)
                        logger.debug('syncing {} files in {}'.format(
                                len(files), root)
                        )
                        base_dir = root[(len(seq_root) + 1):]

                        for file_name in files:
                            s3_key = os.path.join(S3_BCL_DIR, base_dir,
                                                  file_name)
                            if s3_key not in file_set:
                                logger.debug('uploading {} to {}'.format(
                                        file_name, s3_key)
                                )
                                try:
                                    client.upload_file(
                                        Filename=os.path.join(root, file_name),
                                        Bucket=S3_BUCKET,
                                        Key=s3_key)
                                    uploads += 1
                                except IOError:
                                    logger.warning("couldn't read {}".format(
                                            file_name)
                                    )
                                time.sleep(SLEEPY_TIME)

                    if (uploads + len(file_set)) == num_files:
                        logger.info('{} is synced'.format(seq_dir))
                        upload_set.add(seq_dir)
                        logger.debug('added {} for upload_set'.format(seq_dir))

                    total_uploads += uploads

    logger.info("sync complete")
    logger.info("{} files uploaded".format(total_uploads))

    return upload_set


if __name__ == "__main__":
    # check for an existing process running
    maybe_exit_process()

    mainlogger = logging.getLogger(__name__)
    mainlogger.setLevel(logging.DEBUG)

    # create a logging format
    formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )

    info_log_file = '/home/utility/flexo_watcher.log'
    info_handler = TimedRotatingFileHandler(info_log_file, when='W0',
                                            backupCount=10)
    info_handler.setLevel(logging.INFO)

    debug_log_file = '/home/utility/flexo_debug.log'
    debug_handler = TimedRotatingFileHandler(debug_log_file, when='midnight',
                                             backupCount=3)
    debug_handler.setLevel(logging.DEBUG)

    info_handler.setFormatter(formatter)
    debug_handler.setFormatter(formatter)

    # add the handlers to the logger
    mainlogger.addHandler(info_handler)
    mainlogger.addHandler(debug_handler)

    upload_record_file = '/home/utility/flexo_record.txt'
    if os.path.exists(upload_record_file):
        mainlogger.debug('reading record file')
        with open(upload_record_file) as f:
            old_upload_set = {line.strip() for line in f}
    else:
        mainlogger.debug('no record file exists, syncing everything...')
        old_upload_set = set()

    mainlogger.info('{} runs recorded as uploaded'.format(len(old_upload_set)))

    updated_upload_set = main(mainlogger, old_upload_set.copy())

    mainlogger.info('synced {} new runs'.format(
            len(updated_upload_set) - len(old_upload_set))
    )
    with open(upload_record_file, 'w') as OUT:
        print >> OUT, '\n'.join(sorted(updated_upload_set))

    mainlogger.debug('wrote new record file')
