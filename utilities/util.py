import logging
import os
import subprocess
import threading

import multiprocessing as mp


def log_command(logger, command, **kwargs):
    logger.info(' '.join(command))
    output = subprocess.check_output(' '.join(command), **kwargs)
    logger.debug(output)


def log_command_to_queue(log_queue, command, **kwargs):
    log_queue.put((' '.join(command), logging.INFO))
    try:
        output = subprocess.check_output(' '.join(command), **kwargs)
        failed = False
    except subprocess.CalledProcessError:
        output = "Command failed with output {}".format(output)
        failed = True

    log_queue.put((output, logging.DEBUG))
    return failed


def process_logs(q, logger):
    for msg,level in iter(q.get, 'STOP'):
        if level == logging.INFO:
            logger.info(msg)
        else:
            logger.debug(msg)


def get_logger(name, debug=False, dryrun=False):
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)

    # create a logging format
    if dryrun:
        formatter = logging.Formatter(
                '%(asctime)s - %(name)s - %(levelname)s - (DRYRUN) - %(message)s'
        )
    else:
        formatter = logging.Formatter(
                '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        )

    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.DEBUG if debug else logging.INFO)
    stream_handler.setFormatter(formatter)

    logger.addHandler(stream_handler)

    if os.environ.get('AWS_BATCH_JOB_ID'):
        log_file = os.path.abspath(
                '{}.log'.format(os.environ['AWS_BATCH_JOB_ID'])
        )
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)

        # add the handlers to the logger
        logger.addHandler(file_handler)
    else:
        log_file = None
        file_handler = None

    return logger, log_file, file_handler


def get_thread_logger(logger):
    log_queue = mp.Queue()
    log_thread = threading.Thread(target=process_logs,
                                  args=(log_queue, logger))
    log_thread.start()

    return log_queue, log_thread
