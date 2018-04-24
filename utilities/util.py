import logging
import os
import subprocess


def log_command(logger, command, **kwargs):
    logger.info(' '.join(command))
    output = subprocess.check_output(' '.join(command), **kwargs)
    logger.debug(output)


def log_command_to_queue(log_queue, command, **kwargs):
    log_queue.put((' '.join(command), logging.INFO))
    output = subprocess.check_output(' '.join(command), **kwargs)
    log_queue.put((output, logging.DEBUG))


def process_logs(q, logger):
    for msg,level in iter(q.get, 'STOP'):
        if level == logging.INFO:
            logger.info(msg)
        else:
            logger.debug(msg)


def get_logger(name):
    logging.basicConfig(level=logging.DEBUG)
    logger = logging.getLogger(name)

    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)

    # create a logging format
    formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )

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


def maybe_make_directory(folder):
    try:
        os.makedirs(folder)
    except OSError:
        pass
