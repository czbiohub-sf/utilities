import logging
import os
import subprocess

from logging.handlers import TimedRotatingFileHandler


def log_command(logger, command, **kwargs):
    logger.info(" ".join(command))

    proc = subprocess.run(" ".join(command), **kwargs)

    if proc.returncode != 0:
        logger.error("Command failed")
        if proc.stdout:
            logger.error(proc.stdout.decode())

        return True
    else:
        return False


def get_logger(name, debug=False, dryrun=False):
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)

    # create a logging format
    if dryrun:
        formatter = logging.Formatter(
            "%(asctime)s - %(name)s - %(levelname)s - (DRYRUN) - %(message)s"
        )
    else:
        formatter = logging.Formatter(
            "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
        )

    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.DEBUG if debug else logging.INFO)
    stream_handler.setFormatter(formatter)

    logger.addHandler(stream_handler)

    if os.environ.get("AWS_BATCH_JOB_ID"):
        log_file = os.path.abspath("{}.log".format(os.environ["AWS_BATCH_JOB_ID"]))
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)

        # add the handlers to the logger
        logger.addHandler(file_handler)
    else:
        log_file = None
        file_handler = None

    return logger, log_file, file_handler


def get_trfh_logger(name, *args):
    # function to create a rotating-file logger
    # with potentially multiple file handlers

    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)

    # create a logging format
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )

    for file_name, log_level, when, backup_count in args:
        log_handler = TimedRotatingFileHandler(
            file_name, when=when, backupCount=backup_count
        )
        log_handler.setLevel(log_level)
        log_handler.setFormatter(formatter)
        logger.addHandler(log_handler)

    return logger
