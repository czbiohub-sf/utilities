#!/usr/bin/env python

# template for writing new AWS batch scripts

import argparse
import subprocess

from utilities.log_util import get_logger, log_command


# an s3 bucket to upload your logs, if you want them
S3_LOG_DIR = "s3://your-bucket/script_logs/"


def get_parser():
    """
    Required: Define the arguments your script takes in this function,
    so that evros can test them before launching a job
    """
    parser = argparse.ArgumentParser(
        prog="template.py", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument("--message", default="Hello world!")

    return parser


def get_default_requirements():
    """
    Optional: Define the default hardware requirements for this job
    """
    return argparse.Namespace(vcpus=4, memory=8000, storage=500)


def main(logger):
    # get the argument parser and parse args
    parser = get_parser()
    args = parser.parse_args()

    # use the logger
    logger.info("Attempting to echo the message...")

    # run a subprocess and log the attempt
    failed = log_command(logger, "echo {}".format(args.message), shell=True)


if __name__ == "__main__":
    mainlogger, log_file, file_handler = get_logger(__name__)

    try:
        main(mainlogger)
    except:
        mainlogger.info("An exception occurred", exc_info=True)
        raise
    finally:
        # upload the log file no matter what. You can remove this if you don't
        # want to accumulate logs
        if log_file:
            log_cmd = "aws s3 cp --quiet {} {}".format(log_file, S3_LOG_DIR)
            mainlogger.info(log_cmd)

            file_handler.close()
            subprocess.check_output(log_cmd, shell=True)
