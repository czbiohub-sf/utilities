#!/usr/bin/env python3

import argparse
import importlib.util
import json
import os
import subprocess

import utilities.log_util as ut_log


REPO_ADDRESS = "https://github.com/czbiohub/utilities.git"


# helper function to check arguments are within a given range
def resource_range(name, min_val, max_val):
    def range_validator(s):
        value = int(s)
        if value < min_val:
            msg = "{} must be at least".format(name, min_val)
            raise argparse.ArgumentTypeError(msg)
        if value > max_val:
            msg = "{} can be at most".format(name, max_val)
            raise argparse.ArgumentTypeError(msg)
        return value

    return range_validator


def main():
    parser = argparse.ArgumentParser(
        prog="evros",
        description=(
            "Run batch jobs on AWS\n"
            "e.g. evros [options] demux.bcl2fastq [script args...]"
        ),
        epilog="See https://github.com/czbiohub/utilities for more examples",
        add_help=False,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # basic usage
    basic_group = parser.add_argument_group("basic arguments")
    basic_group.add_argument(
        "script_name", help="Local path of the script to run, e.g. demux.bcl2fastq"
    )
    basic_group.add_argument(
        "script_args",
        nargs=argparse.REMAINDER,
        help="Arguments for the script (everything after script_name)",
    )

    # instance requirements
    instance_group = parser.add_argument_group("customize the instance")
    image_group = instance_group.add_mutually_exclusive_group()
    image_group.add_argument(
        "--ecr-image",
        metavar="ECR",
        default="demuxer",
        help="ECR image to use for the job",
    )

    instance_group.add_argument(
        "--queue", default="aegea_batch", help="Queue to submit the job"
    )
    instance_group.add_argument(
        "--vcpus",
        type=resource_range("vcpus", 1, 64),
        help="Number of vCPUs needed, e.g. 16",
    )
    instance_group.add_argument(
        "--memory",
        type=resource_range("memory", 0, 256000),
        help="Amount of memory needed, in MB, e.g. 16000",
    )
    instance_group.add_argument(
        "--storage",
        type=resource_range("storage", 500, 16000),
        help="Request additional storage, in GiB (min 500)",
    )
    instance_group.add_argument(
        "--ulimits",
        metavar="U",
        default=None,
        nargs="+",
        help="Change instance ulimits, e.g. nofile:1000000",
    )
    instance_group.add_argument(
        "--environment",
        metavar="ENV",
        default=None,
        nargs="+",
        help="Set environment variables",
    )

    # other arguments
    other_group = parser.add_argument_group("other options")
    other_group.add_argument(
        "--dryrun",
        action="store_true",
        help="Print the command but don't launch the job",
    )
    other_group.add_argument(
        "--branch", default="master", help="branch of utilities repo to use"
    )
    other_group.add_argument(
        "-d", "--debug", action="store_true", help="Set logging to debug level"
    )
    other_group.add_argument(
        "-h", "--help", action="help", help="show this help message and exit"
    )

    args = parser.parse_args()

    logger = ut_log.get_logger(__name__, args.debug, args.dryrun)[0]

    logger.debug("Importing script as a module")
    if not args.script_name.startswith("."):
        args.script_name = ".{}".format(args.script_name)
    script_module = importlib.import_module(args.script_name, "utilities")
    script_base = os.path.basename(script_module.__file__)

    logger.debug("Checking for script default requirements")

    if hasattr(script_module, "get_default_requirements"):
        script_reqs = script_module.get_default_requirements()
        logger.debug(
            "{} defines default requirements: {}".format(args.script_name, script_reqs)
        )
        args = parser.parse_args(namespace=script_reqs)
    else:
        logger.warning(
            "{} does not define default requirements".format(args.script_name)
        )

    logger.debug("Testing script args")

    if hasattr(script_module, "get_parser"):
        script_parser = script_module.get_parser()
        try:
            script_parser.parse_args(args.script_args)
        except:
            logger.error(
                "{} failed with the given arg string\n\t{}".format(
                    args.script_name, args.script_args
                )
            )
            raise
    else:
        raise NotImplementedError(
            "{} must have a 'get_parser' method to test args".format(args.script_name)
        )

    logger.debug("Script parsed args successfully")

    job_command = "; ".join(
        (
            "PATH=$HOME/anaconda/bin:$PATH",
            "cd utilities",
            "git pull",
            "git checkout {}".format(args.branch),
            "python -m utilities.{} {}".format(
                args.script_name, " ".join(args.script_args)
            ),
        )
    )

    aegea_command = [
        "aegea",
        "batch",
        "submit",
        "--queue",
        args.queue,
        "--vcpus",
        str(args.vcpus),
        "--memory",
        str(args.memory),
        "--ecr-image",
        args.ecr_image,
    ]

    if args.storage:
        aegea_command.extend(["--storage", "/mnt={}".format(args.storage)])

    if args.ulimits:
        aegea_command.extend(["--ulimits", " ".join(args.ulimits)])

    if args.environment:
        aegea_command.extend(["--environment", " ".join(args.environment)])

    aegea_command.extend(["--command", "'{}'".format(job_command)])

    logger.info("executing command:\n\t{}".format(" ".join(aegea_command)))
    if not args.dryrun:
        output = subprocess.check_output(" ".join(aegea_command), shell=True)
        try:
            output = json.loads(output)["jobId"]
            logger.info("Launched job with jobId: {}".format(output))
        except json.decoder.JSONDecodeError:
            logger.info(output)
