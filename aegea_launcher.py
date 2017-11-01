#!/usr/bin/env python

import argparse
import logging
import os
import subprocess

import boto3


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

    return resource_range


def get_logger(debug, dryrun):
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG if debug else logging.INFO)

    # create a logging format
    if dryrun:
        formatter = logging.Formatter(
                '%(asctime)s - %(name)s - %(levelname)s - (DRYRUN) - %(message)s'
        )
    else:
        formatter = logging.Formatter(
                '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        )

    handler = logging.StreamHandler()
    handler.setFormatter(formatter)

    logger.addHandler(handler)

    return logger


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description=(
                "Run any script as a batch job\n"
                "e.g. aegea_launcher.py my_bucket/my_scripts "
                "[script name] [script args...]"
            ),
            epilog="See https://github.com/czbiohub/utilities for more examples"
    )

    # basic usage
    basic_group = parser.add_argument_group('basic')
    basic_group.add_argument('s3_script_path',
                             help=('S3 bucket/path to scripts.\n'
                                   'e.g. jamestwebber-logs/scripts'))
    basic_group.add_argument('script_name',
                             help='Name of the script to run')
    basic_group.add_argument('script_args', nargs='*',
                             help='Script arguments. Wildcards not recommended')


    # instance requirements
    instance_group = parser.add_argument_group('instance')
    instance_group.add_argument('--ecr-image', default='sra_download',
                                help='ECR image to use for the job')
    # instance_group.add_argument('--ami', default=None,
    #                             help='AMI to use for the job')

    instance_group.add_argument('--queue', default='aegea_batch',
                           help='Queue to submit the job')
    instance_group.add_argument('--vcpus', default=1,
                           type=resource_range('vcpus', 1, 64),
                           help='Number of vCPUs needed, e.g. 16')
    instance_group.add_argument('--memory', default=4000,
                           type=resource_range('memory', 0, 256000),
                           help='Amount of memory needed, in MB, e.g. 16000')
    instance_group.add_argument('--storage', default=None,
                           type=resource_range('storage', 500, 16000),
                           help='Request additional storage, in GiB (min 500)')
    instance_group.add_argument('--ulimits', default=None, nargs='+',
                           help='Change instance ulimits, e.g. nofile:1000000')
    instance_group.add_argument('--environment', default=None, nargs='+',
                           help='Set environment variables')

    # other arguments
    other_group = parser.add_argument_group('other')
    other_group.add_argument('--dryrun', action='store_true',
                             help="Print the job command but don't launch it")
    other_group.add_argument('--upload', action='store_true',
                             help="Upload the script to S3 before running")
    other_group.add_argument('--debug', action='store_true',
                             help="Set logging to debug level")

    args = parser.parse_args()

    logger = get_logger(args.debug, args.dryrun)

    script_base = os.path.basename(args.script_name)

    s3_bucket = os.path.split(args.s3_script_path)[0]
    s3_key = os.path.join(os.sep.join(os.path.split(args.s3_script_path)[1:]),
                          script_base)

    if args.upload:
        if not os.path.exists(args.script_name):
            raise ValueError("Can't find script: {}".format(args.script_name))

        logger.debug("Starting S3 client")
        client = boto3.client('s3')

        logger.info("Uploading {} to s3://{}".format(
                args.script_name, os.path.join(s3_bucket, s3_key))
        )
        logger.debug("Filename: {}, Bucket: {}, Key: {}".format(
                args.script_name, s3_bucket, s3_key)
        )
        if not args.dryrun:
            client.upload_file(
                    Filename=args.script_name,
                    Bucket=s3_bucket,
                    Key=s3_key
            )

    job_command = "aws s3 cp s3://{} .; chmod 755 {}; ./{} {}".format(
            os.path.join(s3_bucket, s3_key),
            script_base, script_base, ' '.join(args.script_args)
    )

    aegea_command = ['aegea', 'batch', 'submit',
                     '--queue', args.queue,
                     '--vcpus', str(args.vcpus),
                     '--memory', str(args.memory),
                     '--ecr-image', args.ecr_image]

    if args.storage:
        aegea_command.extend(['--storage', '/mnt={}'.format(args.storage)])

    if args.ulimits:
        aegea_command.extend(['--ulimits', ' '.join(args.ulimits)])

    if args.environment:
        aegea_command.extend(['--environment', ' '.join(args.environment)])

    aegea_command.extend(['--command', '"{}"'.format(job_command)])

    logger.info('executing command:\n\t{}'.format(' '.join(aegea_command)))
    if not args.dryrun:
        output = subprocess.check_output(aegea_command, shell=True)
        logger.debug(output)
