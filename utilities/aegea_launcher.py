#!/usr/bin/env python

import argparse
import logging
import importlib.util
import json
import os
import subprocess

import boto3
import botocore.exceptions


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
            prog='aegea_launcher.py',
            description=(
                "Run any script as a batch job\n"
                "e.g. aegea_launcher.py my_bucket/my_scripts "
                "[script name] [script args...]"
            ),
            epilog="See https://github.com/czbiohub/utilities for more examples",
            add_help=False,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # basic usage
    basic_group = parser.add_argument_group('basic arguments')
    basic_group.add_argument(
            'script_name',
            help='Local path of the script to run, e.g. demux/bcl2fastq.py'
    )
    basic_group.add_argument(
            'script_args',
            help='Script arguments as a string. e.g. "--exp_id 170823_A001..".'
    )


    # instance requirements
    instance_group = parser.add_argument_group('customize the instance')
    image_group = instance_group.add_mutually_exclusive_group()
    image_group.add_argument('--ecr-image', metavar='ECR',
                             help='ECR image to use for the job')
    image_group.add_argument('--ami',
                             help='AMI to use for the job')

    instance_group.add_argument(
            '--queue', default='aegea_batch',
            help='Queue to submit the job'
    )
    instance_group.add_argument(
            '--vcpus', type=resource_range('vcpus', 1, 64),
            help='Number of vCPUs needed, e.g. 16')
    instance_group.add_argument(
            '--memory', type=resource_range('memory', 0, 256000),
            help='Amount of memory needed, in MB, e.g. 16000'
    )
    instance_group.add_argument(
            '--storage', type=resource_range('storage', 500, 16000),
            help='Request additional storage, in GiB (min 500)'
    )
    instance_group.add_argument(
            '--ulimits', metavar='U', default=None, nargs='+',
            help='Change instance ulimits, e.g. nofile:1000000'
    )
    instance_group.add_argument('--environment', metavar='ENV', default=None,
                                nargs='+', help='Set environment variables')

    # other arguments
    other_group = parser.add_argument_group('other options')
    other_group.add_argument('--dryrun', action='store_true',
                             help="Print the command but don't launch the job")
    other_group.add_argument('--s3_script_bucket', default='czbiohub-scripts',
                             help="S3 bucket containing scripts")
    other_group.add_argument('--s3_script_dir',
                             help="Path to scripts in S3 bucket, if any")
    other_group.add_argument('-u', '--upload', action='store_true',
                             help="Upload the script to S3 before running")
    other_group.add_argument('-d', '--debug', action='store_true',
                             help="Set logging to debug level")
    other_group.add_argument('-h', '--help', action='help',
                             help="show this help message and exit")

    args = parser.parse_args()

    if not os.path.exists(args.script_name):
        raise ValueError("Can't find script: {}".format(args.script_name))

    logger = get_logger(args.debug, args.dryrun)

    script_base = os.path.basename(args.script_name)


    logger.debug('Importing script as a module')
    module_name = os.path.splitext(script_base)[0]
    spec = importlib.util.spec_from_file_location(module_name,
                                                  args.script_name)
    script_module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(script_module)

    logger.debug('Checking for script default requirements')

    if hasattr(script_module, 'get_default_requirements'):
        script_reqs = script_module.get_default_requirements()
        logger.debug("{} defines default requirements: {}".format(
                args.script_name, script_reqs)
        )
        args = parser.parse_args(namespace=script_reqs)
    else:
        logger.warning(
                "{} does not define default requirements".format(
                        args.script_name
                )
        )

    if not (args.ecr_image or args.ami):
        args.ecr_image = 'sra_download'

    logger.debug('Testing script args')

    if hasattr(script_module, 'get_parser'):
        script_parser = script_module.get_parser()
        try:
            script_parser.parse_args(args.script_args.split())
        except:
            logger.error(
                    "{} failed with the given arg string\n\t{}".format(
                            args.script_name, args.script_args)
            )
            raise
    else:
        raise NotImplementedError(
                "{} must have a 'get_parser' method to test args".format(
                        args.script_name
                )
        )

    logger.debug('Script parsed args successfully')

    logger.debug("Starting S3 client")
    client = boto3.client('s3')

    s3_key = os.path.join(args.s3_script_dir or '', script_base)
    logger.debug("Filename: {}, Bucket: {}, Key: {}".format(
            args.script_name, args.s3_script_bucket, s3_key)
    )

    if args.upload:
        logger.info("Uploading {} to s3://{}".format(
                args.script_name, os.path.join(args.s3_script_bucket, s3_key))
        )
        if not args.dryrun:
            client.upload_file(
                    Filename=args.script_name,
                    Bucket=args.s3_script_bucket,
                    Key=s3_key
            )
    else:
        logger.info("Checking s3://{} for {}".format(
                os.path.join(args.s3_script_bucket, s3_key), args.script_name)
        )
        try:
            client.head_object(
                    Bucket=args.s3_script_bucket,
                    Key=s3_key
            )
        except botocore.exceptions.ClientError:
            raise ValueError("{} is not on S3, you should upload it.".format(
                    script_base))


    job_command = "aws s3 cp s3://{} .; chmod 755 {}; ./{} {}".format(
            os.path.join(args.s3_script_bucket, s3_key),
            script_base, script_base, args.script_args
    )

    aegea_command = ['aegea', 'batch', 'submit',
                     '--queue', args.queue,
                     '--vcpus', str(args.vcpus),
                     '--memory', str(args.memory)]

    if args.ecr_image:
        aegea_command.extend(['--ecr-image', args.ecr_image])
    elif args.ami:
        aegea_command.extend(['--ami', args.ami])
    else:
        raise ValueError('either --ecr-image or --ami is required')

    if args.storage:
        aegea_command.extend(['--storage', '/mnt={}'.format(args.storage)])

    if args.ulimits:
        aegea_command.extend(['--ulimits', ' '.join(args.ulimits)])

    if args.environment:
        aegea_command.extend(['--environment', ' '.join(args.environment)])

    aegea_command.extend(['--command', '"{}"'.format(job_command)])

    logger.info('executing command:\n\t{}'.format(' '.join(aegea_command)))
    if not args.dryrun:
        output = json.loads(subprocess.check_output(' '.join(aegea_command),
                                                    shell=True))
        logger.info('Launched job with jobId: {}'.format(output['jobId']))
