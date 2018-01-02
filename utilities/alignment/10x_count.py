#!/usr/bin/env python

# Example: TAXON=homo CELL_COUNT=3000 S3_DIR=s3://biohub-spyros/data/10X_data/CK_Healthy/ ./10x_count.py
import argparse
import logging
import os
import sys
import subprocess
import tarfile


CELLRANGER = 'cellranger'
GC_TABLE_GENERATOR = 'generate_gc_table_from_cellranger.py'

S3_RETRY = 5
S3_LOG_DIR = 's3://jamestwebber-logs/10xcount_logs/'


def get_default_requirements():
    return argparse.Namespace(vcpus=64, memory=256000, storage=2000)


def get_parser():
    parser = argparse.ArgumentParser(
            prog='10x_count.py',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--exp_id', required=True)
    parser.add_argument('--taxon', default='homo', choices=('homo', 'mus'))

    parser.add_argument('--s3_input_dir',
                        default='s3://czbiohub-seqbot/fastqs')
    parser.add_argument('--s3_output_dir')
    parser.add_argument('--cell_count', type=int, default=3000)

    parser.add_argument('--make_gctable', action='store_true')
    parser.add_argument('--root_dir', default='/mnt')

    return parser


def log_command(logger, command, **kwargs):
    logger.info(' '.join(command))
    output = subprocess.check_output(' '.join(command), **kwargs)
    logger.debug(output)


def main(logger):
    parser = get_parser()

    args = parser.parse_args()

    if os.environ.get('AWS_BATCH_JOB_ID'):
        args.root_dir = os.path.join(args.root_dir,
                                     os.environ['AWS_BATCH_JOB_ID'])


    # local directories
    sample_id = os.path.basename(args.s3_input_dir)
    result_path = os.path.join(args.root_dir, 'data', 'hca', sample_id)
    fastq_path = os.path.join(result_path, 'fastqs')
    os.makedirs(fastq_path)

    genome_base_dir = os.path.join(args.root_dir, "genome", "cellranger")
    os.makedirs(genome_base_dir)

    # required environmental variable
    # taxon = os.environ.get('TAXON', 'homo')
    # s3_input_dir = os.environ['S3_INPUT_DIR'].rstrip('/')
    # s3_output_dir = os.environ['S3_OUTPUT_DIR'].rstrip('/')
    # cell_count = os.environ['CELL_COUNT']


    if args.taxon == 'homo':
        genome_name = 'HG38-PLUS'
    elif args.taxon == 'mus':
        genome_name = 'MM10-PLUS'
    else:
        raise ValueError("unknown taxon {}".format(args.taxon))




    # files that should be uploaded outside of the massive tgz
    # path should be relative to the run folder
    files_to_upload = [
        'outs/raw_gene_bc_matrices_h5.h5',
        'outs/raw_gene_bc_matrices/{}/genes.tsv'.format(genome_name),
        'outs/raw_gene_bc_matrices/{}/barcodes.tsv'.format(genome_name),
        'outs/raw_gene_bc_matrices/{}/matrix.mtx'.format(genome_name),
        'outs/web_summary.html',
        'outs/metrics_summary.csv'
    ]

    genome_tar_source = os.path.join('s3://czi-hca/ref-genome/cellranger/',
                                     genome_name + '.tgz')
    genome_dir = os.path.join(genome_base_dir, genome_name)

    # download the ref genome data
    command = ["aws", "s3", "cp", genome_tar_source, genome_base_dir]
    log_command(logger, command, shell=True)

    genome_tar_file = os.path.basename(genome_tar_source)
    logger.debug('Extracting {}'.format(genome_tar_file))
    with tarfile.open(os.path.join(genome_base_dir, genome_tar_file)) as tf:
        tf.extractall(path=genome_base_dir)


    sys.stdout.flush()

    # download the fastq files
    command = ['aws', 's3', 'cp', '--recursive', args.s3_input_dir, fastq_path]
    log_command(logger, command, shell=True)


    # Run cellranger
    os.chdir(result_path)
    command = [CELLRANGER, 'count', '--localmem=240',
               '--nosecondary', '--cells={}'.format(args.cell_count),
               '--sample={}'.format(sample_id), '--id={}'.format(sample_id),
               '--fastqs={}'.format(fastq_path),
               '--transcriptome={}'.format(genome_dir)]
    log_command(logger, command, shell=True,
                stderr=subprocess.STDOUT, universal_newlines=True)


    # Move results(websummary, cell-gene table, tarball) data back to S3
    for file_name in files_to_upload:
        command = ['aws', 's3', 'cp',
                   os.path.join(result_path, sample_id, file_name),
                   '{}/'.format(args.s3_output_dir)]
        for i in range(S3_RETRY):
            try:
                log_command(logger, command, shell=True)
                break
            except subprocess.CalledProcessError:
                logger.info("retrying cp {}".format(file_name))
        else:
            raise RuntimeError("couldn't sync {}".format(file_name))

    command = ['tar', 'cvzf',
               '{}.tgz'.format(os.path.join(result_path, sample_id)),
               sample_id]
    log_command(logger, command, shell=True)


    command = ['aws', 's3', 'cp',
               '{}.tgz'.format(os.path.join(result_path, sample_id)),
               '{}/'.format(args.s3_output_dir)]
    for i in range(S3_RETRY):
        try:
            log_command(logger, command, shell=True)
            break
        except subprocess.CalledProcessError:
            logger.info("retrying cp {}.tgz".format(sample_id))
    else:
        raise RuntimeError("couldn't sync {}.tgz".format(sample_id))


    if args.make_gctable:
        command = [GC_TABLE_GENERATOR, '-d',
                   os.path.join(sample_id, 'outs',
                                'raw_gene_bc_matrices', genome_name),
                   '-f', '{}.{}.cell-gene.csv'.format(sample_id, args.taxon),
                   '-m', '500']
        log_command(logger, command, shell=True)

        command = ['aws', 's3', 'cp',
                   '{}.{}.cell-gene.csv'.format(sample_id, args.taxon),
                   '{}/'.format(args.s3_output_dir)]

        for i in range(S3_RETRY):
            try:
                log_command(logger, command, shell=True)
                break
            except subprocess.CalledProcessError:
                logger.info("retrying cp {}.{}.cell-gene.csv".format(
                        sample_id, args.taxon))
        else:
            raise RuntimeError("couldn't sync {}.{}.cell-gene.csv".format(
                    sample_id, args.taxon))


if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)
    mainlogger = logging.getLogger(__name__)

    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)

    # create a logging format
    formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )

    stream_handler.setFormatter(formatter)

    mainlogger.addHandler(stream_handler)

    if os.environ.get('AWS_BATCH_JOB_ID'):
        log_file = '{}.log'.format(os.environ['AWS_BATCH_JOB_ID'])
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)

        # add the handlers to the logger
        mainlogger.addHandler(file_handler)
    else:
        log_file = None

    try:
        main(mainlogger)
    except:
        mainlogger.info("An exception occurred", exc_info=True)
        raise
    finally:
        if log_file:
            log_cmd = 'aws s3 cp {} {}'.format(log_file, S3_LOG_DIR)
            mainlogger.info(log_cmd)

            file_handler.close()
            subprocess.check_output(log_cmd, shell=True)
