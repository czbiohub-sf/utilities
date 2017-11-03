#!/usr/bin/env python
import argparse
import datetime
import os
import re
import subprocess
import tarfile
import thread


import threading
import multiprocessing as mp

import logging

S3_RETRY = 5
S3_LOG_DIR = 's3://jamestwebber-logs/star_logs/'

STAR = "/usr/local/bin/STAR"
HTSEQ = "htseq-count"
SAMTOOLS = "samtools"

COMMON_PARS = [STAR,
               '--outFilterType', 'BySJout',
               '--outFilterMultimapNmax', '20',
               '--alignSJoverhangMin', '8',
               '--alignSJDBoverhangMin', '1',
               '--outFilterMismatchNmax', '999',
               '--outFilterMismatchNoverLmax', '0.04',
               '--alignIntronMin', '20',
               '--alignIntronMax', '1000000',
               '--alignMatesGapMax', '1000000',
               '--outSAMstrandField', 'intronMotif',
               '--outSAMtype', 'BAM', 'Unsorted',
               '--outSAMattributes', 'NH', 'HI', 'NM', 'MD',
               '--genomeLoad', 'LoadAndKeep',
               '--outReadsUnmapped', 'Fastx',
               '--readFilesCommand', 'zcat']

CURR_MIN_VER = datetime.date(2017, 3, 1)


def process_logs(q, logger):
    for msg,level in iter(q.get, 'STOP'):
        if level == logging.INFO:
            logger.info(msg)
        else:
            logger.debug(msg)


def log_command(logger, command, **kwargs):
    logger.info(' '.join(command))
    output = subprocess.check_output(' '.join(command), **kwargs)
    logger.debug(output)


def log_command_to_queue(log_queue, command, **kwargs):
    log_queue.put(' '.join(command), logging.INFO)
    output = subprocess.check_output(' '.join(command), **kwargs)
    log_queue.put(output, logging.DEBUG)


# def run_sample(sample_name, exp_id):
def run_sample(star_queue, htseq_queue, log_queue,
               s3_input_dir, genome_dir, run_dir, n_proc):
    for sample_name,exp_id in iter(star_queue.get, 'STOP'):
        log_queue.put('{} - {}'.format(exp_id, sample_name), logging.INFO)
        dest_dir = os.path.join(run_dir, sample_name)
        os.makedirs(dest_dir)
        os.mkdir(os.path.join(dest_dir, 'rawdata'))
        os.mkdir(os.path.join(dest_dir, 'results'))

        # copy fastq.gz from s3 to local
        s3_source = os.path.join(s3_input_dir, exp_id, 'rawdata', sample_name)
        command = ['aws', 's3', 'cp', '--recursive',
                   s3_source, os.path.join(dest_dir, 'rawdata'),
                   '--exclude', "'*'",
                   '--include', "'*.fastq.gz'"]
        for i in range(S3_RETRY):
            try:
                log_command_to_queue(log_queue, command, shell=True)
                break
            except subprocess.CalledProcessError:
                log_queue.put("retrying data download - {}".format(i),
                              logging.DEBUG)
        else:
            raise RuntimeError("couldn't download fastq.gz files")

        # start running STAR
        # getting input files first

        reads = glob.glob(os.path.join(dest_dir, 'rawdata', '*.fastq.gz'))
        if not reads:
            log_queue.put("Empty reads for %s" % s3_source, logging.INFO)
            return

        os.makedirs(os.path.join(dest_dir, 'results', 'Pass1'))

        command = COMMON_PARS[:]
        command.extend(('--runThreadN', str(n_proc),
                        '--genomeDir', genome_dir,
                        '--readFilesIn', ' '.join(reads)))
        log_command_to_queue(log_queue, command, shell=True,
                             cwd=os.path.join(dest_dir, 'results', 'Pass1'))

        # running sam tools
        command = [SAMTOOLS, 'sort', '-m', '6000000000', '-o',
                   './Pass1/Aligned.out.sorted.bam', './Pass1/Aligned.out.bam']
        log_command_to_queue(log_queue, command, shell=True,
                    cwd=os.path.join(dest_dir, 'results'))

        # running samtools index -b
        command = [SAMTOOLS, 'index', '-b', 'Aligned.out.sorted.bam']
        log_command_to_queue(log_queue, command, shell=True,
                             cwd=os.path.join(dest_dir, 'results', 'Pass1'))


        # remove unsorted bam files
        os.remove(os.path.join(dest_dir, 'results', 'Pass1', 'Aligned.out.bam'))

        # remove fastq files
        for fastq_file in reads:
            os.remove(fastq_file)

        # generating files for htseq-count
        command = [SAMTOOLS, 'sort', '-m', '6000000000', '-n', '-o',
                   './Pass1/Aligned.out.sorted-byname.bam',
                   './Pass1/Aligned.out.sorted.bam']
        log_command_to_queue(log_queue, command, shell=True,
                             cwd=os.path.join(dest_dir, 'results'))

        # ready to be htseq-ed and cleaned up
        htseq_queue.put((exp_id, sample_name, dest_dir))


def run_htseq(htseq_queue, log_queue, s3_input_dir, taxon, sjdb_gtf):
    for exp_id, sample_name, dest_dir in iter(htseq_queue.get, 'STOP'):
        # running htseq

        command = [HTSEQ,
                   '-r', 'name', '-s', 'no', '-f', 'bam',
                   '-m', '-intersection-nonempty',
                   os.path.join(dest_dir, 'results', 'Pass1',
                                'Aligned.out.sorted-byname.bam'),
                   sjdb_gtf, '>', 'htseq-count.txt']
        log_command_to_queue(log_queue, command, shell=True,
                             cwd=os.path.join(dest_dir, 'results'))
        os.remove(os.path.join(dest_dir, 'results', 'Pass1',
                               'Aligned.out.sorted-byname.bam'))


        # compressed the results dir and move it to s3
        command = ['tar', '-cvfz',
                   '{}.{}.tgz',format(sample_name, taxon),
                   'results']
        log_command_to_queue(log_queue, command, shell=True, cwd=dest_dir)

        # copy htseq and log files out to s3
        s3_dest = os.path.join(s3_input_dir, exp_id, 'results/')

        src_files = [
            os.path.join(dest_dir, '{}.{}.tgz'.format(sample_name, taxon)),
            os.path.join(dest_dir, 'results', 'htseq-count.txt'),
            os.path.join(dest_dir, 'results', 'Pass1', 'Log.final.out')
        ]

        dest_names = [
            '{}.{}.tgz'.format(sample_name, taxon),
            '{}.{}.htseq-count.txt'.format(sample_name, taxon),
            '{}.{}.log.final.out'.format(sample_name, taxon)
        ]

        for src_file,dest_name in zip(src_files, dest_names):
            command = ['aws', 's3', 'cp', file_name,
                       os.path.join(s3_dest, dest_name)]

            log_command_to_queue(log_queue, command, shell=True)

        # rm all the files
        command = ['rm', '-rf', dest_dir]
        log_command_to_queue(log_queue, command, shell=True)


def main(logger):
    parser = argparse.ArgumentParser()

    parser.add_argument('--root_dir', default='/mnt')
    parser.add_argument('--taxon', default='homo', choices=('homo', 'mus'))

    parser.add_argument('--s3_input_dir')
    parser.add_argument('--num_partitions', type=int)
    parser.add_argument('--partition_id', type=int)
    parser.add_argument('--exp_ids', nargs='+')

    parser.add_argument('--star_proc', type=int, default=8,
                        help='Number of processes to give to each STAR run')
    parser.add_argument('--htseq_proc', type=int, default=4,
                        help='Number of htseq processes to run')

    args = parser.parse_args()

    if os.environ.get('AWS_BATCH_JOB_ID'):
        args.root_dir = os.path.join(args.root_dir, os.environ['AWS_BATCH_JOB_ID'])

    run_dir = os.path.join(args.root_dir, 'data' 'hca')

    if args.taxon == 'homo':
        genome_dir = os.path.join(args.root_dir, "genome/STAR/HG38-PLUS/") # change
        ref_genome_file = 'hg38-plus.tgz'
        ref_genome_star_file = 'HG38-PLUS.tgz'
        sjdb_gtf = os.path.join(args.root_dir, 'genome', 'hg38-plus',
                                'hg38-plus.gtf')
    elif args.taxon == 'mus':
        genome_dir = os.path.join(args.root_dir, "genome/STAR/MM10-PLUS/")
        ref_genome_file = 'mm10-plus.tgz'
        ref_genome_star_file = 'MM10-PLUS.tgz'
        sjdb_gtf = os.path.join(args.root_dir, 'genome', 'mm10-plus',
                                'mm10-plus.gtf')

    else:
        raise ValueError('Invalid taxon {}'.format(args.taxon))

    if args.star_proc > mp.cpu_count():
        raise ValueError('Not enough CPUs to give {} processes to STAR'.format(
                args.star_proc))


    logger.info(
            '''Run Info: partition {} out of {}
                   star_proc:\t{}
                  htseq_proc:\t{}
                  genome_dir:\t{}
             ref_genome_file:\t{}
        ref_genome_star_file:\t{}
                    sjdb_gtf:\t{}
                       taxon:\t{}
                s3_input_dir:\t{}
                     exp_ids:\t{}'''.format(
                    args.partition_id, args.num_partitions,
                    args.star_proc, args.htseq_proc,
                    genome_dir, ref_genome_file,
                    ref_genome_star_file, sjdb_gtf,
                    args.taxon, args.s3_input_dir,
                    ', '.join(args.exp_ids)
            )
    )

    # download the genome data
    os.makedirs(os.path.join(args.root_dir, 'genome'))
    command = ['aws', 's3', 'cp',
               os.path.join('s3://czi-hca', 'ref-genome', ref_genome_file),
               os.path.join(args.root_dir, 'genome/')]
    log_command(logger, command, shell=True)

    logger.debug('Extracting {}'. format(ref_genome_file))
    with tarfile.open(os.path.join(args.root_dir, 'genome',
                                   ref_genome_file)) as tf:
        tf.extractall(path=os.path.join(args.root_dir, 'genome'))


    # download STAR stuff
    os.makedirs(os.path.join(args.root_dir, 'genome', 'STAR'))
    command = ['aws', 's3', 'cp',
               os.path.join('s3://czi-hca', 'ref-genome', 'STAR',
                            ref_genome_star_file),
               os.path.join(args.root_dir, 'genome', 'STAR/')]
    log_command(logger, command, shell=True)

    logger.debug('Extracting {}'.format(ref_genome_star_file))
    with tarfile.open(os.path.join(args.root_dir, 'genome',
                                   'STAR', ref_genome_star_file)) as tf:
        tf.extractall(path=os.path.join(args.root_dir ,'genome', 'STAR'))


    # Load Genome Into Memory
    command = [STAR, '--genomeDir', genome_dir, '--genomeLoad', 'LoadAndExit']
    log_command(logger, command, shell=True)

    log_queue = mp.Queue()
    log_thread = threading.Thread(target=process_logs,
                                  args=(log_queue, logger))
    log_thread.start()

    star_queue = mp.Queue()
    htseq_queue = mp.Queue()

    n_star_procs = mp.cpu_count() / args.star_proc

    star_args = (star_queue, htseq_queue, log_queue, args.s3_input_dir,
                 genome_dir, run_dir, args.star_proc)
    star_procs = [mp.Process(target=run_sample, args=star_args)
                  for i in range(n_star_procs)]

    for p in star_procs:
        p.start()

    htseq_args = (htseq_queue, log_queue, args.s3_input_dir,
                  args.taxon, sjdb_gtf)
    htseq_procs = [mp.Process(target=run_htseq, args=htseq_args)
                   for i in range(args.htseq_proc)]

    for p in htseq_procs:
        p.start()


    for exp_id in args.exp_ids:
        if exp_id.startswith('Undetermined'):
            logger.info("Skipping file: %s" % exp_id)
            continue

        # Check the exp_id folder for existing runs
        try:
            command = ['aws', 's3', 'ls',
                       os.path.join(args.s3_input_dir, exp_id, 'results/')]
            logger.info(' '.join(command))
            output = subprocess.check_output(' '.join(command),
                                             shell=True).split('\n')
        except subprocess.CalledProcessError:
            logger.info("Nothing in the results directory")
            output = []

        output_files = {(tuple(line[:10].split('-')), line.split()[-1])
                        for line in output
                        if line.strip().endswith('htseq-count.txt')}
        output_files = {fn for dt,fn in output_files
                        if datetime.date(*map(int, dt)) > CURR_MIN_VER}
        logger.info("number of files: {}".format(len(output_files)))

        logger.info("Running partition {} of {} for exp {}".format(
                args.partition_id, args.num_partitions, exp_id)
        )

        command = ['aws', 's3', 'ls',
                   os.path.join(args.s3_input_dir, exp_id, 'rawdata/')]
        logger.info(' '.join(command))
        try:
            output = subprocess.check_output(' '.join(command),
                                             shell=True).split("\n")
        except subprocess.CalledProcessError:
            logger.info("Nothing in the rawdata directory", exc_info=True)
            output = []

        sample_list = []

        for f in output:
            matched = re.search("\s([\d\w\-.]+)/", f)
            if matched:
                sample_list.append(matched.group(1))

        for sample_name in sample_list[args.partition_id::args.num_partitions]:
            if ('{}.{}.htseq-count.txt'.format(sample_name, args.taxon)
                in output_files):
                logger.info("{} already exists, skipping".format(sample_name))
                continue

            logger.info("Adding sample {} to queue".format(sample_name))
            star_queue.put((sample_name, exp_id))

    for i in range(n_star_procs):
        star_queue.put('STOP')

    for p in star_procs:
        p.join()

    for i in range(args.htseq_proc):
        htseq_queue.put('STOP')

    for p in htseq_procs:
        p.join()

    log_queue.put('STOP')
    log_thread.join()

    # Remove Genome from Memory
    command = [STAR, '--genomeDir', genome_dir, '--genomeLoad', 'Remove']
    log_command(logger, command, shell=True)

    logger.info('Job completed')


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
