#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser()

parser.add_argument('taxon', choices=('mus', 'homo'))
parser.add_argument('num_partitions', type=int)
parser.add_argument('exp_ids', nargs='+')


parser.add_argument('--s3_input_dir', default='s3://czbiohub-seqbot/fastqs')
parser.add_argument('--queue', default='aegea_batch')
parser.add_argument('--vcpus', type=int, default=16)
parser.add_argument('--memory', type=int, default=64000)
parser.add_argument('--storage', type=int, default=500)

args = parser.parse_args()

for i in range(args.num_partitions):
    print(' '.join(('python aegea_launcher.py',
                    '--queue {}'.format(args.queue),
                    '--vcpus {}'.format(args.vcpus),
                    '--memory {}'.format(args.memory),
                    '--storage {}'.format(args.storage),
                    'jamestwebber-logs/new_scripts run_star_and_htseq.py',
                    '"--taxon {}'.format(args.taxon),
                    '--s3_input_dir {}'.format(args.s3_input_dir),
                    '--num_partitions {}'.format(args.num_partitions),
                    '--partition_id {}'.format(i),
                    '--exp_ids {}"'.format(' '.join(args.exp_ids))))
          )
    print('sleep 10')
