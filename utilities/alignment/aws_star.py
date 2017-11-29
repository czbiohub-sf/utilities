#!/usr/bin/env python3

import argparse

parser = argparse.ArgumentParser()

parser.add_argument('taxon', choices=('mus', 'homo'))
parser.add_argument('num_partitions', type=int)
parser.add_argument('exp_ids', nargs='+')

args = parser.parse_args()

for i in range(args.num_partitions):
    print(' '.join(('python3 aegea_launcher.py',
                    'alignment/run_star_and_htseq.py',
                    '"--taxon {}'.format(args.taxon),
                    '--num_partitions {}'.format(args.num_partitions),
                    '--partition_id {}'.format(i),
                    '--exp_ids {}"'.format(' '.join(args.exp_ids))))
          )
    print('sleep 10')
