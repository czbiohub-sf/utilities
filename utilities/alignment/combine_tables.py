#!/usr/bin/env python3


import argparse
import csv


def combine_files(fileA, fileB, output_file):
    with open(fileA) as fA, open(fileB) as fB:
        rdrA = csv.DictReader(fA)
        rdrB = csv.DictReader(fB)

        hA = rdrA.fieldnames
        hB = rdrB.fieldnames

        hU = sorted(set(hA[1:]) | set(hB[1:]))

        print('{} cells in file A, {} in file B, {} total'.format(
                len(hA), len(hB), len(hU))
        )

        with open(output_file, 'w') as OUT:
            print('gene,{}'.format(','.join(hU)), file=OUT)

            for i,(rA,rB) in enumerate(zip(rdrA, rdrB)):
                if rA['gene'] != rB['gene']:
                    raise ValueError(
                        'Gene list mismatch: {} != {}'.format(rA['gene'],
                                                              rB['gene'])
                    )

                summed_r = [rA['gene']]
                summed_r.extend(str(int(rA.get(c, 0)) + int(rB.get(c, 0)))
                                for c in hU)
                print(','.join(summed_r), file=OUT)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            prog='combine_tables.py',
            description=(
                "Combine the gene-cell counts from two flow cells\n"
                "e.g. ./combine_tables.py fileA fileB output_file\n"
                "Input files should have genes as rows and cells as columns."
            ),
    )

    parser.add_argument('fileA', help='Gene cell table for first flow cell')
    parser.add_argument('fileB', help='Gene cell table for second flow cell')
    parser.add_argument('output_file', help='File to write combined counts')

    args = parser.parse_args()

    combine_files(args.fileA, args.fileB, args.output_file)
