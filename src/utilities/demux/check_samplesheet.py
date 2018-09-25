#!/usr/bin/env python

from __future__ import print_function

import argparse
import csv
import functools
import string


valid_chars = set(string.ascii_letters + string.digits + "_-")

parser = argparse.ArgumentParser()
parser.add_argument("samplesheet", nargs="+")

args = parser.parse_args()

for samplesheet in args.samplesheet:
    print("Checking {}".format(samplesheet))

    with open(samplesheet) as f:
        rdr = csv.reader(f)
        rows = list(rdr)

    if len(set(map(len, rows))) > 1:
        print("    Rows are not all the same length")

    if rows[0][0] != "[Data]":
        print("    Does not start with [Data] section header")

    all_char = functools.reduce(set.union, (v for r in rows[1:] for v in r), set())
    invalid_chars = all_char - valid_chars

    if invalid_chars:
        print("    Invalid characters in sample sheet: {}".format(invalid_chars))
