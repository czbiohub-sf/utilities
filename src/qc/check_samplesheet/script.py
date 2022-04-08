from __future__ import print_function

import csv
import functools
import string
import sys

## VIASH START
par = {
  'input': [ "resources_test/bs_195891710/bcl/SampleSheet.csv" ]
}
## VIASH END

valid_chars = set(string.ascii_letters + string.digits + "_-")

exit_code = 0

for samplesheet in par["input"]:
    print(f"Checking {samplesheet}")

    with open(samplesheet) as f:
        rdr = csv.reader(f)
        rows = list(rdr)

    if len(set(map(len, rows))) > 1:
        print("    Rows are not all the same length")
        exit_code = 1

    
    data_start = [ ix for ix, val in enumerate(rows) if val[0] == "[Data]" ]

    if len(data_start) != 1:
        print("    Sample sheet should have exactly one [Data] header")
        exit_code = 1

    data_start_ix = data_start[0]+1

    all_char = functools.reduce(
        set.union, 
        (val for row in rows[data_start_ix:] for val in row), 
        set()
    )
    invalid_chars = all_char - valid_chars

    if invalid_chars:
        print(f"    Invalid characters in sample sheet: {invalid_chars}")
        exit_code = 1

if exit_code == 0:
    print("Formatting of all sample sheets are OK!")
sys.exit(exit_code)