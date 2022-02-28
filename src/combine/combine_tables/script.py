
#!/usr/bin/env python3

import csv

## VIASH START
par = {
  'input': [ "path/to/file1", "path/to/file2" ],
  'output': '/path/to/file'
}
## VIASH END

with open(par["input"][0]) as file1, open(par["input"][1]) as file2:
    reader_file1 = csv.DictReader(file1)
    reader_file2 = csv.DictReader(file2)

    header_file1 = reader_file1.fieldnames
    header_file2 = reader_file2.fieldnames

    header_union = sorted(set(header_file1[1:]) | set(header_file2[1:]))

    print(
        "{} cells in file 1, {} in file 1, {} total".format(
            len(header_file1), len(header_file2), len(header_union)
        )
    )

    with open(par["output"], "w") as OUT:
        print("gene,{}".format(",".join(header_union)), file=OUT)

        for i, (value_file1, value_file2) in enumerate(zip(reader_file1, reader_file2)):
            if value_file1["gene"] != value_file2["gene"]:
                raise ValueError(
                    "Gene list mismatch: {} != {}".format(value_file1["gene"], value_file2["gene"])
                )

            summed_r = [value_file1["gene"]]
            summed_r.extend(str(int(value_file1.get(c, 0)) + int(value_file2.get(c, 0))) for c in header_union)
            print(",".join(summed_r), file=OUT)
