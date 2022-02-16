import subprocess
## VIASH START

par = {
    "sample_sheet" : "../resources_test/bcl/sample_sheet.csv",
    "input": "../resources_test/bcl/input/",
    "output": "../resources_test/bcl/output/",
    "reports": "../resources_test/bcl/reports/",
    "skip_undetermined": True,
    "star_structure": True,
}

## VIASH END

# Run command: viash test config.vsh.yaml 

# get some data
input = 'resources_test/bcls'
sample_sheet = 'resources_test/sample_sheet/*.csv'
output = 'resources_test/output'

# construct command args
print("Help")
command = [
    "bcl2fastq",
    "--runfolder-dir", par["input"],
    "--output-dir", par["output"],
]
if par["sample_sheet"] is not None:
    command = command + [ "--sample-sheet", par["sample_sheet"] ]

# run bcl2fastq
with subprocess.Popen(
    command,
    stdout=subprocess.PIPE,
    stderr=subprocess.STDOUT,
) as p:
    for line in p.stdout:
        print(line.decode(), end="")

if p.returncode > 0:
    raise RuntimeError(f"bcl2fastq failed with exit code {p.returncode}")

print("Completed Successfully!")