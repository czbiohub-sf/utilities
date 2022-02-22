
import sys
import glob
import os
import re
import subprocess

## VIASH START
par = {
    "sample_sheet" : None,
    "input": "resources_test/bs_195891710/bcl_data/",
    "output": "resources_test/bs_195891710/fastqs/",
    "reports": "resources_test/bs_195891710/reports/",
    "skip_undetermined": True,
    "star_structure": True,
}
## VIASH END

if not os.path.isdir(par["output"]):
    print(f"creating output directory {par['output']}")
    os.makedirs(par["output"])

# construct command args
# hard-code path to bcl2fastq for now; change to
# 'bcl2fastq' when viash is updated to viash 0.5.8
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

# Fix directory structure of the files
fastqgz_files = glob.glob(os.path.join(par['output'], "*fastq.gz"))
print("found fastq.gz files:")
print("\n".join(fastqgz_files))

for fastq_file in fastqgz_files:
    fq_basename = os.path.basename(fastq_file)

    if par["skip_undetermined"] and fq_basename.startswith("Undetermined"):
        print(f"removing {fq_basename}")
        os.remove(fastq_file)
    elif par["star_structure"]:
        m = re.match("(.+)(_R[12]_001.fastq.gz)", fq_basename)
        if m:
            sample = m.group(1)
            sample_dir = os.path.join(par['output'], sample)
            if not os.path.exists(sample_dir):
                print(f"creating {sample_dir}")
                os.mkdir(sample_dir)
            print(f"moving {fastq_file}")
            os.rename(
                fastq_file,
                os.path.join(sample_dir, fq_basename),
            )
        else:
            print(f"Warning: regex didn't match {fastq_file}")

print("demux/bcl2fastq run completed.")
sys.stdout.flush()