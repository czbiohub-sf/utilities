import subprocess

## VIASH START
meta = {
    "functionality_name": "bcl2fastq",
    "resources_dir": "resources_test"
}
## VIASH END

# get some data
bcl_data = f"{meta['resources_dir']}/sbs_195891710/bcl"
fastq_output = 'test_output'

# construct command args
command = [
    "./" + meta["functionality_name"],
    "--input", bcl_data,
    "--output", fastq_output,
]

print(f"> Running {meta['functionality_name']}")
out = subprocess.run(
    command,
    stdout=subprocess.PIPE,
    stderr=subprocess.STDOUT,
    check=True
)

print("Completed Successfully!")