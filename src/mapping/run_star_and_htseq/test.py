import subprocess

## VIASH START
meta = {
    "functionality_name": "run_star_and_htseq",
    "resources_dir": "resources_test"
}
## VIASH END

# get some data
fastq_data = f"{meta['resources_dir']}/path/to/fastqs"
aligned_output = 'test_output'

# construct command args
command = [
    "./" + meta["functionality_name"],
    "--input", fastq_data,
    "--output", aligned_output,
]

print(f"> Running {meta['functionality_name']}")
out = subprocess.run(
    command,
    stdout=subprocess.PIPE,
    stderr=subprocess.STDOUT,
    check=True
)

print("Completed Successfully!")