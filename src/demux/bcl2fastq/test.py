
# get some data
input = meta["resources_dir"] + "/path/to/data"
output = "..."
executable = "./" + meta["functionality_name"]

# Run bcl2 fastq
command = [
    executable,
    "--sample-sheet", par["sample_sheet"],
    "--runfolder-dir", par["input"],
    "--output-dir", par["output"],
]
failed = log_command(
    logger, command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True
)

