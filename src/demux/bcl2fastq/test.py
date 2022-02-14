import subprocess
## VIASH START

par = {
    "sample_sheet" : "resources_test/bcl/sample_sheet.csv",
    "input": "resources_test/bcl/input/",
    "output": "resources_test/bcl/output/",
    "reports": "resources_test/bcl/reports/",
    "skip_undetermined": True,
    "star_structure": True,
}

## VIASH END

# Run command: 

# get some data
input = 'resources_test/bcls'
sample_sheet = 'resources_test/sample_sheet/20220125_FS10000331_179_BRB11620-3029.csv'
output = 'resources_test/output'
executable = "bcl2fastq"

# Run memory monitor
command = (
    "while true;"
    ' do echo "memory usage" `cat /sys/fs/cgroup/memory/memory.usage_in_bytes`;'
    ' echo "disk usage" `df -h | grep "/mnt"`;'
    " sleep 300;"
    " done"
)
p = subprocess.Popen([command], shell=True)

# Run bcl2 fastq
command = [
    "bcl2fastq",
    "--sample-sheet", sample_sheet,
    "--runfolder-dir", input,
    "--output-dir", output,
]
pipe = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
if pipe.communicate()[1]:
    p.kill()
    raise RuntimeError("bcl2fastq failed: " + str(pipe[1]))

print("Completed Successfully!")