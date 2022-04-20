import subprocess
import os
import shutil

## VIASH START
par = {
    "input_fastqs": "resources_test/demuxed/fastqs/",
    "input_tif": "resources_test/demuxed/*.tif",
    "output": "aligned",
    "transcriptome": "resources/refdata-gex-mm10-2020-A.tar.gz",
    "id": "output",
    "sample": None,
    "lanes": None,
    "slide": None
}
## VIASH END

command = [
      "spaceranger", "count",
      "--fastqs", par["input_fastqs"],
      "--image", par["input_tif"],
      "--id", par["id"],
      "--transcriptome", par["reference_genome"],
]

if par["sample"] is not None:
      command = command + ["--sample", par["sample"]]

if par["lanes"] is not None:
      command = command + ["--lanes", par["lanes"]]
if par["slide"] is not None:
    command = command + ["--slide", par["slide"]]
else:
    command = command + ["--unknown-slide"]

# run bcl2fastq
with subprocess.Popen(
    command,
    stdout=subprocess.PIPE,
    stderr=subprocess.STDOUT,
) as p:
    for line in p.stdout:
        print(line.decode(), end="", flush=True)

if p.returncode > 0:
    raise RuntimeError(f"spaceranger count failed with exit code {p.returncode}")

if os.path.exists(f"${par['id']}/outs/"):
    os.makedirs(par["output"])
    shutil.move(f"${par['id']}/outs/", par["output"])