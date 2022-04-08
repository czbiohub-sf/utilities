import subprocess

## VIASH START
meta = {
    "functionality_name": "check_samplesheet",
    "resources_dir": "resources_test"
}
## VIASH END

print("> Check whether the sample sheet of bs_195891710 passes successfully")
out = subprocess.run(
    [
        # "viash", "run", "src/qc/check_samplesheet/config.vsh.yaml", "--", 
        "./" + meta["functionality_name"],
        "--input", f"{meta['resources_dir']}/bs_195891710/bcl/SampleSheet.csv"
    ],
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE
)

assert out.returncode == 0, "Exit code should be 0"
assert "Formatting of all sample sheets are OK!" in out.stdout.decode(), "Stdout does not contain expected output"



print("> Check whether invalid csv fails")

lines = ['[Data]', 'one,two,three']
with open('samplesheet_fail.csv', 'w') as f:
    f.writelines(lines)

out = subprocess.run(
    [
        # "viash", "run", "src/qc/check_samplesheet/config.vsh.yaml", "--", 
        "./" + meta["functionality_name"],
        "--input", "samplesheet_fail.csv"
    ],
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE
)

assert out.returncode != 0, "Exit code should not be 0"
assert "Formatting of all sample sheets are OK!" not in out.stdout.decode(), "Stdout should not pretend everything is OK"


print("> Completed Successfully!")