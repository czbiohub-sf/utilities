# utilities
A collection of scripts for common data management and processing tasks

## Quick Reference

| Task | Command | Description |
| ----------- | -------- | ----------- |
| Demux a standard sequencing run | `evros demux.bcl2fastq --exp_id YYMMDD_EXP_ID` | Assumes your sample sheet is uploaded to S3. If planning to run the alignment script, use `--star_structure` |
| Demux a 10X run | `evros demux.10x_mkfastq --exp_id YYMMDD_EXP_ID` | Again, assumes a sample sheet is present on S3 | 
| Align using STAR and htseq | `aws_star [mus or homo] [# partitions] YYMMDD_EXP_ID > your_script.sh` | Creates a shell script locally to launch many alignments using `source your_script.sh` |
| Align a 10X run | `evros alignment.10x_count --taxon [mus or homo] --s3_input_dir s3://czbiohub-seqbot/fastqs/YYMMDD_EXP_ID/SAMPLE --s3_output_dir s3://output-bucket/` | Run once for each channel of the run. Very slow! |
| Create a download token | `aws_access fastqs/YYMMDD_EXP_ID > download_instructions.txt` | Currently only works for paths within `s3://czbiohub-seqbot` |



## Installation

Clone this report and make a new environment, then install the package:

Create a conda environment called `utilities-env` specific to this Python package:

**Note: this package is written for Python 3, if you are using anaconda2 this will not work (and you should update)**

```
# Create an environment called `utilities-env` with Python and pip
➜  conda create --name utilities-env python pip
➜  source activate utilities-env
(utilities-env) ➜
```
The `(utilities-env)` indicates that the `utilities-env` environment is active.

Now let's clone the environment into a `code` folder:

```
(utilities-env) ➜ cd code
(utilities-env) ➜ git clone https://github.com/czbiohub/utilities.git

```

Now we'll change to that folder with `cd` and install the program with `pip install -e .`, where the `.` means right here.

```
(utilities-env) ➜ cd utilities
(utilities-env) ➜ pip install -e .
```

## Usage

As of the latest version of this repository, you *do not* need to be in the `utilities` folder to run scripts, but you *must* activate the environment:

```
➜  source activate utilities-env
(utilities-env) ➜
```

The `evros` command (and others) will be on your path, and it will be able to find the scripts to run on AWS.

### How to run a script on AWS, in general:

The `evros` script wraps aegea (which can be tricky to use). It has sane defaults so you don't need to choose most of them.

```
(utilities-env) ➜ evros --help
usage: evros [--ecr-image ECR] [--queue QUEUE] [--vcpus VCPUS]
             [--memory MEMORY] [--storage STORAGE] [--ulimits U [U ...]]
             [--environment ENV [ENV ...]] [--dryrun]
             [--s3_script_bucket S3_SCRIPT_BUCKET]
             [--s3_script_dir S3_SCRIPT_DIR] [-u] [-d] [-h]
             script_name ...

Run batch jobs on AWS e.g. evros [options] demux.bcl2fastq [script args...]

basic arguments:
  script_name           Local path of the script to run, e.g. demux.bcl2fastq
  script_args           Arguments for the script (everything after
                        script_name)

customize the instance:
  --ecr-image ECR       ECR image to use for the job (default: sra_download)
  --queue QUEUE         Queue to submit the job (default: aegea_batch)
  --vcpus VCPUS         Number of vCPUs needed, e.g. 16 (default: None)
  --memory MEMORY       Amount of memory needed, in MB, e.g. 16000 (default:
                        None)
  --storage STORAGE     Request additional storage, in GiB (min 500) (default:
                        None)
  --ulimits U [U ...]   Change instance ulimits, e.g. nofile:1000000 (default:
                        None)
  --environment ENV [ENV ...]
                        Set environment variables (default: None)

other options:
  --dryrun              Print the command but don't launch the job (default:
                        False)
  --s3_script_bucket S3_SCRIPT_BUCKET
                        S3 bucket containing scripts (default: czbiohub-
                        scripts)
  --s3_script_dir S3_SCRIPT_DIR
                        Path to scripts in S3 bucket, if any (default: None)
  -u, --upload          Upload the script to S3 before running (default:
                        False)
  -d, --debug           Set logging to debug level (default: False)
  -h, --help            show this help message and exit

See https://github.com/czbiohub/utilities for more examples
```

Scripts are located inside this repository. Scripts are referred to using *Python import syntax* relative to the utilities package&mdash;e.g. `demux.bcl2fastq` or `alignment.10x_count`. Options for `evros` go before the script name, and options for the script itself go after. `evros` will read that script and test the arguments you are passing (by calling a function named `get_parser` in the script) to make sure they are acceptable.

If the script defines a function named `get_default_requirements` it will call that function to set instance requirements for your job, so you do not need to specify them.

If you write custom scripts that follow these conventions, you can put them in `utilities/utilities/custom` and run them with `evros custom.your_script`. A template script is included as an example.

### How to demux a thing using the standard workflow:

If your BCLs were uploaded to czbiohub-seqbot/bcl and your sample sheet is in czbiohub-seqbot/sample-sheets, you can do this:

```
(utilities-env) ➜ evros demux.bcl2fastq --exp_id YYMMDD_EXP_ID
```

If you want to stick the results somewhere other than czbiohub-seqbot/fastqs, you can change that option:

```
(utilities-env) ➜ evros demux.bcl2fastq --exp_id YYMMDD_EXP_ID --s3_output_dir s3://my-special-bucket
```

### How to align some stuff:

```
(utilities-env) ➜ aws_star mus 10 YYMMDD_EXP_ID > my_star_jobs.sh
(utilities-env) ➜ cat my_star_jobs.sh
evros alignment.run_star_and_htseq --taxon mus --num_partitions 10 --partition_id 0 --exp_ids YYMMDD_EXP_ID
sleep 10
[...lots more...]
(utilities-env) ➜ source my_star_jobs.sh
```

#### How to check for failed alignment jobs:

For some reason, a fraction of alignment jobs fail to start because of AWS problems. It happens enough that there's a script to help with the problem:

```
(utilities-env) ➜ starfails my_star_jobs.sh
8d920e9f-313a-465a-ae4d-df77bdbe990d
(utilities-env) ➜ cat my_star_jobs_failed_jobs.sh 
evros alignment.run_star_and_htseq --taxon mus --num_partitions 10 --partition_id 4 --exp_ids YYMMDD_EXP_ID
sleep 20
(utilities-env) ➜ source my_star_jobs_failed_jobs.sh 
```

This new file contains the command to re-try the failed jobs.

### How to make a gene-cell table from an alignment:

This one runs on your local machine--it'll download alignment results from S3 and make a table out of it.

```
(utilities-env) ➜ gene_cell_table --help
usage: gene_cell_table [--s3_bucket S3_BUCKET] [--dryrun] [--debug] [-h]
                       s3_path output_file

Construct the gene-cell table for an experiment e.g. gene_cell_table
--s3_bucket czbiohub-maca data/YYMMDD_EXP_ID path/to/output.csv

basic arguments:
  s3_path               Path to experiment. e.g.
                        fastqs/YYMMDD_EXP_ID
  output_file           File to save the output, e.g. my_gc_table.csv

other options:
  --s3_bucket S3_BUCKET
                        S3 bucket. e.g. czbiohub-seqbot (default: czbiohub-
                        seqbot)
  --dryrun              Don't actually download any files (default: False)
  --debug               Set logging to debug level (default: False)
  -h, --help            Show this help message and exit

See https://github.com/czbiohub/utilities for more examples

(utilities-env) ➜ gene_cell_table fastqs/YYMMDD_EXP_ID YYMMDD_EXP_ID.csv --dryrun
2017-11-08 18:19:22,494 - __main__ - INFO - (DRYRUN) - Starting
2017-11-08 18:19:22,494 - __main__ - INFO - (DRYRUN) - Starting S3 client
2017-11-08 18:19:22,594 - __main__ - INFO - (DRYRUN) - Getting htseq file list
2017-11-08 18:19:23,176 - __main__ - INFO - (DRYRUN) - 19 htseq files found
2017-11-08 18:19:23,176 - __main__ - INFO - (DRYRUN) - Downloaded 19 files
2017-11-08 18:19:23,177 - __main__ - INFO - (DRYRUN) - Writing to YYMMDD_EXP_ID.csv
2017-11-08 18:19:23,177 - __main__ - INFO - (DRYRUN) - Downloaded 19 files
2017-11-08 18:19:23,177 - __main__ - INFO - (DRYRUN) - Writing to YYMMDD_EXP_ID.log.csv
2017-11-08 18:19:23,177 - __main__ - INFO - (DRYRUN) - Done!

```

### How to share data on s3://czbiohub-seqbot with an outside collaborator

They don't need an AWS account but they _do_ need to install the [AWS CLI](aws.amazon.com/cli).

This script generates an access token for your collaborator to download from our S3 storage, and writes instructions to the console (the command below will redirect the instructions to a text file you can email to your collaborator). The token is good for 36 hours&mdash;if they need more time, just generate another one.

To use the script, just give it the path to the folder you want to share (not including the bucket name, as shown below)

**Note:** This script _only_ applies to the `czbiohub-seqbot` bucket. In the future we'll hopefully have more general tools for sharing our data.

```
(utilities-env) ➜ aws_access fastqs/YYMMDD_EXP_ID > download_token.txt
```
