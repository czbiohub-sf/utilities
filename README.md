# czb-util
A collection of scripts for common data management and processing tasks

## Quick Reference

| Task | Command | Description |
| ----------- | -------- | ----------- |
| Demux a standard sequencing run | `evros demux.bcl2fastq --exp_id YYMMDD_EXP_ID` | Assumes your sample sheet is uploaded to S3.|
| Demux a 10X run | `evros demux.10x_mkfastq --exp_id YYMMDD_EXP_ID` | Again, assumes a sample sheet is present on S3 |
| Align using STAR and htseq | ```aws_star [see --help for taxons] [# partitions] --s3_input_path s3://input-bucket/path/to/fastqs --s3_output_path s3://output-bucket/path/for/results > your_script.sh``` | Creates a shell script locally to launch many alignments using `source your_script.sh` |
| Align a 10X run | ```evros alignment.10x_count --taxon [see --help for options] --s3_input_dir s3://input-bucket/path/to/fastqs --s3_output_dir s3://output-bucket/path/for/results``` | Run once for each channel of the run. Very slow! |
| Run Velocyto | ```aws_velocyto [see --help for taxon] [# partitions] --s3_input_path s3://input-bucket/path/to/star_output --s3_output_path s3://output-bucket/path/for/results > your_script.sh``` | Creates a shell script locally to run velocyto on STAR output |
| Create a download token | `aws_access fastqs/YYMMDD_EXP_ID [optional bucket] > download_instructions.txt` | Defaults to the `czb-seqbot` bucket |


## Installation

Clone this repo and install it, preferably in a fresh conda environment.

**Note: this package is written for Python 3, it is not compatible with 2.7**

Create a conda environment called `utilities-env` specific to this Python package:

```zsh
# Create an environment called `utilities-env` with Python and pip
➜  conda create --name utilities-env python=3.6 pip
➜  source activate utilities-env
(utilities-env) ➜
```

The `(utilities-env)` indicates that the `utilities-env` environment is active.

Now let's clone the environment into a `code` folder:

```zsh
(utilities-env) ➜ cd code
(utilities-env) ➜ git clone https://github.com/czbiohub/utilities.git
```

Now we'll change to that folder with `cd` and install the package by running the `setup.py` script.

```zsh
(utilities-env) ➜ cd utilities
(utilities-env) ➜ python setup.py install
```

**Note** This does not install the repo as an editable package (unlike previous versions). When you update from Github you should reinstall by running the command again. The old way was more trouble than it was worth. 

## Usage

As of the latest version of this repository, you *do not* need to be in the `utilities` folder to run scripts, but you *must* activate the environment:

```zsh
➜  source activate utilities-env
(utilities-env) ➜
```

The `evros` command (and others) will be on your path, and it will be able to find the scripts to run on AWS.

### How to run a script on AWS, in general:

The `evros` script wraps aegea (which can be tricky to use). It has sane defaults so you don't need to choose most of them.

```zsh
(utilities-env) ➜ evros --help
usage: evros [--ecr-image ECR] [--queue QUEUE] [--vcpus VCPUS]
             [--memory MEMORY] [--storage STORAGE] [--ulimits U [U ...]]
             [--environment ENV [ENV ...]] [--dryrun] [--branch BRANCH] [-d]
             [-h]
             script_name ...

Run batch jobs on AWS e.g. evros [options] demux.bcl2fastq [script args...]

basic arguments:
  script_name           Local path of the script to run, e.g. demux.bcl2fastq
  script_args           Arguments for the script (everything after
                        script_name)

customize the instance:
  --ecr-image ECR       ECR image to use for the job (default: demuxer)
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
  --branch BRANCH       branch of utilities repo to use (default: master)
  -d, --debug           Set logging to debug level (default: False)
  -h, --help            show this help message and exit

See https://github.com/czbiohub/utilities for more examples
```

Scripts are located inside this repository. Scripts are referred to using *Python import syntax* relative to the utilities package&mdash;e.g. `demux.bcl2fastq` or `alignment.10x_count`. Options for `evros` go before the script name, and options for the script itself go after. `evros` will read that script and test the arguments you are passing (by calling a function named `get_parser` in the script) to make sure they are acceptable.

If the script defines a function named `get_default_requirements` it will call that function to set instance requirements for your job, so you do not need to specify them.

If you write custom scripts that follow these conventions, `evros` will be able to run them. A template script is included as an example. To run custom scripts, use the `--branch` option. First, create a new branch of the repo, then write your script (or modify an existing one). Once you've committed your changes, push them back to this repo. The batch job will run `git checkout [branch]` at runtime.


```zsh
(utilities-env) ➜ git checkout -b my_custom_branch
Switched to a new branch 'my_custom_branch'

...[ make your changes to the code, e.g. create custom.my_custom_script ]...

(utilities-env) ➜ git commit -am "here are all my changes"
(utilities-env) ➜ git push --set-upstream origin my_custom_branch
Total #### (delta ###), reused ### (delta ####)
To github.com:czbiohub/utilities.git
 * [new branch]      my_custom_branch -> my_custom_branch
Branch 'my_custom_branch' set up to track remote branch 'my_custom_branch' from 'origin'.
(utilities-env) ➜ evros --branch my_custom_branch custom.my_custom_script --arg1 --arg2
```

### How to demux something:

*First consider: maybe don't? Demuxing should happen automatically for most runs* 

If your BCLs were uploaded to czb-seqbot/bcl and your sample sheet is in czb-seqbot/sample-sheets, you can do this:

```zsh
(utilities-env) ➜ evros demux.bcl2fastq --exp_id YYMMDD_EXP_ID
```

If you want to stick the results somewhere other than czb-seqbot/fastqs, you can change that option:

```zsh
(utilities-env) ➜ evros demux.bcl2fastq --exp_id YYMMDD_EXP_ID --s3_output_dir s3://my-special-bucket
```

### How to align some stuff:

*Important change: you need to explicitly specify the input and output paths for your alignment.*

To the run the first of ten partitions:

```zsh
(utilities-env) ➜ evros alignment.run_star_and_htseq --taxon mm10-plus --num_partitions 10 --partition_id 0 --s3_input_path s3://input-bucket/path/to/fastqs --s3_output_path s3://output-bucket/path/for/results
```

This will find all `.fastq.gz` files under the given input path and align them to to the MM10-PLUS genome.

You can use this helper script to create a bunch of commands:

```zsh
(utilities-env) ➜ aws_star mm10-plus 10 --s3_input_path s3://input-bucket/path/to/fastqs --s3_output_path s3://output-bucket/path/for/results > my_star_jobs.sh
(utilities-env) ➜ cat my_star_jobs.sh
evros --branch master alignment.run_star_and_htseq --taxon mus --num_partitions 10 --partition_id 0 --s3_input_path s3://input-bucket/path/to/fastqs --s3_output_path s3://output-bucket/path/for/results
sleep 10
[...lots more...]
(utilities-env) ➜ source my_star_jobs.sh
```

#### How to check for failed alignment jobs:

For some reason, a fraction of alignment jobs fail to start because of AWS problems. It happens enough that there's a script to help with the problem:

```zsh
(utilities-env) ➜ starfails my_star_jobs.sh
8d920e9f-313a-465a-ae4d-df77bdbe990d
(utilities-env) ➜ cat my_star_jobs_failed_jobs.sh 
evros --branch master alignment.run_star_and_htseq --taxon mus --num_partitions 10 --partition_id 4 --s3_input_path s3://input-bucket/path/to/fastqs --s3_output_path s3://output-bucket/path/for/results
sleep 20
(utilities-env) ➜ source my_star_jobs_failed_jobs.sh 
```

This new file contains the commands to re-try the failed jobs.

### How to make a gene-cell table from an alignment:

This one runs on your local machine&mdash;it'll download alignment results from S3 and make a table out of it.

```zsh
(utilities-env) ➜ gene_cell_table --help
usage: gene_cell_table [--no_log] [--dryrun] [--debug] [-h]
                       s3_input_path output_file

Construct the gene-cell table for an experiment e.g. gene_cell_table
s3://bucket-name/path/to/results path/to/output.csv

basic arguments:
  s3_input_path  Location of data on S3
  output_file    File to save the output, e.g. my_gc_table[.csv,.h5ad]

other options:
  --no_log       Don't try to download log files (default: False)
  --dryrun       Don't actually download any files (default: False)
  --debug        Set logging to debug level (default: False)
  -h, --help     Show this help message and exit

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

### *New!* How to run Velocyto on some alignments

This script will use the BAM files from a STAR alignment and create loom files using Velocyto. Currently supports hg38-plus and mm10-plus.

To run the first of ten partitions:

```zsh
(utilities-env) ➜ evros alignment.velocyto --taxon mm10-plus --s3_input_path s3://input-bucket/path/to/star_output --s3_output_path s3://output-bucket/path/to/velocyto_loom_files --num_partitions 10 --partition_id 0 --input_dirs YYMMDD_EXP_ID
```

Or use this helper script:

```zsh
(utilities-env) ➜ aws_velocyto mm10-plus 10 --s3_input_path s3://input-bucket/path/to/star_output --s3_output_path s3://output-bucket/path/to/velocyto_loom_files > my_velocyto_jobs.sh
(utilities-env) ➜ cat my_velocyto_jobs.sh
evros --branch master alignment.velocyto --taxon mm10-plus --num_partitions 10 --partition_id 0 --s3_input_path s3://input-bucket/path/to/star_output --s3_output_path s3://output-bucket/path/to/velocyto_loom_files
sleep 10
[...lots more...]
(utilities-env) ➜ source my_velocyto_jobs.sh
```

This will run Velocyto on every BAM file under the input paths 

### How to share data with an outside collaborator

They don't need an AWS account but they _do_ need to install the [AWS CLI](aws.amazon.com/cli).

This script generates an access token for your collaborator to download from our S3 storage, and writes instructions to the console. The command below will redirect the instructions to a shell script you can email to your collaborator, which they can execute with `source`. The token is good for 36 hours&mdash;if they need more time, just generate another one.

To use the script, just give it the path to the folder you want to share (not including the bucket name, as shown below)

By default this will share data in the `czb-seqbot` bucket. If you want to share data from somewhere else, give it the bucket name:

```zsh
(utilities-env) ➜ aws_access fastqs/YYMMDD_EXP_ID > download_script.sh
(utilities-env) ➜ aws_access some/data/here my-own-bucket > download_script_for_my_bucket.sh
```

If you need someone to _upload_ data to our S3 storage, there is a separate script for that, called `aws_upload`, which works similarly.
