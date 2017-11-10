# utilities
A collection of scripts for common data management and processing tasks

## Installation

Clone this report and make a new environment, then install the package:

Create a conda environment called `utilities-env` specific to this Python package:

```
➜  conda create --name utilities-env python pip
➜  ~ source activate utilities-env
(utilities-env) ➜  ~ 
```
The `(utilities-env)` indicates that the `utilities-env` environment is active.

Now let's clone the environment into a `code` folder:

```
(utilities-env) ➜  ~ cd code 
(utilities-env) ➜ code git clone https://github.com/czbiohub/utilities.git
# Create an environment called `utilities-env` with Python and pip
```

Now we'll change to that folder with `cd` and install the program with `pip install -e .`, where the `.` means right here.

```
(utilities-env) ➜  cd utilities
(utilities-env) ➜  utilities git:(master) pip install -e .
```

## Usage

### How to run a script on AWS, in general:

The `aegea_launcher.py` script just wraps aegea (which can be tricky to use). It has sane defaults so you don't need to choose most of them.

```
➜  utilities git:(master) ./aegea_launcher.py --help 
usage: aegea_launcher.py [--ecr-image ECR | --ami AMI] [--queue QUEUE]
                         [--vcpus VCPUS] [--memory MEMORY] [--storage STORAGE]
                         [--ulimits U [U ...]] [--environment ENV [ENV ...]]
                         [--dryrun] [-u] [-d] [-h]
                         s3_script_path script_name script_args

Run any script as a batch job e.g. aegea_launcher.py my_bucket/my_scripts
[script name] [script args...]

[...more stuff...]
```

### How to demux a thing using the standard workflow:

If your BCLs were uploaded to czbiohub-seqbot/bcl and your sample sheet is in czbiohub-seqbot/sample-sheets, you can do this:

```
➜  utilities git:(master) ./aegea_launcher.py --queue aegea_batch_demux --vcpus 16 --memory 128000 --storage 1000 jamestwebber-logs/new_scripts bcl2fastq.py "--exp_id 171103_M05295_0051_000000000-BDHHB"
```

If you want to stick the results somewhere other than czbiohub-seqbot/fastqs, you can change that option:

```
➜  utilities git:(master) ./aegea_launcher.py --queue aegea_batch_demux --vcpus 16 --memory 128000 --storage 1000 jamestwebber-logs/new_scripts bcl2fastq.py "--exp_id 171103_M05295_0051_000000000-BDHHB --s3_output_dir s3://my-special-bucket"
```


### How to align some stuff:

```
➜  utilities git:(master) ./alignment/aws_star.py mus 10 171101_NB501961_0026_AHL33MBGX3 > giana_star.sh
➜  utilities git:(master) cat giana_star.sh 
python aegea_launcher.py --queue aegea_batch --vcpus 16 --memory 64000 --storage 500 jamestwebber-logs/new_scripts run_star_and_htseq.py "--taxon mus --num_partitions 10 --partition_id 0 --exp_ids 171101_NB501961_0026_AHL33MBGX3"
sleep 10
[...lots more...]
➜  utilities git:(master) source giana_star.sh 
```


### How to make a gene-cell table from an alignment:

This one runs on your local machine--it'll download alignment results from S3 and make a table out of it.

```
➜  alignment git:(master) ./gene_cell_table.py --help
usage: gene_cell_table.py [--s3_bucket S3_BUCKET] [--dryrun] [--debug] [-h]
                          s3_path output_file

Construct the gene-cell table for an experiment e.g.
./generate_gene_cel_table.py my_bucket/my_scripts [script name] [script
args...]

basic arguments:
  s3_path               Path to experiment. e.g.
                        fastqs/171101_NB501961_0026_AHL33MBGX3
  output_file           File to save the output, e.g. my_gc_table.csv

other options:
  --s3_bucket S3_BUCKET
                        S3 bucket. e.g. czbiohub-seqbot
  --dryrun              Don't actually download any files
  --debug               Set logging to debug level
  -h, --help            Show this help message and exit

See https://github.com/czbiohub/utilities for more examples

➜  alignment git:(master) ./gene_cell_table.py fastqs/171101_NB501961_0026_AHL33MBGX3 171101_NB501961_0026_AHL33MBGX3.csv --dryrun 
2017-11-08 18:19:22,494 - __main__ - INFO - (DRYRUN) - Starting
2017-11-08 18:19:22,494 - __main__ - INFO - (DRYRUN) - Starting S3 client
2017-11-08 18:19:22,594 - __main__ - INFO - (DRYRUN) - Getting htseq file list
2017-11-08 18:19:23,176 - __main__ - INFO - (DRYRUN) - 19 htseq files found
2017-11-08 18:19:23,176 - __main__ - INFO - (DRYRUN) - Downloaded 19 files
2017-11-08 18:19:23,177 - __main__ - INFO - (DRYRUN) - Writing to 171101_NB501961_0026_AHL33MBGX3.csv
2017-11-08 18:19:23,177 - __main__ - INFO - (DRYRUN) - Downloaded 19 files
2017-11-08 18:19:23,177 - __main__ - INFO - (DRYRUN) - Writing to 171101_NB501961_0026_AHL33MBGX3.log.csv
2017-11-08 18:19:23,177 - __main__ - INFO - (DRYRUN) - Done!

```

