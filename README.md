# czb-util
A collection of scripts for common data management and processing tasks

## Quick Reference

| Task | Command | Description |
| ----------- | -------- | ----------- |
| Demux a standard sequencing run | `evros demux.bcl2fastq --exp_id YYMMDD_EXP_ID` | Assumes your sample sheet is uploaded to S3.|
| Demux a 10X run | `evros demux.10x_mkfastq --exp_id YYMMDD_EXP_ID` | Again, assumes a sample sheet is present on S3 |
| Align using 10X (slow) | ```evros alignment.run_10x_count [see --help for input format] --taxon [see --help for options] --s3_input_path s3://input-bucket/path/to/fastqs --s3_output_path s3://output-bucket/path/for/results``` | Run once for each channel of the run. Very slow! |
| Align using 10X (fast) | ```aws_10x [see --help for input format] --taxon TAXON --s3_input_path s3://input-bucket/path/to/sample/folders --s3_output_path s3://output-bucket/path/for/results > your_script.sh``` | Creates a shell script locally to launch many alignments using `source your_script.sh` |
| Align using STAR and htseq (slow) | ```evros alignment.run_star_and_htseq [see --help for input format] --taxon [see --help for options] --num_partitions NUM_PARTITIONS --partition_id PARTITION_ID --s3_input_path s3://input-bucket/path/to/fastqs --s3_output_path s3://output-bucket/path/for/results``` | Run once for each channel of the run. Very slow! |
| Align using STAR and htseq (fast) | ```aws_star [see --help for input format] --taxon TAXON --num_partitions NUM_PARTITIONS --s3_input_path s3://input-bucket/path/to/sample/folders --s3_output_path s3://output-bucket/path/for/results > your_script.sh``` | Creates a shell script locally to launch many alignments using `source your_script.sh` |
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
(utilities-env) ➜ pip install .
```

If you are planning to use `evros` to launch AWS batch jobs, you will need additional dependencies:

```zsh
(utilities-env) ➜ cd utilities
(utilities-env) ➜ pip install .'[evros]'
```

The quotes are necessary if you are using `zsh` as it uses square brackets for pattern matching.

**Note** `pip install .` and `pip install .'[evros]'` do not install the repo as an editable package (unlike previous versions). When you update from Github locally you should reinstall by running the two commands again. The old way was more trouble than it was worth. If you want to install an editable package for the sake of debugging, run the following lines instead:

```zsh
(utilities-env) ➜ cd utilities
(utilities-env) ➜ pip install -e .'[evros]'
```

You can edit existing files and debug correspondingly in terminal with the editable package. However, if you add a new file locally to the repo and want to run it directly in terminal, you should update `entry_points` in `setup.py` and reinstall the package to make the new file run.

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

Scripts are located inside this repository. Scripts are referred to using *Python import syntax* relative to the utilities package&mdash;e.g. `demux.bcl2fastq` or `alignment.run_10x_count`. Options for `evros` go before the script name, and options for the script itself go after. `evros` will read that script and test the arguments you are passing (by calling a function named `get_parser` in the script) to make sure they are acceptable.

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

#### How to choose the right genome (--taxon argument in running the alignment job):

To choose the right genome as the --taxon input, first decide on the species of the sample, then select from the genomes of that species.
The following is our reference_genomes dictionary, with input argument in alignment jobs as keys, and genome name as values:

```zsh
reference_genomes = {
    "homo": "HG38-PLUS",
    "hg38-plus": "HG38-PLUS",
    "homo.gencode.v30.ERCC.chrM": "homo.gencode.v30.annotation.ERCC92",
    "mus": "MM10-PLUS",
    "mm10-plus": "MM10-PLUS",
    "mm10-1.2.0": "mm10-1.2.0",
    "mus-premrna": "mm10-1.2.0-premrna",
    "mm10-1.2.0-premrna": "mm10-1.2.0-premrna",
    "hg19-mm10-3.0.0": "hg19-mm10-3.0.0",
    "microcebus": "MicMur3-PLUS",
    "gencode.vM19": "gencode.vM19",
    "GRCh38_premrna": "GRCh38_premrna",
    "zebrafish-plus": "danio_rerio_plus_STAR2.6.1d"
}
```  

The genomes above span four species: human, mouse, mouse lemur, and zebrafish.
Human genomes are shown below:
```zsh
human_genomes = {
    "homo": "HG38-PLUS"
    "hg38-plus": "HG38-PLUS"
    "homo.gencode.v30.ERCC.chrM": "homo.gencode.v30.annotation.ERCC92"
    "GRCh38_premrna": "GRCh38_premrna"
}
```

The first and second keys point to the same genome, and that's because the first one - "homo" - is deprecated. Other than that, genome "HG38-PLUS" is taken from the homo sapiens genome from UCSC, available on the website of [Illumina iGenomes](https://support.illumina.com/sequencing/sequencing_software/igenome.html). The third genome, "homo.gencode.v30.annotation.ERCC92", is taken from [GENCODE Human Release 30](https://www.gencodegenes.org/human/release_30.html). "chrM" in the key means mitochondria chromosomes - thus genes - are included in this genome, and "[ERCC92](https://www.thermofisher.com/order/catalog/product/4456740#/4456740)" means external controls of RNA variations are applied, using 92 transcripts. The last genome, "GRCh38_premrna", is taken from [Ensembl](https://uswest.ensembl.org/Homo_sapiens/Info/Index), where "premrna" means the genome includes unspliced RNA.

Then we move to mouse genomes:
```zsh
mouse_genomes = {
    "mus": "MM10-PLUS"
    "mm10-plus": "MM10-PLUS"
    "mm10-1.2.0": "mm10-1.2.0"
    "mus-premrna": "mm10-1.2.0-premrna"
    "mm10-1.2.0-premrna": "mm10-1.2.0-premrna"
    "gencode.vM19": "gencode.vM19"
}
```

Again, the first and fourth keys are deprecated. "mm10-plus" is now the key to get genome "MM10-PLUS", and "mm10-1.2.0-premrna" for genome "mm10-1.2.0-premrna". "MM10-PLUS" is taken from the *Mus musculus* (Mouse) genome from UCSC, available on the website of [Illumina iGenomes](https://support.illumina.com/sequencing/sequencing_software/igenome.html). Both "mm10-1.2.0" and "mm10-1.2.0-premrna" are taken from [Ensembl Release 84](https://support.10xgenomics.com/single-cell-gene-expression/software/release-notes/build#mm10_1.2.0), and they're available for download following the instructions on 10x Genomics website. Differences among the above genomes are that "MM10-PLUS" includes a few transgenes (genes inserted and transcribed), "mm10-1.2.0" has better genomic and mitochondria annotation, and "mm10-1.2.0-premrna" further includes unspliced RNA. The last genome "gencode.vM19" is taken from [GENCODE Mouse Release 19](https://www.gencodegenes.org/mouse/release_M19.html).

To wrap up human and mouse sections, `"[hg19-mm10-3.0.0": "hg19-mm10-3.0.0](https://support.10xgenomics.com/single-cell-gene-expression/software/release-notes/build#hg19mm10_3.0.0)"` is the human and mouse genome, downloadable also on 10X Genomics website.

Now we have two genomes left. "danio_rerio_plus_STAR2.6.1d" is the danio rerio (zebrafish) genome taken from [Ensembl](https://uswest.ensembl.org/Danio_rerio/Info/Index). "MicMur3-PLUS" is the mouse lemur genome taken from [Ensembl](https://uswest.ensembl.org/Microcebus_murinus/Info/Index), and you can browse the genome on the website of [UCSC](http://genome-preview.cse.ucsc.edu/cgi-bin/hgTracks?db=micMur3&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr17%3A30720951%2D30776185&hgsid=394238320_4HFkLIDhzsNba96BNgHUlsvY5Ldc).

After selecting a proper genome, you can go ahead running the alignment job.

#### How to launch an alignment job for one sample folder:

*Important change: you need to explicitly specify the input and output paths for your alignment.*

The following shows an example to align a sample using STAR and htseq. Don't forget to include "--taxon", "--num_partitions", "--partition_id", "--s3_input_path", and "--s3_output_path" in the command since the arguments are not positional, and their sequence doesn't matter. When running only one sample folder, just enter 10 for --num_partitions and 0 for --partition_id as default values, since these two arguments are only of interest in running multiple sample folders all together. See documentation in "How to launch alignment job for multiple sample folders" below for more details:

```zsh
(utilities-env) ➜ evros alignment.run_star_and_htseq --taxon homo.gencode.v30.ERCC.chrM --num_partitions 10 --partition_id 0 --s3_input_path s3://tabula-sapiens/Pilot1/fastqs/smartseq2/pilot/B107809_A10_S130 --s3_output_path s3://output-bucket/path/for/results
```

This will find all `.fastq.gz` files under the given input path and align them to to the homo.gencode.v30.ERCC.chrM genome.

Likewise, to align a 10x run, take the following code as an example:

```zsh
(utilities-env) ➜ evros alignment.run_10x_count --taxon homo.gencode.v30.ERCC.chrM --num_partitions 10 --partition_id 0 --s3_input_path s3://tabula-sapiens/Pilot2/fastqs/10X/NovaSeq_ReRun_Folders/TSP2_BM_vertebralbody_10X_1_1 --s3_output_path s3://output-bucket-name/path/for/results
```

Running either `run_star_and_htseq` or `run_10x_count` launches an alignment job on AWS Batch, and you'll get the jobID from the terminal. Copy and paste that jobID into the search box on AWS Batch Job page, and the job you just launched will pop out. You don't need to do anything else there - the job status will turn from runnable to starting and then running automatically. Usually alignment using `run_star_and_htseq` takes around half an hour, and alignment using `run_10x_count` takes hours.

#### How to launch alignment jobs for multiple sample folders:

You can use the helper script below, replacing inputs with your values, to create a bunch of commands to efficiently align using STAR and htseq, instead of running one sample folder at a time. Argument `num_partitions` is the number of jobs you want to launch to align all the sample reads under the input folder, and `partition_id` tags each job with integers from 0 to `num_partitions - 1`. For example, if there're 200 sample folders under the input folder, and we enter 10 for `num_partitions`, then each job runs 20 sample folders. The job with `partition_id` 0 runs 1st, 11th, 21st, ..., 191th sample folders, and likewise for jobs with other `partition_id`.

```zsh
(utilities-env) ➜ aws_star --taxon homo.gencode.v30.ERCC.chrM --num_partitions 10 --s3_input_path s3://tabula-sapiens/Pilot1/fastqs/smartseq2/pilot --s3_output_path s3://output-bucket/path/for/results > my_star_jobs.sh
(utilities-env) ➜ cat my_star_jobs.sh
evros --branch master alignment.run_star_and_htseq --taxon homo.gencode.v30.ERCC.chrM --num_partitions 10 --partition_id 0 --s3_input_path s3://tabula-sapiens/Pilot1/fastqs/smartseq2/pilot --s3_output_path s3://output-bucket/path/for/results
sleep 10
[...lots more, with increasing partition_id totaling num_partitions...]
(utilities-env) ➜ source my_star_jobs.sh
```

Likewise, to align multiple 10x runs, take the following code as an example:

```zsh
(utilities-env) ➜ aws_10x --taxon homo.gencode.v30.ERCC.chrM --s3_input_path s3://tabula-sapiens/Pilot2/fastqs/10X/NovaSeq_ReRun_Folders --s3_output_path s3://output-bucket/path/for/results > my_10x_jobs.sh
(utilities-env) ➜ cat my_10x_jobs.sh
evros --branch master alignment.run_10x_count --taxon homo.gencode.v30.ERCC.chrM --num_partitions 35 --partition_id 0 --s3_input_path s3://tabula-sapiens/Pilot2/fastqs/10X/NovaSeq_ReRun_Folders/TSP2_BM_vertebralbody_10X_1_1/ --s3_output_path s3://output-bucket/path/for/results/TSP2_BM_vertebralbody_10X_1_1
sleep 10
[...lots more, with increasing partition_id totaling num_partitions, which is the total number of sample folders under the input path as `run_10x_count` takes one sample folder at a time...]
(utilities-env) ➜ source my_10x_jobs.sh
```

Running either `aws_star` or `aws_10x` will launch multiple jobs all at once for all the sample folders under the input folder, which saves much time. Again, copy the jobIDs printed in the terminal and paste them on AWS Batch Job page will give you information about the status of all those jobs.


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

## Troubleshooting

### AWS Configure

It's necesary to configure the AWS Command Line Interface (AWS CLI) to launch AWS Batch jobs. Otherwise, you may get the following error:

``` zsh
ERROR:aegea:Failed to install Lambda helper:
ERROR:aegea:NoRegionError: You must specify a region.
ERROR:aegea:Aegea will be unable to look up logs for old Batch jobs.
The AWS CLI is not configured. Please configure it using instructions at http://docs.aws.amazon.com/cli/latest/userguide/cli-chap-getting-started.html
```
To fix the error, go to the website for detailed instructions. Sample codes from the website is pasted below. Replace the sample values with your own values to configure the AWS CLI:

```zsh
(utilities-env) ➜ aws configure
AWS Access Key ID [None]:
AWS Secret Access Key [None]:
Default region name [None]: us-west-2
Default output format [None]: json
```
### Aegea Version

If you get errors and exceptions similar to those below, there may be problems with aegea version conflicts:

```zsh
botocore.errorfactory.ResourceNotFoundException: An error occurred (ResourceNotFoundException) when calling the GetFunction operation: Function is not found: arn:aws:lambda:us-west-2:[a number series]:function:aegea-dev-process_batch_event
botocore.exceptions.ClientError: An error occurred (AccessDenied) when calling the PutRolePolicy operation: User: arn:aws:iam::[a number series]:user/[AWS IAM user name] is not authorized to perform: iam-PutRolePolicy on resource: role aegea-dev-process_batch_event
chalice.deploy.deployer.ChaliceDeploymentError: ERROR - While deploying your chalice application, received the following error:

  An error occurred (AccessDenied) when calling the PutRolePolicy operation: 
  User: arn:aws:iam::[a number series]:user/[AWS IAM user name] is not authorized to perform: 
  iam-PutRolePolicy on resource: role aegea-dev-process_batch_event
```

To fix the error, try installing aegea version 2.6.9 which is compatible with the codes of this repo:

```zsh
(utilities-env) ➜ pip install aegea==2.6.9
```