# utilities 0.1.3

* `mapping/process_smartseq2_auto`: Strip trailing slashes from id.

* `mapping/process_10x_auto`: Strip trailing slashes from id.

* `mapping/process_10x_auto`: Fix regex not detecting fastq files not containing `_L001`.

* Switch to development build of OpenPipelines (pre0.9.0).

# utilities 0.1.2

* `operations/create_runner_script`: Add workaround for submodule being unintentionally removed when 'nextflow pull' is not run beforehand.

* `mapping/process_smartseq2_auto`: Added above workaround to helper script.

* `mapping/process_10x_auto`: Added above workaround to helper script.

# utilities 0.1.1

* `mapping/process_smartseq2_auto`: Allow grouping input fastqs per folder or via a custom regex statement.

* `mapping/process_smartseq2_auto`: Add more documentation to the arguments to clarify regex usage.

# utilities 0.1.0

First release of the `utilities` as a Viash+Nextflow pipeline.

* `mapping/process_10x`: Mapping 10x samples
* `mapping/process_10x_auto`: Map all 10x Fastq files in a directory.
* `mapping/process_smartseq2`: Mapping smartseq2 samples
* `mapping/process_smartseq2_auto`: Map all smartseq2 Fastq files in a directory.