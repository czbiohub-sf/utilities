#!/usr/bin/env python

import argparse
import csv
import os

rc_d = {"A": "T", "G": "C", "C": "G", "T": "A"}


def batch_samplesheet(
    samplesheet_file,
    run_prefix,
    exp_id,
    n,
    reverse_comp_i7,
    reverse_comp_i5,
    s3_input_dir,
    s3_output_dir,
    s3_report_dir,
    s3_sample_sheet_dir,
    star_structure,
):
    """
    samplesheet_file - the giant samplesheet (ideally with the right indexes now)
    run_prefix - shorthand for the run, I usually use something like YYMMDD_A00111
    exp_id - the sequencing run(s) to demux
    n - number of samples to put in each batch. 300 seems to run and doesn't blow up
    reverse_comp_i7 - whether to reverse-complement the first index (it happens)
    reverse_comp_i5 - whether to reverse-complement the second index (for NextSeq runs)
    s3_*_dir - parameters for bcl2fastq.py
    """

    with open(samplesheet_file) as f:
        rdr = csv.reader(f)
        rows = list(rdr)
        hdr = "\n".join(",".join(r) for r in rows[:2])
        if reverse_comp_i7:
            i7_c = list(map(str.lower, rows[1])).index("index")
        if reverse_comp_i5:
            i5_c = list(map(str.lower, rows[1])).index("index2")

        rows = rows[2:]

    print(len(rows), "rows")

    for r in rows:
        if reverse_comp_i7:
            r[i7_c] = "".join(rc_d[nt] for nt in r[i7_c][::-1])
        if reverse_comp_i5:
            r[i5_c] = "".join(rc_d[nt] for nt in r[i5_c][::-1])

    run_prefix_dir = os.path.join(os.path.dirname(samplesheet_file), run_prefix)

    if not os.path.exists(run_prefix_dir):
        os.mkdir(run_prefix_dir)

    # need to upload this whole folder to the sample sheets directory on S3
    for i in range(0, len(rows) + int(len(rows) % n > 0), n):
        with open(os.path.join(run_prefix_dir, f"batch_{i}.csv"), "w") as OUT:
            print(hdr, file=OUT)
            for r in rows[i : i + n]:
                print(",".join(r), file=OUT)

    # print all the run commands to a file so you can just source it
    with open(f"{os.path.dirname(samplesheet_file)}/{run_prefix}.sh", "w") as OUT:
        for i in range(0, len(rows) + int(len(rows) % n > 0), n):
            for run in exp_id:
                print(
                    (
                        f"evros demux.bcl2fastq"
                        f" --exp_id {run}"
                        f" --s3_input_dir {s3_input_dir}"
                        f" --s3_output_dir {s3_output_dir}"
                        f" --s3_report_dir {s3_report_dir}/{run}/batch_{i}"
                        f" --s3_sample_sheet_dir {s3_sample_sheet_dir}/{run_prefix}"
                        f" --sample_sheet_name batch_{i}.csv"
                        " --skip_undetermined"
                        f' {"--star_structure" if star_structure else ""}'
                    ),
                    file=OUT,
                )

    print(
        f"""To run the batch:
    1. Upload the samplesheet folder to S3
        aws s3 cp --recursive {run_prefix_dir} {s3_sample_sheet_dir}/{run_prefix}
    2. Run the shell script of commands
        source {os.path.dirname(samplesheet_file)}/{run_prefix}.sh
    """
    )


def main():
    parser = argparse.ArgumentParser(
        prog="batched_samplesheet",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "samplesheet_file", help="Local copy of a samplesheet to batch up"
    )

    parser.add_argument("--exp_id", nargs="+", help="Run ID(s) to demux")

    parser.add_argument(
        "--run_prefix",
        default=None,
        help="Name for this batch (defaults to [Date]_[SeqId]",
    )
    parser.add_argument(
        "--n", type=int, default=300, help="Number of samples per batch"
    )
    parser.add_argument(
        "--reverse_comp_i7",
        action="store_true",
        help="Reverse-complement the i7 indexes",
    )
    parser.add_argument(
        "--reverse_comp_i5",
        action="store_true",
        help="Reverse-complement the i5 indexes",
    )

    bcl2fastq_options = parser.add_argument_group("bcl2fastq options")
    bcl2fastq_options.add_argument("--s3_input_dir", default="s3://czbiohub-seqbot/bcl")
    bcl2fastq_options.add_argument("--s3_output_dir", default="s3://czb-seqbot/fastqs")
    bcl2fastq_options.add_argument("--s3_report_dir", default="s3://czb-seqbot/reports")
    bcl2fastq_options.add_argument(
        "--s3_sample_sheet_dir", default="s3://czb-seqbot/sample-sheets"
    )
    bcl2fastq_options.add_argument("--star_structure", action="store_true")

    args = parser.parse_args()

    if args.run_prefix is None:
        args.run_prefix = args.exp_id[0].rsplit("_", 2)[0]

    batch_samplesheet(
        args.samplesheet_file,
        args.run_prefix,
        args.exp_id,
        args.n,
        args.reverse_comp_i7,
        args.reverse_comp_i5,
        args.s3_input_dir,
        args.s3_output_dir,
        args.s3_report_dir,
        args.s3_sample_sheet_dir,
        args.star_structure,
    )
