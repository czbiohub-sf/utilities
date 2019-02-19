#!/usr/bin/env python

import argparse
import os

import aegea.util.aws.clients as clients


def main():
    parser = argparse.ArgumentParser(
        prog="starfails",
        description="Search for failed STAR jobs and create a new file to retry them",
    )

    parser.add_argument("job_file")

    args = parser.parse_args()

    with open(args.job_file) as f:
        jobs = {
            line.strip().split("run_star_and_htseq")[-1]: line.strip()
            for line in f
            if line.find("evros") > -1
        }

    failed_job_ids = [
        job_info["jobId"]
        for job_info in clients.batch.list_jobs(
            jobQueue="aegea_batch", jobStatus="FAILED"
        )["jobSummaryList"]
    ]

    job_cmds = {
        job_desc["container"]["command"][-1].split("run_star_and_htseq")[-1]
        for job_desc in clients.batch.describe_jobs(jobs=failed_job_ids)["jobs"]
    }

    failed_cmds = []
    for job_cmd in job_cmds:
        if job_cmd in jobs:
            failed_cmds.append(jobs[job_cmd])

    if failed_cmds:
        print(f"{len(failed_cmds)} jobs failed :(")
        with open("_failed".join(os.path.splitext(args.job_file)), "w") as out:
            for failed_cmd in failed_cmds:
                print(failed_cmd, file=out)
                print("sleep 20", file=out)
    else:
        print("Looks like nothing failed!")
