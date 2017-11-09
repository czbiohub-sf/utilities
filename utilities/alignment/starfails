#!/bin/sh
if [ $# -lt 1 ]; then
  echo "Usage: $0 <job file>"
  exit 1
fi

for job in `aegea batch ls --status FAILED -c jobId | grep "minutes ago" | cut -c 4-39 | grep -a "-"`
    do echo $job
    aegea batch describe $job | grep "PARTITION_ID=" | awk '{print $13}' >> ${1}_failed_jobs
done

if [ -f ${1}_failed_jobs ]
then
    grep -f ${1}_failed_jobs $1| awk '{print $0, "\nsleep 20"}' > ${1}_failed_jobs.sh
    rm ${1}_failed_jobs
else
    echo "looks like nothing failed!"
fi
