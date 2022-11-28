#!/bin/bash

aws s3 sync --profile czb "resources_test" "s3://czbiohub-pipelines/utilities/" --delete --dryrun