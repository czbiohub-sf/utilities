import os

import boto3

from boto3.s3.transfer import TransferConfig

import multiprocessing


# cribbed from https://github.com/chanzuckerberg/s3mi/blob/master/scripts/s3mi
def s3_bucket_and_key(s3_uri):
    prefix = "s3://"
    assert s3_uri.startswith(prefix)

    return s3_uri[len(prefix) :].split("/", 1)


def prefix_gen(bucket, prefix, fn=None):
    """Generic generator of fn(result) from an S3 paginator"""
    client = boto3.client("s3")
    paginator = client.get_paginator("list_objects")

    response_iterator = paginator.paginate(Bucket=bucket, Prefix=prefix)

    for result in response_iterator:
        if "Contents" in result:
            yield from (fn(r) for r in result["Contents"])


def get_files(bucket="czb-seqbot", prefix=None):
    """Generator of keys for a given S3 prefix"""
    yield from prefix_gen(bucket, prefix, lambda r: r["Key"])


def get_size(bucket="czb-seqbot", prefix=None):
    """Generator of (key,size) for a given S3 prefix"""
    yield from prefix_gen(bucket, prefix, lambda r: (r["Key"], r["Size"]))


def get_status(file_list, bucket_name="czb-seqbot"):
    """Print the storage/restore status for a list of keys"""
    s3res = boto3.resource("s3")

    for fn in file_list:
        obj = s3res.Object(bucket_name, fn)
        print(obj.key, obj.storage_class, obj.restore)


def restore_file(k):
    obj = s3r.Object("czbiohub-seqbot", k)
    if obj.storage_class == "GLACIER" and not obj.restore:
        bucket.meta.client.restore_object(
            Bucket="czbiohub-seqbot", Key=k, RestoreRequest={"Days": 7}
        )


def restore_files(file_list, *, n_proc=16):
    """Restore a list of files from czbiohub-seqbot in parallel"""

    global s3r
    s3r = boto3.resource("s3")
    global bucket
    bucket = s3r.Bucket("czbiohub-seqbot")

    print("creating pool")
    p = multiprocessing.Pool(processes=n_proc)

    print("restoring files...")
    p.map(restore_file, file_list, chunksize=64)

    p.close()
    p.join()


def copy_file(key, new_key):
    s3c.copy(
        CopySource={"Bucket": bucket, "Key": key},
        Bucket=new_bucket,
        Key=new_key,
        Config=TransferConfig(use_threads=False),
    )


def copy_files(src_list, dest_list, *, b, nb, force_copy=False, n_proc=16):
    """
    Copy a list of files from src_list to dest_list.
    b - original bucket
    nb - destination bucket
    """

    global s3c
    s3c = boto3.client("s3")

    global bucket
    bucket = b
    global new_bucket
    new_bucket = nb

    if not force_copy:
        existing_files = set(
            get_files(bucket=nb, prefix=os.path.commonprefix(dest_list))
        ) & set(dest_list)
    else:
        existing_files = set()

    print("creating pool")
    p = multiprocessing.Pool(processes=n_proc)

    print("copying files")
    p.starmap(
        copy_file,
        (
            (src, dest)
            for src, dest in zip(src_list, dest_list)
            if dest not in existing_files
        ),
        chunksize=64,
    )

    p.close()
    p.join()


def remove_file(k):
    s3c.delete_object(Bucket=bucket, Key=k)


def remove_files(file_list, *, b, really=False, n_proc=16):
    """Remove a list of file keys from S3"""

    assert really

    global s3c
    s3c = boto3.client("s3")

    global bucket
    bucket = b

    print("creating pool")
    p = multiprocessing.Pool(processes=n_proc)

    print("Removing {} files!".format(len(file_list)))
    p.map(remove_file, file_list, chunksize=64)

    p.close()
    p.join()


def download_file(key, dest):
    s3c.download_file(
        Bucket=bucket, Key=key, Filename=dest, Config=TransferConfig(use_threads=False)
    )


def download_files(src_list, dest_list, *, b, force_download=False, n_proc=16):
    """Download a list of file to local storage"""

    global s3c
    s3c = boto3.client("s3")

    global bucket
    bucket = b

    if not force_download:
        existing_files = set(fn for fn in dest_list if os.path.exists(fn))
    else:
        existing_files = set()

    p = multiprocessing.Pool(processes=n_proc)

    # downloading files
    p.starmap(
        download_file,
        (
            (src, dest)
            for src, dest in zip(src_list, dest_list)
            if dest not in existing_files
        ),
        chunksize=64,
    )

    p.close()
    p.join()
