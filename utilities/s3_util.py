import boto3
import botocore.exceptions

from boto3.s3.transfer import TransferConfig

import multiprocessing


# cribbed from https://github.com/chanzuckerberg/s3mi/blob/master/scripts/s3mi
def s3_bucket_and_key(s3_uri):
    prefix = "s3://"
    assert s3_uri.startswith(prefix)

    return s3_uri[len(prefix):].split("/", 1)


# function to init a pool process with an s3 client and bucket
def init_resource(b, nb=None):
    global s3c
    s3c = boto3.client('s3')

    global bucket
    bucket = b

    if nb is not None:
        global new_bucket
        new_bucket = nb


# function to init a pool process with an s3 resource and bucket
def init_client(b):
    global s3r
    s3r = boto3.resource('s3')

    global bucket
    bucket = b


def prefix_gen(bucket, prefix, fn=None):
    """Generic generator of fn(result) from an S3 paginator"""
    client = boto3.client('s3')
    paginator = client.get_paginator('list_objects')

    response_iterator = paginator.paginate(
            Bucket=bucket, Prefix=prefix
    )

    for result in response_iterator:
        if 'Contents' in result:
            yield from (fn(r) for r in result['Contents'])


def get_files(bucket='czbiohub-seqbot', prefix=None):
    """Generator of keys for a given S3 prefix"""
    yield from prefix_gen(bucket, prefix, lambda r: r['Key'])


def get_size(bucket='czbiohub-seqbot', prefix=None):
    """Generator of (key,size) for a given S3 prefix"""
    yield from prefix_gen(bucket, prefix, lambda r: (r['Key'], r['Size']))


def get_status(file_list, bucket_name='czbiohub-seqbot'):
    """Print the storage/restore status for a list of keys"""
    s3r = boto3.resource('s3')

    for fn in file_list:
        obj = s3r.Object(bucket_name, fn)
        print(obj.key, obj.storage_class, obj.restore)


def restore_file(k):
    obj = s3r.Object('czbiohub-seqbot', k)
    if not obj.restore:
        bucket.meta.client.restore_object(
            Bucket='czbiohub-seqbot',
            Key=k,
            RestoreRequest={'Days': 3}
        )


def restore_files(file_list, *, b='czbiohub-seqbot', n_proc=16):
    """Restore a list of files from czbiohub-seqbot in parallel"""

    print('creating pool')
    p = multiprocessing.Pool(processes=n_proc,
                             initializer=init_resource,
                             initargs=(b,))

    try:
        print('restoring files...')
        p.map(restore_file, file_list, chunksize=100)
    finally:
        p.close()
        p.join()


def copy_file(k):
    key, new_key = k
    try:
        s3c.head_object(Bucket=new_bucket, Key=new_key)
    except botocore.exceptions.ClientError:
        s3c.copy(CopySource={'Bucket': bucket, 'Key': key},
                 Bucket=new_bucket,
                 Key=new_key,
                 Config=TransferConfig(use_threads=False))


def copy_files(src_list, dest_list, *, b, nb, n_proc=16):
    """
    Copy a list of files from src_list to dest_list.
    b - original bucket
    nb - destination bucket
    """

    p = multiprocessing.Pool(processes=n_proc,
                             initializer=init_client,
                             initargs=(b, nb))

    try:
        print('copying files')
        p.map(copy_file, zip(src_list, dest_list), chunksize=100)
    finally:
        p.close()
        p.join()


def remove_file(k):
    s3c.delete_object(Bucket=bucket, Key=k)


def remove_files(file_list, *, b, really=False, n_proc=16):
    """Remove a list of file keys from S3"""

    assert really

    p = multiprocessing.Pool(processes=n_proc,
                             initializer=init_client,
                             initargs=(b,))

    try:
        print("Removing {} files!".format(len(file_list)))
        p.map(remove_file, file_list, chunksize=100)
    finally:
        p.close()
        p.join()


def download_file(k):
    key, dest = k
    s3c.download_file(Bucket=bucket,
                      Key=key,
                      Filename=dest,
                      Config=TransferConfig(use_threads=False))


def download_files(src_list, dest_list, *, b, n_proc=16):
    """Download a list of file to local storage"""

    p = multiprocessing.Pool(processes=n_proc,
                             initializer=init_client,
                             initargs=(b,))

    try:
        p.map(download_file, zip(src_list, dest_list), chunksize=100)
    finally:
        p.close()
        p.join()
