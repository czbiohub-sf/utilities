import argparse
import csv
import posixpath

from utilities.references import reference_genomes
from utilities.s3_util import s3_sync


def get_base_parser(prog, description):
    parser = argparse.ArgumentParser(
        prog=prog,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=description,
    )

    parser.add_argument(
        "--taxon",
        required=True,
        choices=list(reference_genomes.keys()),
        help="Reference genome for the alignment run",
    )

    parser.add_argument(
        "--run_id",
        required=True,
        help="Name of the folder to write results to"
    )

    parser.add_argument(
        "--s3_libraries_csv_path",
        required=True,
        help="The csv with the s3 paths and metadata needed for cellranger arc count"
    )

    parser.add_argument(
        "--s3_output_path",
        required=True,
        help="The s3 path to store the alignment results",
    )

    parser.add_argument("--root_dir", default="/mnt")

    return parser


#TODO(neevor): Clean up the number of args this function takes.
def process_libraries_file(original_libraries_path, libraries_path, data_dir, logger):
    with open(original_libraries_path, newline='') as csvfile, \
            open(libraries_path, 'w') as new_csv:
        headers = next(csvfile)
        new_csv.write(f"{headers}")

        for row in csv.reader(csvfile):
            s3_path_of_fastqs = row[0]
            sample_id = row[1]
            method = row[-1].replace(" ", "_")
            local_path = data_dir / posixpath.basename(s3_path_of_fastqs) / sample_id / method
            s3_sync(logger, s3_path_of_fastqs, str(local_path))
            row[0] = str(local_path)
            row_values = ",".join(row)
            new_csv.write(f"{row_values}\n")


def get_default_requirements():
    return argparse.Namespace(
        vcpus=64, memory=256000, storage=2000, ecr_image="multiomics"
    )
