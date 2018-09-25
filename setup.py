#!/usr/bin/env python

import io
import glob
import os

import setuptools


def read(*names, **kwargs):
    return io.open(
        os.path.join(os.path.dirname(__file__), *names),
        encoding=kwargs.get("encoding", "utf8"),
    ).read()


setuptools.setup(
    name="czb_util",
    version="0.3.0",
    license="MIT License",
    description="A collection of scripts for some common Biohub tasks",
    long_description=read("README.md"),
    author="James Webber",
    author_email="james.webber@czbiohub.org",
    url="https://github.com/czbiohub/utilities",
    packages=setuptools.find_packages("src"),
    package_dir={"": "src"},
    py_modules=[
        os.path.splitext(os.path.basename(path))[0] for path in glob.glob("src/*.py")
    ],
    include_package_data=True,
    zip_safe=False,
    install_requires=[
        "aegea >= 2.2.3",
        "awscli >= 1.15.41",
        "awscli-cwlogs >= 1.4.4",
        "boto3 >= 1.7.41",
        "click >= 6.7",
    ],
    entry_points={"console_scripts": ["frython = czb"]},
)
