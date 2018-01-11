#!/usr/bin/env python

import glob
import os

from setuptools import setup

version = '0.3.0'

required = open('requirements.txt').read().split('\n')
with open(os.path.join(os.path.dirname(__file__), 'requirements.txt')) as f:
    install_requires = [line.rstrip() for line in f]


setup(
    name='utilities',
    version=version,
    description='A collection of scripts for some common Biohub tasks',
    author='James Webber',
    author_email='james.webber@czbiohub.org',
    url='https://github.com/czbiohub/utilities',
    packages=['utilities'],
    install_requires=install_requires,
    scripts=glob.glob('scripts/*'),
    long_description='See https://github.com/czbiohub/utilities',
    license=open("LICENSE").readline().strip(),
)
