#!/usr/bin/env python

from setuptools import setup

version = '0.3.0'

required = open('requirements.txt').read().split('\n')

setup(
    name='utilities',
    version=version,
    description='A collection of scripts for some common Biohub tasks',
    author='James Webber',
    author_email='james.webber@czbiohub.org',
    url='https://github.com/czbiohub/utilities',
    packages=['utilities'],
    install_requires=required,
    long_description='See ' + 'https://github.com/czbiohub/utilities',
    license=open("LICENSE").readline().strip(),
)
