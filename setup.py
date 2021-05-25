#!/usr/bin/env python
# -*- coding: utf-8

from setuptools import find_packages
from setuptools import setup

NAME = "profile_binr"
VERSION = "0.1.0"

setup(
    name=NAME,
    version=VERSION,
    author="Gustavo Magaña López",
    author_email="gustavo.magana-lopez@u-psud.fr",
    url="https://github.com/bnediction/profile_binr",
    packages=[NAME],
    # packages=find_packages() <- should we use this ?,
    license="BSD-3",
    classifiers=[
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    keywords='',
    description="PROFILE methodology for the binarisation and normalisation of RNA-seq data",
    long_description=open('README.md').read(),
    long_description_content_type="text/markdown",
    install_requires=[
        "pandas",
        "rpy2",
    ],
    package_data={'_R': ['_R/*.R']}
)

