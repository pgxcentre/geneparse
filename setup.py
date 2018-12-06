#!/usr/bin/env python

# How to build source distribution
#   - python setup.py sdist --format bztar
#   - python setup.py sdist --format gztar
#   - python setup.py sdist --format zip
#   - python setup.py bdist_wheel


import os
import sys

from setuptools import setup, find_packages


MAJOR = 0
MINOR = 7
MICRO = 5
VERSION = "{0}.{1}.{2}".format(MAJOR, MINOR, MICRO)


def check_python_version():
    """Checks the python version, exits if < 3.4."""
    python_major, python_minor = sys.version_info[:2]

    if python_major != 3 or python_minor < 4:
        sys.stderr.write("geneparse requires python 3 "
                         "(version 3.4 or higher)\n")
        sys.exit(1)


def write_version_file(fn=None):
    if fn is None:
        fn = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            os.path.join("geneparse", "version.py"),
        )

    content = ("\n# THIS FILE WAS GENERATED AUTOMATICALLY\n"
               'geneparse_version = "{version}"\n')

    a = open(fn, "w")
    try:
        a.write(content.format(version=VERSION))
    finally:
        a.close()


def setup_package():
    # Checking the python version prior to installation
    check_python_version()

    # Saving the version into a file
    write_version_file()

    setup(
        name="geneparse",
        version=VERSION,
        description="A suite of parse for genotype formats.",
        url="https://github.com/pgxcentre/geneparse",
        license="MIT",
        test_suite="geneparse.tests.test_suite",
        zip_safe=False,
        install_requires=["numpy >= 1.11.0", "pandas >= 0.19.0",
                          "pyplink >= 1.3.4", "setuptools >= 26.1.0",
                          "biopython >= 1.68", "pybgen >= 0.5.0"],
        packages=find_packages(),
        package_data={"geneparse.tests": ["data/*", "data/*/*"]},
        classifiers=["Development Status :: 4 - Beta",
                     "Intended Audience :: Science/Research",
                     "License :: Free for non-commercial use",
                     "Operating System :: Unix",
                     "Operating System :: POSIX :: Linux",
                     "Operating System :: MacOS :: MacOS X",
                     "Operating System :: Microsoft",
                     "Programming Language :: Python",
                     "Programming Language :: Python :: 3.4",
                     "Programming Language :: Python :: 3.5",
                     "Programming Language :: Python :: 3.6",
                     "Topic :: Scientific/Engineering :: Bio-Informatics"],
        keywords="bioinformatics genetics statistics",
    )


if __name__ == "__main__":
    setup_package()
