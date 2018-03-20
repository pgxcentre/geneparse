"""Simple script to index genotype files of different format."""

# This file is part of geneparse.
#
# The MIT License (MIT)
#

# Copyright (c) 2017 Pharmacogenomics Centre
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.


import sys
import logging
import argparse
import subprocess

from .impute2 import generate_index as impute2_index


# Logging configuration
logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s %(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger("geneparse-indexer")


def main():
    args = parse_args()

    # IMPUTE2
    if args.impute2:
        for fn in args.impute2:
            index_impute2(fn)

    # BGEN
    if args.bgen:
        for fn in args.bgen:
            index_bgen(fn, legacy=args.legacy)


def index_impute2(fn):
    """Indexes an IMPUTE2 file.

    Args:
        fn (str): The name of the IMPUTE2 file.

    """
    logger.info("Indexing {} (IMPUTE2)".format(fn))
    impute2_index(fn, cols=[0, 1, 2], names=["chrom", "name", "pos"], sep=" ")
    logger.info("Index generated")


def index_bgen(fn, legacy=False):
    """Indexes a BGEN file.

    Args:
        fn (str): The name of the BGEN file.

    """
    logger.info("Indexing {} (BGEN) using 'bgenix'{}".format(
        fn, " (legacy mode)" if legacy else "",
    ))
    command = ["bgenix", "-g", fn, "-index"]
    if legacy:
        command.append("-with-rowid")
    try:
        logger.info("Executing '{}'".format(" ".join(command)))
        subprocess.Popen(command).communicate()
    except FileNotFoundError:
        logger.error("Cannot find 'bgenix', impossible to index {}".format(fn))
        sys.exit(1)
    logger.info("Index generated")


def parse_args():
    """Parses the arguments and options."""
    parser = argparse.ArgumentParser(
        prog="geneparse-indexer",
        description="Genotype file indexer."
    )

    # IMPUTE2 files
    group = parser.add_argument_group("IMPUTE2 index")
    group.add_argument(
        "--impute2", metavar="IMPUTE2", type=str, nargs="+",
        help="Index an IMPUTE2 genotype file format. The file can be plain "
             "text or bgzipped.",
    )

    # BGEN files
    group = parser.add_argument_group("BGEN index")
    group.add_argument(
        "--bgen", metavar="BGEN", type=str, nargs="+",
        help="Index a BGEN genotype file. This requires 'bgenix' to be in the "
             "PATH.",
    )
    group.add_argument(
        "--legacy", action="store_true",
        help="Index the file using the '-with-rowid' option. This flag "
             "enables compatibility with SQLITE prior to version 3.8.2. See "
             "https://bitbucket.org/gavinband/bgen/wiki/bgenix for more "
             "information.",
    )

    return parser.parse_args()


if __name__ == "__main__":
    main()
