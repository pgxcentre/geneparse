"""Simple script to extract markers from genotype files of different format."""

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
from datetime import datetime

import numpy as np

from .. import parsers
from .. import __version__
from .extractor import Extractor


# Logging configuration
logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s %(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger("geneparse-extractor")


VCF_HEADER = """##fileformat=VCFv4.3
##fileDate={date}
##source=geneparseV{version}
##INFO=<ID=AF,Number=A,Type=Float,Description="Alternative allele frequency in the initial population">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DS,Number=1,Type=Float,Description="Alternate allele dosage">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{samples}
"""


VCF_GT_MAP = {0: "0/0", 1: "0/1", 2: "1/1"}


def main():
    args = parse_args()

    try:
        # We need to clean the parser's arguments
        parser_args = {}
        for argument in args.parser_args:
            key, value = argument.split(":")

            # Trying to infer the type
            if value.isdigit():
                # Type is integer
                value = int(value)
            else:
                try:
                    # Type is float?
                    value = float(value)
                except ValueError:
                    pass

            parser_args[key] = value

        # Markers to extract
        extract = None
        if args.extract is not None:
            extract = set(args.extract.read().splitlines())
            logger.info("Extracting {:,d} markers".format(len(extract)))

        # Samples to keep
        keep = None
        if args.keep is not None:
            # Reading the samples to keep
            keep = set(args.keep.read().splitlines())
            logger.info("Keeping {:,d} samples".format(len(keep)))

        # Executing the extraction
        GenoParser = parsers[args.input_format]
        with GenoParser(**parser_args) as parser:
            # Getting the samples
            samples = np.array(parser.get_samples(), dtype=str)
            k = np.ones_like(samples, dtype=bool)
            if keep is not None:
                k = np.array([s in keep for s in samples], dtype=bool)

            # Writing the VCF header
            sys.stdout.write(VCF_HEADER.format(
                date=datetime.today().strftime("%Y%m%d"),
                version=__version__,
                samples="\t".join(samples[k]),
            ))

            # The data generator
            generator = parser
            if extract is not None:
                generator = Extractor(parser, names=extract)

            for data in generator.iter_genotypes():
                # Keeping only the required genotypes
                genotypes = data.genotypes[k]

                # Computing the alternative allele frequency
                af = np.nanmean(data.genotypes) / 2

                print(data.variant.chrom, data.variant.pos, data.variant.name,
                      data.reference, data.coded, ".", "PASS",
                      "AF={}".format(af), "GT:DS", sep="\t", end="",
                      file=sys.stdout)

                for geno in genotypes:
                    if np.isnan(geno):
                        sys.stdout.write("\t./.:.")
                    else:
                        rounded_geno = int(round(geno, 0))
                        sys.stdout.write("\t{}:{}".format(
                            VCF_GT_MAP[rounded_geno], geno,
                        ))

                sys.stdout.write("\n")

    finally:
        # Closing the output file
        args.output.close()

        # Closing the extract file
        if args.extract is not None:
            args.extract.close()

        # Closing the keep file
        if args.keep is not None:
            args.keep.close()


def parse_args():
    """Parses the arguments and options."""
    parser = argparse.ArgumentParser(
        prog="geneparse-extractor",
        description="Genotype file extractor. This tool will extract markers "
                    "according to names or to genomic locations."
    )

    # The input file format
    group = parser.add_argument_group("Input Options")
    group.add_argument(
        "-f", "--format", metavar="FORMAT", required=True, type=str,
        dest="input_format", choices=set(parsers.keys()),
        help="The input file format.",
    )

    group.add_argument(
        nargs="+", dest="parser_args", type=str,
        help="The arguments that will be passed to the genotype parsers.",
    )

    # The extract options
    group = parser.add_argument_group("Extract Options")
    group.add_argument(
        "-e", "--extract", metavar="FILE", type=argparse.FileType("r"),
        help="The list of markers to extract (one per line, no header).",
    )
    group.add_argument(
        "-k", "--keep", metavar="FILE", type=argparse.FileType("r"),
        help="The list of samples to keep (one per line, no header).",
    )

    # The output options
    group = parser.add_argument_group("Output Options")
    group.add_argument(
        "-o", "--output", metavar="FILE", type=argparse.FileType("w"),
        required=True, help="The output file (can be '-' for STDOUT).",
    )

    return parser.parse_args()


if __name__ == "__main__":
    main()
