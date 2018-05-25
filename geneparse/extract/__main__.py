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

from pyplink import PyPlink


# Logging configuration
logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s %(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger("geneparse-extractor")


# Streamable output
_streamable_format = {"vcf", "csv"}


# VCF utilities
_VCF_HEADER = """##fileformat=VCFv4.3
##fileDate={date}
##source=geneparseV{version}
##INFO=<ID=AF,Number=A,Type=Float,Description="Alternative allele frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DS,Number=1,Type=Float,Description="Alternate allele dosage">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{samples}
"""
_VCF_GT_MAP = {0: "0/0", 1: "0/1", 2: "1/1"}


# Plink utilities
_PLINK_CHROM_ENCODE = {"X": "23", "Y": "24", "XY": "25", "M": "26", "MT": "26"}


def main():
    # Getting and checking the arguments and options
    args = parse_args()
    check_args(args)

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

        # The writer
        writer = None
        if args.output_format == "vcf":
            writer = vcf_writer
        elif args.output_format == "plink":
            writer = bed_writer
        elif args.output_format == "csv":
            writer = csv_writer

        # Executing the extraction
        GenoParser = parsers[args.input_format]
        with GenoParser(**parser_args) as parser:
            writer(parser=parser, keep=keep, extract=extract, args=args)

    finally:
        # Closing the extract file
        if args.extract is not None:
            args.extract.close()

        # Closing the keep file
        if args.keep is not None:
            args.keep.close()


def vcf_writer(parser, keep, extract, args):
    """Writes the data in VCF format."""
    # The output
    output = sys.stdout if args.output == "-" else open(args.output, "w")

    try:
        # Getting the samples
        samples = np.array(parser.get_samples(), dtype=str)
        k = _get_sample_select(samples=samples, keep=keep)

        # Writing the VCF header
        output.write(_VCF_HEADER.format(
            date=datetime.today().strftime("%Y%m%d"),
            version=__version__,
            samples="\t".join(samples[k]),
        ))

        # The data generator
        generator = _get_generator(parser=parser, extract=extract, keep=k,
                                   check_maf=args.maf)

        # The number of markers extracted
        nb_extracted = 0

        for data in generator:
            # Keeping only the required genotypes
            genotypes = data.genotypes

            # Computing the alternative allele frequency
            af = np.nanmean(genotypes) / 2

            print(data.variant.chrom, data.variant.pos, data.variant.name,
                  data.reference, data.coded, ".", "PASS", "AF={}".format(af),
                  "GT:DS", sep="\t", end="", file=output)

            for geno in genotypes:
                if np.isnan(geno):
                    output.write("\t./.:.")
                else:
                    rounded_geno = int(round(geno, 0))
                    output.write("\t{}:{}".format(
                        _VCF_GT_MAP[rounded_geno], geno,
                    ))

            output.write("\n")
            nb_extracted += 1

        if nb_extracted == 0:
            logger.warning("No markers matched the extract list")

    finally:
        output.close()


def csv_writer(parser, keep, extract, args):
    """Writes the data in CSV format."""
    # The output
    output = sys.stdout if args.output == "-" else open(args.output, "w")

    try:
        # Getting the samples
        samples = np.array(parser.get_samples(), dtype=str)
        k = _get_sample_select(samples=samples, keep=keep)

        # Writing the CSV header
        print("sample_id", "variant_id", "chromosome", "position", "reference",
              "coded", "dosage", "hard_call", sep=",", file=output)

        # The data generator
        generator = _get_generator(parser=parser, extract=extract, keep=k,
                                   check_maf=args.maf)

        # The number of markers extracted
        nb_extracted = 0

        for data in generator:
            # Keeping only the required genotypes
            genotypes = data.genotypes

            # The hard call mapping
            hard_call_mapping = {
                0: "{ref}/{ref}".format(ref=data.reference),
                1: "{ref}/{alt}".format(ref=data.reference, alt=data.coded),
                2: "{alt}/{alt}".format(alt=data.coded),
            }

            for sample, geno in zip(samples[k], genotypes):
                # Is the genotype missing
                is_missing = np.isnan(geno)

                # Hard coding (NaN values are empty string)
                hard_coded = None
                if is_missing:
                    geno = ""
                    hard_coded = ""
                else:
                    hard_coded = hard_call_mapping[int(round(geno, 0))]

                print(sample, data.variant.name, data.variant.chrom,
                      data.variant.pos, data.reference, data.coded,
                      geno, hard_coded, sep=",", file=output)

            nb_extracted += 1

        if nb_extracted == 0:
            logger.warning("No markers matched the extract list")

    finally:
        output.close()


def bed_writer(parser, keep, extract, args):
    """Writes BED/BIM/FAM files."""
    # The output bed and bim file
    bim_fn = args.output + ".bim"
    with open(bim_fn, "w") as bim, PyPlink(args.output, "w") as bed:
        # Getting the samples
        samples = np.array(parser.get_samples(), dtype=str)
        k = _get_sample_select(samples=samples, keep=keep)

        # Writing the FAM file
        with open(args.output + ".fam", "w") as fam:
            for sample in samples[k]:
                print(sample, sample, "0", "0", "0", "-1", sep=" ", file=fam)

        # Getting the data generator
        generator = _get_generator(parser=parser, extract=extract, keep=k,
                                   check_maf=args.maf)

        # The number of markers extracted
        nb_extracted = 0

        for data in generator:
            # Keeping only the required genotypes, changing NaN to -1 and
            # rounding to get a hard call
            genotypes = data.genotypes
            genotypes[np.isnan(genotypes)] = -1
            genotypes = np.round(genotypes, 0)

            # Writing the genotypes and the BIM file
            bed.write_genotypes(genotypes)
            print(
                _PLINK_CHROM_ENCODE.get(str(data.variant.chrom),
                                        data.variant.chrom),
                data.variant.name, "0", data.variant.pos, data.coded,
                data.reference, sep="\t", file=bim,
            )
            nb_extracted += 1

        if nb_extracted == 0:
            logger.warning("No markers matched the extract list")


def _get_sample_select(samples, keep):
    """Returns a vector of True/False to keep samples."""
    k = np.ones_like(samples, dtype=bool)
    if keep is not None:
        k = np.array([s in keep for s in samples], dtype=bool)
        if np.sum(k) == 0:
            logger.warning("No samples matched the keep list")
    return k


def _get_generator(parser, extract, keep, check_maf):
    """Generates the data (with extract markers and keep, if required."""
    if extract is not None:
        parser = Extractor(parser, names=extract)

    for data in parser.iter_genotypes():
        data.genotypes = data.genotypes[keep]

        # Checking the MAF, if required
        if check_maf:
            data.code_minor()

        yield data


def check_args(args):
    """Checks the arguments and options."""
    # Checking that only VCF can have a - (stdout) as output
    if args.output_format not in _streamable_format and args.output == "-":
        logger.error("{} format cannot be streamed to standard output"
                     "".format(args.output_format))
        sys.exit(1)

    # Checking the file extensions
    if args.output_format == "vcf" and args.output != "-":
        if not args.output.endswith(".vcf"):
            args.output += ".vcf"

    elif args.output_format == "plink":
        if args.output.endswith(".bed"):
            args.output = args.output[:-4]


def parse_args():
    """Parses the arguments and options."""
    parser = argparse.ArgumentParser(
        prog="geneparse-extractor",
        description="Genotype file extractor. This tool will extract markers "
                    "according to names or to genomic locations.",
        epilog="The parser arguments (PARSER_ARGS) are the same as the one in "
               "the API. For example, the arguments for the Plink parser is "
               "'prefix:PREFIX' (where PREFIX is the prefix of the "
               "BED/BIM/FAM files).",
    )

    # The input file format
    group = parser.add_argument_group("Input Options")
    group.add_argument(
        "-f", "--format", metavar="FORMAT", required=True, type=str,
        dest="input_format", choices=set(parsers.keys()),
        help="The input file format.",
    )

    group.add_argument(
        nargs="+", dest="parser_args", type=str, metavar="PARSER_ARGS",
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
    group.add_argument(
        "--maf", action="store_true",
        help="Check MAF and flip the allele coding if the MAF is higher "
             "than 50%%.",
    )

    # The output options
    group = parser.add_argument_group("Output Options")
    group.add_argument(
        "-o", "--output", metavar="FILE", type=str, required=True,
        help="The output file (can be '-' for STDOUT when using VCF or CSV as "
             "output format).",
    )
    group.add_argument(
        "--output-format", metavar="FORMAT", default="vcf", type=str,
        choices={"vcf", "plink", "csv"},
        help="The output file format. Note that the extension will be added "
             "if absent. Note that CSV is a long format (hence it might take "
             "more disk space).",
    )

    return parser.parse_args()


if __name__ == "__main__":
    main()
