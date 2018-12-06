"""A module to parse genetics file formats."""

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


import re

from .readers import plink, impute2, dataframe, bgen, dict_based, vcf
from .core import (Genotypes, Variant, ImputedVariant, SplitChromosomeReader,
                   Chromosome)
from .extract.extractor import Extractor

try:
    from .version import geneparse_version as __version__
except ImportError:
    __version__ = None


__author__ = "Marc-Andre Legault"
__copyright__ = "Copyright 2014, Beaulieu-Saucier Pharmacogenomics Centre"
__credits__ = ["Louis-Philippe Lemieux Perreault", "Marc-Andre Legault"]
__license__ = "MIT"
__maintainer__ = "Louis-Philippe Lemieux Perreault"
__email__ = "louis-philippe.lemieux.perreault@statgen.org"
__status__ = "Development"


# TODO:
# 1. Warn and show last exception if no reader correctly initialized.
# 2. Could also make it async to load faster.
class _SplitChromosomeReaderFactory(object):
    def __init__(self, reader_class):
        self.reader_class = reader_class

    def __call__(self, pattern, *args, **kwargs):
        if "{chrom}" not in pattern:
            raise ValueError("Expected '{chrom}' as a placeholder in the "
                             "pattern.")

        # Explode the path for every possible chromosome.
        chrom_to_reader = {}
        last_exception = None
        for chrom in list(range(1, 23)) + ["X", "Y", "XY", "MT"]:
            chrom = str(chrom)
            cur = re.sub("{chrom}", chrom, pattern)
            try:
                # Instantiate the reader.
                chrom_to_reader[chrom] = self.reader_class(
                    cur, *args, **kwargs
                )
            except Exception as e:
                last_exception = e

        if len(chrom_to_reader) == 0:
            raise ValueError(
                "Could not initialize any genotype reader for chromosomes 1 "
                "to 22 or X, Y, XY, MT.\nLast exception was: {}."
                "".format(last_exception)
            )

        return SplitChromosomeReader(chrom_to_reader)


parsers = {
    "plink": plink.PlinkReader,
    "bgen": bgen.BGENReader,
    "vcf": vcf.VCFReader,
    "chrom-split-plink": _SplitChromosomeReaderFactory(plink.PlinkReader),
    "impute2": impute2.Impute2Reader,
    "chrom-split-impute2": _SplitChromosomeReaderFactory(
        impute2.Impute2Reader
    ),
    "chrom-split-bgen": _SplitChromosomeReaderFactory(bgen.BGENReader),
    "dataframe": dataframe.DataFrameReader,
    "dict-based": dict_based.DictBasedReader,
    "pickle": dict_based.PickleBasedReader,
}
