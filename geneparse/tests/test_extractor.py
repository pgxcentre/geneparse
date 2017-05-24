"""
"""

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


import unittest
import logging

import numpy as np

from .test_plink import PLINK_PREFIX
from .truth import genotypes as truth_genotypes
from .truth import variants as truth_variants
from ..extractor import Extractor
from ..readers.plink import PlinkReader


logging.disable(logging.CRITICAL)


class TestExtractor(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.reader_f = lambda x: PlinkReader(PLINK_PREFIX)

    def test_extract_by_name(self):
        to_extract = {"rs785467", "rs140543381"}
        reader = self.reader_f()
        extractor = Extractor(reader, names=to_extract)

        seen = set()
        for genotype in extractor:
            name = genotype.variant.name
            truth = truth_genotypes[name]
            seen.add(name)

            self.assertEqual(genotype.variant, truth.variant)
            self.assertEqual(genotype.reference, truth.reference)
            self.assertEqual(genotype.coded, truth.coded)
            np.testing.assert_array_equal(genotype.genotypes, truth.genotypes)

        self.assertEqual(seen, to_extract)
        reader.close()

    def test_extract_by_variant(self):
        to_extract = {
            truth_variants["rs785467"],
            truth_variants["rs140543381"]
        }
        reader = self.reader_f()
        extractor = Extractor(reader, variants=to_extract)

        seen = set()
        for genotype in extractor:
            truth = truth_genotypes[genotype.variant.name]
            seen.add(genotype.variant)

            self.assertEqual(genotype.variant, truth.variant)
            self.assertEqual(genotype.reference, truth.reference)
            self.assertEqual(genotype.coded, truth.coded)
            np.testing.assert_array_equal(genotype.genotypes, truth.genotypes)

        self.assertEqual(seen, to_extract)
        reader.close()
