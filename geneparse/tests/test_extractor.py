"""Test the extractor."""

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
from ..extract.extractor import Extractor
from ..core import complement_alleles
from ..readers.plink import PlinkReader
from ..readers.dict_based import DictBasedReader
from ..testing.simulation import simulate_genotypes


logging.disable(logging.CRITICAL)


class TestExtractor(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.reader_f = lambda x: PlinkReader(PLINK_PREFIX)

    def test_extract_by_name(self):
        """Tests the extraction using variants name."""
        to_extract = {"rs785467", "rs140543381"}
        reader = self.reader_f()
        extractor = Extractor(reader, names=to_extract)

        seen = set()
        for genotype in extractor.iter_genotypes():
            name = genotype.variant.name
            truth = truth_genotypes[name]
            seen.add(name)

            self.assertEqual(genotype.variant, truth.variant)
            self.assertEqual(genotype.reference, truth.reference)
            self.assertEqual(genotype.coded, truth.coded)
            np.testing.assert_array_equal(genotype.genotypes, truth.genotypes)

        self.assertEqual(seen, to_extract)
        reader.close()

    def test_extract_missing_variant(self):
        """Tests extracting a missing variant."""
        n = 100
        ambiguous = None
        non_ambiguous = None
        while ambiguous is None or non_ambiguous is None:
            g = simulate_genotypes(0.1, n, 0.95)
            if ambiguous is None and g.variant.alleles_ambiguous():
                ambiguous = g
            elif non_ambiguous is None and not g.variant.alleles_ambiguous():
                non_ambiguous = g

        reader = DictBasedReader({}, ["s{}".format(i + 1) for i in range(n)])

        extractor = Extractor(reader, variants=[
            ambiguous.variant,
            non_ambiguous.variant,
            non_ambiguous.variant.complementary_strand_copy()
        ])

        results = list(extractor.iter_genotypes())
        self.assertTrue(results == [])

    def test_extract_by_variant(self):
        """Tests the extraction using variants object."""
        to_extract = {
            truth_variants["rs785467"],
            truth_variants["rs140543381"]
        }
        reader = self.reader_f()
        extractor = Extractor(reader, variants=to_extract)

        seen = set()
        for genotype in extractor.iter_genotypes():
            truth = truth_genotypes[genotype.variant.name]
            seen.add(genotype.variant)

            self.assertEqual(genotype.variant, truth.variant)
            self.assertEqual(genotype.reference, truth.reference)
            self.assertEqual(genotype.coded, truth.coded)
            np.testing.assert_array_equal(genotype.genotypes, truth.genotypes)

        self.assertEqual(seen, to_extract)
        reader.close()

    def test_extract_by_variant_other_strand(self):
        """Test extracting a variant with the 'wrong' strand."""
        # DictBasedReader(name_to_genotypes, samples)
        # Simulate a non-ambiguous variant.
        n = 100
        while True:
            g = simulate_genotypes(0.1, n, 0.95)
            if not g.variant.alleles_ambiguous():
                break

        reader = DictBasedReader(
            {"simulated": g}, ["s{}".format(i + 1) for i in range(n)]
        )

        # Create a complemented copy of the variant.
        complemented = g.variant.complementary_strand_copy()

        extractor = Extractor(reader, variants=[complemented])

        results = list(extractor.iter_genotypes())
        self.assertEqual(len(results), 1)
        complemented_g = results[0]

        # Make sure the coding is correct.
        np.testing.assert_array_equal(g.genotypes, complemented_g.genotypes)

        self.assertEqual(complemented_g.reference,
                         complement_alleles(g.reference))

        self.assertEqual(
            complemented_g.coded, complement_alleles(g.coded)
        )

    def test_multiple_extract(self):
        """Tests extracting twice (simulating a subgroup analysis)."""
        to_extract = {"rs785467", "rs140543381"}
        reader = self.reader_f()
        extractor = Extractor(reader, names=to_extract)

        for i in range(2):
            seen = set()
            for genotype in extractor.iter_genotypes():
                name = genotype.variant.name
                truth = truth_genotypes[name]
                seen.add(name)

                self.assertEqual(genotype.variant, truth.variant)
                self.assertEqual(genotype.reference, truth.reference)
                self.assertEqual(genotype.coded, truth.coded)
                np.testing.assert_array_equal(genotype.genotypes,
                                              truth.genotypes)

            self.assertEqual(seen, to_extract)

        reader.close()
