"""Tests for the core data structures."""

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

import numpy as np

from ..core import Genotypes, Variant


class TestGenotypes(unittest.TestCase):
    def setUp(self):
        # Creating a genotype object
        self.g = Genotypes(
            variant=Variant("marker_1", 1, 1234, ["ACGCT", "C"]),
            genotypes=np.array([2, 2, 1, 2, 0, 1, 1, 1, np.nan, 0]),
            reference="ACGCT",
            coded="C",
            multiallelic=False,
        )

    def test_flip_coded(self):
        """Tests flipping the coding."""
        # Flipping the encoding
        self.g.flip_coded()

        # Checking the reference and coded allele did change
        self.assertEqual(self.g.reference, "C")
        self.assertEqual(self.g.coded, "ACGCT")

        # Checking the genotypes changed
        np.testing.assert_array_equal(
            np.array([0, 0, 1, 0, 2, 1, 1, 1, np.nan, 2]),
            self.g.genotypes,
        )

        # The variant and multiallelic boolean should be unchanged
        self.assertEqual("marker_1", self.g.variant.name)
        self.assertEqual(1, self.g.variant.chrom)
        self.assertEqual(1234, self.g.variant.pos)
        self.assertEqual(("ACGCT", "C"), self.g.variant.alleles)
        self.assertFalse(self.g.multiallelic)

    def test_flip_strand(self):
        """Tests flipping the strand."""
        # Flipping the strand
        self.g.flip_strand()

        # Checking the reference and coded allele did change
        self.assertEqual(self.g.reference, "AGCGT")
        self.assertEqual(self.g.coded, "G")

        # Checking the genotypes did not changed
        np.testing.assert_array_equal(
            np.array([2, 2, 1, 2, 0, 1, 1, 1, np.nan, 0]),
            self.g.genotypes,
        )

        # The variant's allele should have changed
        self.assertEqual(("AGCGT", "G"), self.g.variant.alleles)

        # The rest of the information should be the same
        self.assertEqual("marker_1", self.g.variant.name)
        self.assertEqual(1, self.g.variant.chrom)
        self.assertEqual(1234, self.g.variant.pos)
        self.assertFalse(self.g.multiallelic)

    def test_maf(self):
        """Tests getting the MAF."""
        # Flipping the strand
        self.assertAlmostEqual(self.g.maf(), 8 / 18)

    def test_coded_freq(self):
        """Tests getting the coded allele frequency."""
        self.assertAlmostEqual(self.g.coded_freq(), 10 / 18)

    def test_code_minor_has_effect(self):
        """Tests coding the genotypes to the minor allele (need to flip)."""
        # Flipping the encoding
        self.g.code_minor()

        # Checking the reference and coded allele did change
        self.assertEqual(self.g.reference, "C")
        self.assertEqual(self.g.coded, "ACGCT")

        # Checking the genotypes changed
        np.testing.assert_array_equal(
            np.array([0, 0, 1, 0, 2, 1, 1, 1, np.nan, 2]),
            self.g.genotypes,
        )

        # The variant and multiallelic boolean should be unchanged
        self.assertEqual("marker_1", self.g.variant.name)
        self.assertEqual(1, self.g.variant.chrom)
        self.assertEqual(1234, self.g.variant.pos)
        self.assertEqual(("ACGCT", "C"), self.g.variant.alleles)
        self.assertFalse(self.g.multiallelic)

    def test_code_minor_has_no_effect(self):
        """Tests coding the genotypes to the minor allele (no need to flip)."""
        other_g = Genotypes(
            variant=Variant("marker_1", 1, 1234, ["ACGCT", "C"]),
            genotypes=np.array([0, 0, 1, 0, 2, 1, 1, 1, np.nan, 2]),
            reference="ACGCT",
            coded="C",
            multiallelic=False,
        )
        other_g.code_minor()

        # Checking the reference and coded allele did not change
        self.assertEqual(other_g.reference, "ACGCT")
        self.assertEqual(other_g.coded, "C")

        # Checking the genotypes changed
        np.testing.assert_array_equal(
            np.array([0, 0, 1, 0, 2, 1, 1, 1, np.nan, 2]),
            other_g.genotypes,
        )

        # The variant and multiallelic boolean should be unchanged
        self.assertEqual("marker_1", other_g.variant.name)
        self.assertEqual(1, other_g.variant.chrom)
        self.assertEqual(1234, other_g.variant.pos)
        self.assertEqual(("ACGCT", "C"), other_g.variant.alleles)
        self.assertFalse(other_g.multiallelic)
