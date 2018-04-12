"""Tests the Variant object."""

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

from ..exceptions import InvalidChromosome
from ..core import Variant, Chromosome, UNKNOWN_CHROMOSOME


logging.disable(logging.CRITICAL)


class TestVariant(unittest.TestCase):
    def test_valid_variants(self):
        """Tests the creation of valid variants."""
        var = Variant("rs1234", 3, 1234, "AG")
        self.assertEqual("rs1234", var.name)
        self.assertEqual(Chromosome(3), var.chrom)
        self.assertEqual(1234, var.pos)
        self.assertEqual({"A", "G"}, var.alleles_set)

        var = Variant("rs12345", "3", 23456, "TGA")
        self.assertEqual("rs12345", var.name)
        self.assertEqual(Chromosome(3), var.chrom)
        self.assertEqual(23456, var.pos)
        self.assertEqual({"T", "G", "A"}, var.alleles_set)

        var = Variant("rs12345", "3", 23456, None)
        self.assertEqual("rs12345", var.name)
        self.assertEqual(Chromosome(3), var.chrom)
        self.assertEqual(23456, var.pos)
        self.assertTrue(var.alleles_set is None)

        var = Variant("rs12345", "chr3", 23456, "TGA")
        self.assertEqual("rs12345", var.name)
        self.assertEqual(Chromosome(3), var.chrom)
        self.assertEqual(23456, var.pos)
        self.assertEqual({"T", "G", "A"}, var.alleles_set)

        var = Variant("rs12346", Chromosome("Unknown"), 123456, "TA")
        self.assertEqual("rs12346", var.name)
        self.assertEqual(Chromosome("Unknown"), var.chrom)
        self.assertEqual(123456, var.pos)
        self.assertEqual({"T", "A"}, var.alleles_set)

        var = Variant("rs12346", None, 123456, "TA")
        self.assertEqual("rs12346", var.name)
        self.assertTrue(var.chrom is UNKNOWN_CHROMOSOME)
        self.assertEqual(123456, var.pos)
        self.assertEqual({"T", "A"}, var.alleles_set)

    def test_invalid_variant(self):
        """Tests the creation of invalid variants (e.g. invalid chrom)."""
        with self.assertRaises(InvalidChromosome) as cm:
            Variant("rs1234", "Unknown", 123435, "CG")
        self.assertEqual("UNKNOWN", str(cm.exception))

    def test_variant_comparison(self):
        """Tests comparison between variants."""
        self.assertEqual(
            Variant("rs1234", 3, 1234, "AG"),
            Variant("rs1234", 3, 1234, "AG"),
        )

        self.assertEqual(
            Variant("rs1234", 3, 1234, None),
            Variant("rs1234", 3, 1234, "AG"),
        )

        self.assertNotEqual(
            Variant("rs1234", 3, 1234, "AG"),
            Variant("rs1234", 3, 1234, "TG"),
        )

        self.assertNotEqual(
            Variant("rs1234", None, 1234, "AG"),
            Variant("rs1234", None, 1234, "AG"),
        )

    def test_variant_complement(self):
        v = Variant("rs1234", 1, 1234, ["CGTA", "ATTGCGC"])
        v.complement_alleles()
        expected = set(["GCGCAAT", "TACG"])
        self.assertEqual(v.alleles_set, expected)

        v = Variant(None, 1, 1234, ["A", "G", "T", "Z"])
        v.complement_alleles()
        expected = set(["T", "C", "A", "Z"])
        self.assertEqual(v.alleles_set, expected)
