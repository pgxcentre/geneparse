"""Tests for the VCF implementation."""

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


import os
import unittest
import logging

from pkg_resources import resource_filename

from .generic_tests import TestContainer
from ..readers import vcf
from . import truth


logging.disable(logging.CRITICAL)


VCF_FILE = resource_filename(
    __name__,
    os.path.join("data", "vcf", "test.vcf.gz"),
)

VARIANT_NAME_FIX = {
    ("rs9628434", "T"): "subal_2_rs9628434",
    ("rs9628434", "A"): "subal_3_rs9628434",
}


@unittest.skipIf(not vcf.CYVCF2_AVAILABLE, "cyvcf2 is not installed")
class TestVCF(TestContainer, unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.reader_f = lambda x: vcf.VCFReader(VCF_FILE)

    def test_iter_variants(self):
        """Test that all variants are iterated over"""
        # We expect the variants in the same order as the BIM.
        expected = [
            truth.variants["rs785467"],
            truth.variants["rs146589823"],
            truth.variants["rs9628434"],
            truth.variants["rs140543381"],
        ]

        with self.reader_f() as f:
            for i, v in enumerate(f.iter_variants()):
                self.assertEqual(v, expected[i])

        self.assertEqual(i, len(expected) - 1)  # All variants were iterated.

    def test_iter_genotypes(self):
        """Test that the genotypes are read correctly"""
        with self.reader_f() as f:
            for g in f.iter_genotypes():
                variant_name = VARIANT_NAME_FIX.get(
                    (truth.variant_to_key[g.variant], g.coded),
                    truth.variant_to_key[g.variant],
                )

                expected = truth.genotypes[variant_name]
                self.assertEqual(expected, g)

    def test_multiallelic_identifier(self):
        """Test that the multiallelic flag gets set when iterating"""
        with self.reader_f() as f:
            for g in f.iter_genotypes():
                variant_name = VARIANT_NAME_FIX.get(
                    (truth.variant_to_key[g.variant], g.coded),
                    truth.variant_to_key[g.variant],
                )

                expected = truth.genotypes[variant_name]
                self.assertEqual(expected.multiallelic, g.multiallelic)

    def test_get_multiallelic_variant_by_locus(self):
        """Test getting a multiallelic variant using a locus."""
        v = truth.variants["locus_rs9628434"]
        # Coded allele to expected genotypes.
        expected = {
            "T": truth.genotypes["subal_2_rs9628434"],
            "A": truth.genotypes["subal_3_rs9628434"]
        }
        with self.reader_f() as f:
            for g in f.get_variant_genotypes(v):
                self.assertEqual(
                    expected[g.coded], g
                )
                del expected[g.coded]

        self.assertEqual(len(expected), 0)

    @unittest.skip("Not implemented")
    def test_get_variant_by_name(self):
        pass

    @unittest.skip("Not implemented")
    def test_get_multiallelic_variant_by_name(self):
        pass

    @unittest.skip("Not implemented")
    def test_get_variant_by_name_invalid(self):
        pass
