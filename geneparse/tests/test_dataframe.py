"""Tests for the dataframe implementation."""

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
import pandas as pd

from .generic_tests import TestContainer
from ..readers import dataframe


logging.disable(logging.CRITICAL)


class TestDataFrame(TestContainer, unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        genotypes = pd.DataFrame(
            {"rs785467": [0, 1, 2, 0, 0],
             "rs146589823": [2, 1, 0, 0, 0],
             "rs9628434": [1, 0, np.nan, 1, 0],
             "rs140543381": [1, 2, 0, 0, 1]},
            index=["SAMPLE{}".format(_+1) for _ in range(5)],
        )

        mapping_info = pd.DataFrame(
            {"chrom": ["1", "2", "22", "X"],
             "pos": [46521559, 74601606, 16615065, 89932529],
             "a1": ["T", "C", "A", "T"],
             "a2": ["A", "CAGG", "G", "A"]},
            index=["rs785467", "rs146589823", "rs9628434", "rs140543381"],
        )

        cls.reader_f = lambda x: dataframe.DataFrameReader(
            dataframe=genotypes,
            map_info=mapping_info,
        )

    @unittest.skip("Not implemented")
    def test_iter_variants(self):
        """Test that all variants are iterated over"""
        pass

    @unittest.skip("Not implemented")
    def test_multiallelic_identifier(self):
        """Test that the multiallelic flag gets set when iterating"""
        pass

    @unittest.skip("Not implemented")
    def test_get_biallelic_variant(self):
        """Test simplest possible case of variant accession."""
        pass

    @unittest.skip("Not implemented")
    def test_get_na_biallelic_variant(self):
        """Test asking for an unavailable biallelic variant."""
        pass

    @unittest.skip("Not implemented")
    def test_get_multiallelic_variant_by_locus(self):
        """Test getting a multiallelic variant using a locus."""
        pass

    @unittest.skip("Not implemented")
    def test_get_multiallelic_variant_by_specific(self):
        """Find a biallelic variant at a multiallelic locus."""
        pass

    @unittest.skip("Not implemented")
    def test_get_multiallelic_variant_by_unavailable(self):
        """Test asking for an unavailable biallelic variant at a multiallelic
        locus.

        """
        pass

    @unittest.skip("Not implemented")
    def test_get_multiallelic_variant_by_multiallelic(self):
        """Test asking for a multiallelic variant."""
        pass

    @unittest.skip("Not implemented")
    def test_get_variant_in_region(self):
        """Test getting a variant by region."""
        pass

    @unittest.skip("Not implemented")
    def test_get_variant_in_empty_region(self):
        """Test getting an empty region."""
        pass

    @unittest.skip("Not implemented")
    def test_get_multiallelic_variant_by_name(self):
        """Find a biallelic variant at a multiallelic locus by name."""
        pass
