"""
Tests for the dataframe implementation.
"""

import unittest
import logging

import numpy as np
import pandas as pd

from .generic_tests import TestContainer
from .. import dataframe


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
