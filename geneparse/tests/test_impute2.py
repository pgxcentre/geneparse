"""
Tests for the plink implementation.
"""

import os
import unittest
import logging

from pkg_resources import resource_filename

from .generic_tests import TestContainer
from .. import impute2


logging.disable(logging.CRITICAL)


IMPUTE2_FN = resource_filename(
    __name__,
    os.path.join("data", "impute2", "impute2_test.impute2.gz"),
)
IMPUTE2_SAMPLE_FN = resource_filename(
    __name__,
    os.path.join("data", "impute2", "impute2_test.sample"),
)


class TestImpute2(TestContainer, unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.reader_f = lambda x: impute2.Impute2Reader(
            filename=IMPUTE2_FN,
            sample_filename=IMPUTE2_SAMPLE_FN,
        )

    @unittest.skip("Not implemented")
    def test_get_multiallelic_variant_by_specific(self):
        pass

    @unittest.skip("Not implemented")
    def test_multiallelic_identifier(self):
        pass

    @unittest.skip("Not implemented")
    def test_get_na_biallelic_variant(self):
        pass

    @unittest.skip("Not implemented")
    def test_get_multiallelic_variant_by_unavailable(self):
        pass

    @unittest.skip("Not implemented")
    def test_get_multiallelic_variant_by_multiallelic(self):
        pass

    @unittest.skip("Not implemented")
    def test_get_multiallelic_variant_by_locus(self):
        pass

    @unittest.skip("Not implemented")
    def test_get_biallelic_variant(self):
        pass
