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


# TODO: Add tests for actual dosage value (not just 100% probability)
class TestImpute2(TestContainer, unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.reader_f = lambda x: impute2.Impute2Reader(
            filename=IMPUTE2_FN,
            sample_filename=IMPUTE2_SAMPLE_FN,
        )
