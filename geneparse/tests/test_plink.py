"""
Tests for the plink implementation.
"""

import os
import unittest
import logging

from .generic_tests import TestContainer
from .. import plink


logging.disable(logging.CRITICAL)


# TODO use the proper mechanism to get the file.
PLINK_PREFIX = os.path.abspath(os.path.join(
    os.path.dirname(__file__), "..", "..", "data", "plink", "btest"
))


class TestPlink(TestContainer, unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.reader_f = lambda x: plink.PlinkReader(PLINK_PREFIX)
