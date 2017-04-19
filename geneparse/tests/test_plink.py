"""
Tests for the plink implementation.
"""

import os
import unittest
import logging

from pkg_resources import resource_filename

from .generic_tests import TestContainer
from .. import plink


logging.disable(logging.CRITICAL)


PLINK_PREFIX = resource_filename(
    __name__,
    os.path.join("data", "plink", "btest"),
)


class TestPlink(TestContainer, unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.reader_f = lambda x: plink.PlinkReader(PLINK_PREFIX)
