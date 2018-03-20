"""Tests for the plink implementation."""

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
from ..readers import impute2


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
