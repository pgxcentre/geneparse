"""
Tests for the plink implementation.
"""

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
from . import truth
from ..readers import dict_based


logging.disable(logging.CRITICAL)


class TestDictBasedReader(TestContainer, unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.reader_f = lambda x: dict_based.DictBasedReader(truth.genotypes, truth.samples)

    @unittest.skip("Unsupported because keys need to be unique.")
    def test_get_multiallelic_variant_by_name(self):
        pass

    def test_iter_variants(self):
        expected = set([
            truth.variants["rs785467"],
            truth.variants["rs146589823"],
            truth.variants["subal_2_rs9628434"],
            truth.variants["subal_3_rs9628434"],
            truth.variants["rs140543381"],
        ])

        with self.reader_f() as f:
            observed = set(f.iter_variants())

        self.assertEqual(expected, observed)
