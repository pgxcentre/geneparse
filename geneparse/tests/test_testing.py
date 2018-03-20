"""Tests the testing module."""

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

from .. import testing
from ..core import Variant, Genotypes


class TestTesting(unittest.TestCase):
    def test_simulate_genotypes(self):
        g = testing.simulate_genotypes(0.1, 1000, call_rate=0.99)

        self.assertTrue(isinstance(g, Genotypes))

    def test_simulate_genotypes_for_variant(self):
        v = Variant("my_variant", 1, 1234, "AT")
        g = testing.simulate_genotypes_for_variant(v, "T", 0.2, 1000)

        self.assertEqual(g.variant, v)
        self.assertEqual(g.variant.name, v.name)
        self.assertTrue(isinstance(g, Genotypes))

    def test_simulate_genotypes_for_variant_no_allele(self):
        v = Variant("my_variant", 1, 1234, None)
        with self.assertRaises(ValueError):
            testing.simulate_genotypes_for_variant(v, "T", 0.2, 1000)
