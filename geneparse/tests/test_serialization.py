"""Test the serialization of Variant and Genotypes."""

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


import pickle
import tempfile
import unittest

import numpy as np

from ..core import ImputedVariant, Variant, Genotypes


class TestSerializer(unittest.TestCase):
    def test_genotypes_serialization(self):
        v = Variant("rs1235", 1, 123456, "TC")
        genotypes = np.array([0, 1, np.nan, 2, 0, 1, np.nan, 0, 0])
        g = Genotypes(v, genotypes, "C", "T", multiallelic=False)

        with tempfile.TemporaryFile() as f:
            pickle.dump(g, f)

            f.seek(0)
            recovered = pickle.load(f)

        self.assertEqual(v, recovered.variant)
        self.assertEqual(v.name, recovered.variant.name)

        self.assertEqual(g.reference, recovered.reference)
        self.assertEqual(g.coded, recovered.coded)
        self.assertEqual(g.multiallelic, recovered.multiallelic)

        np.testing.assert_array_almost_equal(genotypes, recovered.genotypes)

    def test_genotypes_serialization_no_name(self):
        v = Variant(None, 1, 123456, "TC")
        genotypes = np.array([0, 1, np.nan, 2, 0, 1, np.nan, 0, 0])
        g = Genotypes(v, genotypes, "C", "T", multiallelic=False)

        with tempfile.TemporaryFile() as f:
            pickle.dump(g, f)

            f.seek(0)
            recovered = pickle.load(f)

        self.assertTrue(recovered.variant.name is None)

    def test_imputed_genotypes_serialization(self):
        v = ImputedVariant("rs12345", 1, 123456, "TC", 0.87651)
        genotypes = np.random.binomial(2, 0.1, size=100000)
        g = Genotypes(v, genotypes, "T", "C", multiallelic=False)

        with tempfile.TemporaryFile() as f:
            pickle.dump(g, f)

            f.seek(0)
            recovered = pickle.load(f)

        self.assertEqual(recovered.variant.quality, v.quality)
        np.testing.assert_array_almost_equal(genotypes, recovered.genotypes)
