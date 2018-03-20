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


import tempfile
import pickle
import unittest
import logging

import numpy as np

from .generic_tests import TestContainer
from . import truth
from ..readers import dict_based
from ..testing.simulation import simulate_genotypes


logging.disable(logging.CRITICAL)


class TestDictBasedReader(TestContainer, unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.reader_f = lambda x: dict_based.DictBasedReader(
            truth.genotypes, truth.samples,
        )

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

    def test_pickle_reader(self):
        # Create a valid pickle.
        n = 100
        samples = ["sample_{}".format(i + 1) for i in range(n)]
        genotypes = {
            "samples": samples,
            "a": simulate_genotypes(0.21, n, 0.91),
            "b": simulate_genotypes(0.15, n, 0.92),
        }

        with tempfile.NamedTemporaryFile("wb") as f:
            pickle.dump(genotypes, f)
            f.seek(0)

            reader = dict_based.PickleBasedReader(f.name)

            # The reader should be able to iterate over the genotypes.
            g_a = reader.get_variant_genotypes(genotypes["a"].variant)
            self.assertEqual(len(g_a), 1)
            g_a = g_a[0]
            self.assertEqual(g_a.variant, genotypes["a"].variant)
            np.testing.assert_array_equal(
                g_a.genotypes, genotypes["a"].genotypes
            )

            g_b = reader.get_variant_genotypes(genotypes["b"].variant)
            self.assertEqual(len(g_b), 1)
            g_b = g_b[0]
            self.assertEqual(g_a.variant, genotypes["a"].variant)
            np.testing.assert_array_equal(
                g_b.genotypes, genotypes["b"].genotypes
            )

            self.assertEqual(len(list(reader.iter_genotypes())), 2)
