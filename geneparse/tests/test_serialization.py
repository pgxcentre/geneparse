"""
Test the serialization of Variant and Genotypes.
"""


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
