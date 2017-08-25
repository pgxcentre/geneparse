"""
Tests for the utilities.
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

import random
import unittest

import numpy as np
import pandas as pd

from ..testing.simulation import simulate_genotypes
from .. import utils
from ..core import Variant, Genotypes


class TestUtils(unittest.TestCase):
    def test_genotype_to_df(self):
        g = simulate_genotypes(0.2, 100, call_rate=1)
        samples = ["sample{}".format(i + 1) for i in range(100)]

        # Convert to dataframe.
        df = utils.genotype_to_df(g, samples)

        self.assertEqual(df.columns, [g.variant.name])
        np.testing.assert_almost_equal(g.genotypes, df[g.variant.name])
        self.assertEqual(samples, list(df.index))

    def test_genotype_to_df_no_name(self):
        g = simulate_genotypes(0.2, 100, call_rate=1)
        samples = ["sample{}".format(i + 1) for i in range(100)]
        g.variant.name = None

        df = utils.genotype_to_df(g, samples)

        self.assertEqual(df.columns, ["genotypes"])

    def test_genotype_to_df_alleles(self):
        g = simulate_genotypes(0.2, 100, call_rate=1)
        ref_allele = g.reference
        alt_allele = g.coded

        samples = ["sample{}".format(i + 1) for i in range(100)]

        homo_ref = [
            sample for sample, geno in zip(samples, g.genotypes) if geno == 0
        ]
        hetero = [
            sample for sample, geno in zip(samples, g.genotypes) if geno == 1
        ]
        homo_coded = [
            sample for sample, geno in zip(samples, g.genotypes) if geno == 2
        ]

        df = utils.genotype_to_df(g, samples, as_string=True)

        rr = "{}/{}".format(ref_allele, ref_allele)
        rc = "{}/{}".format(ref_allele, alt_allele)
        cc = "{}/{}".format(alt_allele, alt_allele)

        homo_maj_alleles = set(df.loc[homo_ref, g.variant.name])
        self.assertEqual(len(homo_maj_alleles), 1)
        self.assertEqual(list(homo_maj_alleles)[0], rr)

        hetero_alleles = set(df.loc[hetero, g.variant.name])
        self.assertEqual(len(hetero_alleles), 1)
        self.assertEqual(list(hetero_alleles)[0], rc)

        homo_min_alleles = set(df.loc[homo_coded, g.variant.name])
        self.assertEqual(len(homo_min_alleles), 1)
        self.assertEqual(list(homo_min_alleles)[0], cc)

    def test_genotype_to_df_alleles_dosage(self):
        v = Variant(
            "dummy", random.randint(1, 22), random.randint(12345, 99999999),
            random.sample("ATGC", 2)
        )

        coded_allele = v.alleles_set.pop()

        geno = np.array([0.1, 0.8, 1.1, np.nan, 1.6, 2.0])
        samples = ["sample{}".format(i + 1) for i in range(len(geno))]

        geno = Genotypes(
            variant=v,
            genotypes=geno,
            coded=coded_allele,
            reference=(v.alleles_set - {coded_allele}).pop(),
            multiallelic=False
        )

        r = geno.reference
        c = geno.coded

        rr = "{}/{}".format(r, r)
        rc = "{}/{}".format(r, c)
        cc = "{}/{}".format(c, c)

        expected = pd.DataFrame(
            [rr, rc, rc, np.nan, cc, cc],
            index=samples,
            columns=[v.name]
        )

        observed = utils.genotype_to_df(geno, samples, as_string=True)

        self.assertEqual(expected.columns, observed.columns)
        self.assertTrue(expected.equals(observed))
