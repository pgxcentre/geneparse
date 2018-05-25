"""Tests for the utilities."""

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
import random
import unittest

import numpy as np
import pandas as pd

from pkg_resources import resource_filename

from ..testing.simulation import simulate_genotypes
from .. import utils
from ..core import Variant, Genotypes
from ..readers.plink import PlinkReader


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
        g = simulate_genotypes(0.4, 100, call_rate=1)
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


class TestUtilsLD(unittest.TestCase):
    @staticmethod
    def _read_data(prefix):
        rs1800775 = None
        others = []
        with PlinkReader(prefix) as bed:
            for data in bed.iter_genotypes():
                if data.variant.name == "rs1800775":
                    rs1800775 = data
                else:
                    others.append(data)

        return rs1800775, others

    @staticmethod
    def _read_ld(fn):
        return pd.read_csv(fn, delim_whitespace=True).set_index(
            "SNP_B", verify_integrity=True,
        ).rename(columns={"R2": "expected_r2"}).expected_r2

    def setUp(self):
        # The prefix of the two datasets
        prefix = resource_filename(
            __name__,
            os.path.join("data", "ld", "common_extracted_1kg"),
        )

        # Retrieving the two datasets
        self.rs1800775, self.others = self._read_data(prefix)
        self.rs1800775_missing, self.others_missing = self._read_data(
            prefix + ".missing",
        )

        # The prefix of the two LD files
        prefix = resource_filename(
            __name__, os.path.join("data", "ld", "plink_rs1800775_pairs"),
        )

        # Retrieving the two LD files
        self.exp_ld = self._read_ld(prefix + ".ld")
        self.exp_ld_missing = self._read_ld(prefix + ".missing.ld")

    def test_ld_computation(self):
        # Compute the LD
        ld_vector = utils.compute_ld(self.rs1800775, self.others, r2=True)

        # Checking the size is appropriate
        self.assertEqual(ld_vector.shape[0], len(self.others))

        # Checking we have a value for each of the markers
        self.assertEqual(
            set(ld_vector.index), {g.variant.name for g in self.others},
        )

        # Merging with the expected values
        df = pd.DataFrame(
            {"r2": ld_vector, "expected_r2": self.exp_ld},
        )
        self.assertEqual(df.shape[0] - 1, ld_vector.shape[0])

        # Computing the square error (should be close to 0)
        squared_error = np.nanmean((df.r2 - df.expected_r2) ** 2)
        self.assertAlmostEqual(squared_error, 0)

    def test_ld_computation_with_na_values(self):
        # Computing the LD
        ld_vector = utils.compute_ld(
            self.rs1800775_missing, self.others_missing, r2=True,
        )

        # Checking the size is appropriate
        self.assertEqual(ld_vector.shape[0], len(self.others))

        # Checking we have a value for each of the markers
        self.assertEqual(
            set(ld_vector.index), {g.variant.name for g in self.others},
        )

        # Merging with the expected values
        df = pd.DataFrame(
            {"r2": ld_vector, "expected_r2": self.exp_ld_missing},
        )
        self.assertEqual(df.shape[0] - 1, ld_vector.shape[0])

        # Computing the square error (should be close to 0)
        squared_error = np.nanmean((df.r2 - df.expected_r2) ** 2)
        self.assertAlmostEqual(squared_error, 0, places=5)
