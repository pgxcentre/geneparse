"""Abstract tests for container implementations."""

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

from pkg_resources import resource_filename

import numpy as np

from pybgen.tests.truths import truths

from ..readers import bgen
from ..core import Variant
from .generic_tests import TestContainer


BGEN_FILE = resource_filename(
    "pybgen.tests",
    os.path.join("data", "example.8bits.bgen"),
)


class TestBGEN(TestContainer, unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.reader_f = lambda x: bgen.BGENReader(BGEN_FILE)

        # Using truths from pybgen
        cls.truth = truths["dosage"]["example.8bits.truths.txt.bz2"]

        # The expected variant object
        cls.expected_variants = {
            name: Variant(
                name=name,
                chrom=int(info["variant"].chrom),
                pos=info["variant"].pos,
                alleles=[info["variant"].a1, info["variant"].a2],
            ) for name, info in cls.truth["variants"].items()
        }

    def test_get_samples(self):
        with self.reader_f() as f:
            self.assertEqual(
                list(self.truth["samples"]), f.get_samples()
            )

    def test_iter_variants(self):
        """Test that all variants are iterated over"""
        seen = set()
        with self.reader_f() as f:
            for v in f.iter_variants():
                seen.add(v.name)
                self.assertEqual(v, self.expected_variants[v.name])

        self.assertEqual(seen, set(self.expected_variants.keys()))

    def test_iter_genotypes(self):
        """Test that the genotypes are read correctly"""
        seen = set()
        with self.reader_f() as f:
            for g in f.iter_genotypes():
                # Checking the variant
                variant = g.variant
                self.assertEqual(variant, self.expected_variants[variant.name])
                seen.add(variant.name)

                # Checking the genotypes
                expected = self.truth["variants"][variant.name]
                self.assertEqual(g.reference, expected["variant"].a1)
                self.assertEqual(g.coded, expected["variant"].a2)
                np.testing.assert_array_almost_equal(
                    g.genotypes, expected["data"],
                )

        self.assertEqual(seen, set(self.expected_variants.keys()))

    @unittest.skip("Not implemented")
    def test_multiallelic_identifier(self):
        """Test that the multiallelic flag gets set when iterating"""
        pass

    def test_get_biallelic_variant(self):
        """Test simplest possible case of variant accession."""
        random_variant = random.choice(list(self.expected_variants.values()))

        v = self.expected_variants[random_variant.name]
        with self.reader_f() as f:
            # Getting the results
            results = f.get_variant_genotypes(v)
            if v.name in {"RSID_10", "RSID_100"}:
                # Those have the same location and alleles
                self.assertEqual(2, len(results))
            else:
                # The remaining variants are unique
                self.assertEqual(1, len(results))

            for g in results:
                # Checking the variant is the same
                self.assertEqual(g.variant, random_variant)

                # Checking the genotypes
                expected = self.truth["variants"][g.variant.name]
                self.assertEqual(g.reference, expected["variant"].a1)
                self.assertEqual(g.coded, expected["variant"].a2)
                np.testing.assert_array_almost_equal(
                    g.genotypes, expected["data"],
                )

    def test_get_all_biallelic_variant(self):
        """Test simplest possible case of variant accession."""
        for random_variant in self.expected_variants.values():
            v = self.expected_variants[random_variant.name]
            with self.reader_f() as f:
                results = f.get_variant_genotypes(v)
                if v.name in {"RSID_10", "RSID_100"}:
                    self.assertEqual(2, len(results))
                else:
                    self.assertEqual(1, len(results))

                for g in results:
                    # Checking the variant is the same
                    self.assertEqual(g.variant, random_variant)

                    # Checking the genotypes
                    expected = self.truth["variants"][g.variant.name]
                    self.assertEqual(g.reference, expected["variant"].a1)
                    self.assertEqual(g.coded, expected["variant"].a2)
                    np.testing.assert_array_almost_equal(
                        g.genotypes, expected["data"],
                        err_msg="Difference for {}/{}".format(
                            g.variant.name, random_variant.name,
                        ),
                    )

    def test_get_na_biallelic_variant(self):
        """Test asking for an unavailable biallelic variant."""
        v = random.choice(list(self.expected_variants.values())).copy()
        v.alleles = v._encode_alleles(["NO", "ALLELE"])
        with self.reader_f() as f:
            g = f.get_variant_genotypes(v)
            self.assertEqual([], g)

    @unittest.skip("Not implemented")
    def test_get_multiallelic_variant_by_locus(self):
        """Test getting a multiallelic variant using a locus."""
        pass

    @unittest.skip("Not implemented")
    def test_get_multiallelic_variant_by_specific(self):
        """Find a biallelic variant at a multiallelic locus."""
        pass

    @unittest.skip("Not implemented")
    def test_get_multiallelic_variant_by_unavailable(self):
        """Test asking for an unavailable biallelic variant at a multiallelic
        locus.

        """
        pass

    @unittest.skip("Not implemented")
    def test_get_multiallelic_variant_by_multiallelic(self):
        """Test asking for a multiallelic variant."""
        pass

    def test_get_variant_in_region(self):
        """Test getting a variant by region."""
        seen = set()
        with self.reader_f() as f:
            for g in f.get_variants_in_region("1", 67000, 70999):
                # Checking the variant
                self.assertEqual(
                    g.variant, self.expected_variants[g.variant.name],
                )

                # Checking the genotypes
                expected = self.truth["variants"][g.variant.name]
                self.assertEqual(g.reference, expected["variant"].a1)
                self.assertEqual(g.coded, expected["variant"].a2)
                np.testing.assert_array_almost_equal(
                    g.genotypes, expected["data"],
                )

                seen.add(g.variant.name)

        expected = set()
        for name in self.truth["variant_set"]:
            v = self.truth["variants"][name]["variant"]
            if v.pos >= 67000 and v.pos <= 70999:
                expected.add(name)

        self.assertEqual(seen, expected)

    def test_get_variant_in_empty_region(self):
        """Test getting an empty region."""
        with self.reader_f() as f:
            g = list(f.get_variants_in_region("2", 46521000, 46521005))
            self.assertEqual([], g)

    def test_get_variant_by_name(self):
        """Test getting a variant by name."""
        random_variant = random.choice(list(self.expected_variants.values()))
        with self.reader_f() as f:
            g = f.get_variant_by_name(random_variant.name)
            self.assertEqual(len(g), 1)

            g = g.pop()

            # Checking the variant
            self.assertEqual(g.variant, random_variant)

            # Checking the genotypes
            expected = self.truth["variants"][random_variant.name]
            self.assertEqual(g.reference, expected["variant"].a1)
            self.assertEqual(g.coded, expected["variant"].a2)
            np.testing.assert_array_almost_equal(
                g.genotypes, expected["data"],
            )

    def test_iter_variants_by_names(self):
        """Tests getting the variations."""
        # Getting 10 random variants
        random_variants = random.sample(
            list(self.expected_variants.values()), 10,
        )

        # Generating a map of variants
        variant_map = {v.name: v for v in random_variants}

        # Reading the file
        seen_variants = set()
        with self.reader_f() as f:
            names = [v.name for v in random_variants]
            for g in f.iter_variants_by_names(names):
                # Getting the original variant
                ori_variant = variant_map[g.variant.name]

                # Checking the variant
                self.assertEqual(g.variant, ori_variant)

                # Checking the genotypes
                expected = self.truth["variants"][g.variant.name]
                self.assertEqual(g.reference, expected["variant"].a1)
                self.assertEqual(g.coded, expected["variant"].a2)
                np.testing.assert_array_almost_equal(
                    g.genotypes, expected["data"],
                )

                seen_variants.add(g.variant.name)

        self.assertEqual(set(variant_map.keys()), seen_variants)

    def test_get_variant_by_name_invalid(self):
        """Test getting an invalid variant by name."""
        with self.reader_f() as f:
            g = f.get_variant_by_name("invalid_variant_name")
            self.assertEqual(len(g), 0)

    @unittest.skip("Not implemented")
    def test_get_multiallelic_variant_by_name(self):
        """Find a biallelic variant at a multiallelic locus by name."""
        pass


class TestBGENParallel(TestBGEN, unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.reader_f = lambda x: bgen.BGENReader(BGEN_FILE, cpus=2)

        # Using truths from pybgen
        cls.truth = truths["dosage"]["example.8bits.truths.txt.bz2"]

        # The expected variant object
        cls.expected_variants = {
            name: Variant(
                name=name,
                chrom=int(info["variant"].chrom),
                pos=info["variant"].pos,
                alleles=[info["variant"].a1, info["variant"].a2],
            ) for name, info in cls.truth["variants"].items()
        }
