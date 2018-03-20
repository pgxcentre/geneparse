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


from . import truth


class TestContainer(object):
    @classmethod
    def setUpClass(cls):
        raise NotImplementedError(
            "Subclasses should define a callable that returns a instantiated "
            "GenotypesReader as '.reader_f'."
        )

    def test_get_samples(self):
        with self.reader_f() as f:
            self.assertEqual(
                truth.samples, f.get_samples()
            )

    def test_iter_variants(self):
        """Test that all variants are iterated over"""
        # We expect the variants in the same order as the BIM.
        expected = [
            truth.variants["rs785467"],
            truth.variants["rs146589823"],
            truth.variants["subal_3_rs9628434"],
            truth.variants["subal_2_rs9628434"],
            truth.variants["rs140543381"],
        ]

        with self.reader_f() as f:
            for i, v in enumerate(f.iter_variants()):
                self.assertEqual(v, expected[i])

        self.assertEqual(i, len(expected) - 1)  # All variants were iterated.

    def test_iter_genotypes(self):
        """Test that the genotypes are read correctly"""
        with self.reader_f() as f:
            for g in f.iter_genotypes():
                expected = truth.genotypes[truth.variant_to_key[g.variant]]
                self.assertEqual(expected, g)

    def test_multiallelic_identifier(self):
        """Test that the multiallelic flag gets set when iterating"""
        with self.reader_f() as f:
            for g in f.iter_genotypes():
                expected = truth.genotypes[truth.variant_to_key[g.variant]]
                self.assertEqual(expected.multiallelic, g.multiallelic)

    def test_get_biallelic_variant(self):
        """Test simplest possible case of variant accession."""
        v = truth.variants["rs785467"]
        expected = truth.genotypes["rs785467"]
        with self.reader_f() as f:
            g = f.get_variant_genotypes(v)[0]
            self.assertEqual(expected, g)

    def test_get_na_biallelic_variant(self):
        """Test asking for an unavailable biallelic variant."""
        v = truth.variants["rs785467"].copy()
        v.alleles = v._encode_alleles(["A", "G"])
        with self.reader_f() as f:
            g = f.get_variant_genotypes(v)
            self.assertEqual([], g)

    def test_get_multiallelic_variant_by_locus(self):
        """Test getting a multiallelic variant using a locus."""
        v = truth.variants["locus_rs9628434"]
        # Coded allele to expected genotypes.
        expected = {
            "T": truth.genotypes["subal_2_rs9628434"],
            "A": truth.genotypes["subal_3_rs9628434"]
        }
        with self.reader_f() as f:
            for g in f.get_variant_genotypes(v):
                self.assertEqual(
                    expected[g.coded], g
                )
                del expected[g.coded]

        self.assertEqual(len(expected), 0)

    def test_get_multiallelic_variant_by_specific(self):
        """Find a biallelic variant at a multiallelic locus."""
        # subal 2 is the G/T variant.
        v = truth.variants["subal_2_rs9628434"]
        expected = truth.genotypes["subal_2_rs9628434"]
        with self.reader_f() as f:
            g = f.get_variant_genotypes(v)
            self.assertEqual(len(g), 1)

            g = g[0]
            self.assertEqual(expected, g)

    def test_get_multiallelic_variant_by_unavailable(self):
        """Test asking for an unavailable biallelic variant at a multiallelic
        locus.

        """
        v = truth.variants["subal_1_rs9628434"]
        with self.reader_f() as f:
            self.assertEqual([], f.get_variant_genotypes(v))

    def test_get_multiallelic_variant_by_multiallelic(self):
        """Test asking for a multiallelic variant."""
        v = truth.variants["rs9628434"]

        expected = {
            "T": truth.genotypes["subal_2_rs9628434"],
            "A": truth.genotypes["subal_3_rs9628434"],
        }
        with self.reader_f() as f:
            for g in f.get_variant_genotypes(v):
                # Find the right expected genotypes.
                self.assertEqual(
                    expected[g.coded], g
                )
                del expected[g.coded]

        self.assertEqual(len(expected), 0)

    def test_get_variant_in_region(self):
        """Test getting a variant by region."""
        expected = truth.genotypes["rs785467"]
        with self.reader_f() as f:
            genotypes = list(f.get_variants_in_region("1", 46521558, 46521560))
            self.assertEqual(len(genotypes), 1)

            genotypes = genotypes[0]
            self.assertEqual(expected, genotypes)

    def test_get_variant_in_empty_region(self):
        """Test getting an empty region."""
        with self.reader_f() as f:
            g = list(f.get_variants_in_region("1", 46521000, 46521005))
            self.assertEqual([], g)

    def test_get_variant_by_name(self):
        """Test getting a variant by name."""
        with self.reader_f() as f:
            g = f.get_variant_by_name("rs146589823")
            self.assertEqual(len(g), 1)

            g = g[0]
            self.assertEqual(g, truth.genotypes["rs146589823"])

    def test_get_variant_by_name_invalid(self):
        """Test getting an invalid variant by name."""
        with self.reader_f() as f:
            g = f.get_variant_by_name("invalid_variant_name")
            self.assertEqual(len(g), 0)

    def test_get_multiallelic_variant_by_name(self):
        """Find a biallelic variant at a multiallelic locus by name."""
        # Inverted, because inverted in the BIM file
        expected_1 = truth.genotypes["subal_3_rs9628434"]
        expected_2 = truth.genotypes["subal_2_rs9628434"]
        with self.reader_f() as f:
            r = f.get_variant_by_name("rs9628434")
            self.assertEqual(len(r), 2)

            for g, e in zip(r, (expected_1, expected_2)):
                self.assertEqual(e, g)
