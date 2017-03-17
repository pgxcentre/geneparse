"""
Abstract tests for container implementations.
"""

from . import truth


class TestContainer(object):
    @classmethod
    def setUpClass(cls):
        raise NotImplementedError(
            "Subclasses should define a callable that returns a instantiated "
            "GenotypesReader as '.reader_f'."
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
