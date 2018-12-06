"""Reader for VCF files."""

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


from ..core import Variant, ImputedVariant, Genotypes, GenotypesReader

import numpy as np

try:
    from cyvcf2 import VCF
    CYVCF2_AVAILABLE = True
except ImportError:
    CYVCF2_AVAILABLE = False


class VCFReader(GenotypesReader):
    def __init__(self, filename, quality_field=None):
        if not CYVCF2_AVAILABLE:
            raise RuntimeError(
                "cyvcf2 is not installed. Install cyvcf2 (from source or "
                "using `pip install cyvcf2`) to use the VCF reader."
            )

        self.get_vcf = lambda: VCF(filename)
        self.quality_field = quality_field

    def __repr__(self):
        # Impossible to know the number of variants without reading the
        # complete file... so we only show the number of samples...
        return "<{} {:,d} samples>".format(
            self.__class__.__name__,
            self.get_number_samples(),
        )

    def iter_genotypes(self):
        """Iterates on available markers.

        Returns:
            Genotypes instances.

        """
        for v in self.get_vcf():
            alleles = {v.REF} | set(v.ALT)

            if self.quality_field:
                variant = ImputedVariant(v.ID, v.CHROM, v.POS, alleles,
                                         getattr(v, self.quality_field))
            else:
                variant = Variant(v.ID, v.CHROM, v.POS, alleles)

            for coded_allele, g in self._make_genotypes(v.ALT, v.genotypes):
                yield Genotypes(variant, g, v.REF, coded_allele,
                                multiallelic=len(v.ALT) > 1)

    @staticmethod
    def _make_genotypes(alleles, genotypes):
        out = []
        for allele_symbol, allele in enumerate(alleles):
            variant_geno = np.empty(len(genotypes))

            for i, g in enumerate(genotypes):
                variant_geno[i] = 0
                a1, a2, phase = g

                if a1 == -1 or a2 == -1:
                    variant_geno[i] = np.nan
                    continue

                if a1 == allele_symbol + 1:
                    variant_geno[i] += 1

                if a2 == allele_symbol + 1:
                    variant_geno[i] += 1

            out.append((allele, variant_geno))

        return out

    def iter_variants(self):
        """Iterate over marker information."""
        for v in self.get_vcf():
            yield Variant(v.ID, v.CHROM, v.POS, {v.REF} | set(v.ALT))

    def get_variant_genotypes(self, variant):
        region = self.get_vcf()(
            "{}:{}-{}".format(variant.chrom, variant.pos, variant.pos)
        )
        genotypes = []

        for v in region:
            for coded_allele, g in self._make_genotypes(v.ALT, v.genotypes):
                alleles = {v.REF, coded_allele}
                match = (
                    variant.alleles is None or
                    alleles.issubset(set(variant.alleles))
                )

                if match:
                    genotypes.append(Genotypes(
                        Variant(v.ID, v.CHROM, v.POS, alleles),
                        g, v.REF, coded_allele,
                        multiallelic=len(v.ALT) > 1,
                    ))

        return genotypes

    def get_variants_in_region(self, chrom, start, end):
        """Iterate over variants in a region."""
        region = self.get_vcf()(
            "{}:{}-{}".format(chrom, start, end)
        )
        for v in region:
            for coded_allele, g in self._make_genotypes(v.ALT, v.genotypes):
                variant = Variant(
                    v.ID, v.CHROM, v.POS, [v.REF, coded_allele]
                )
                yield Genotypes(variant, g, v.REF, coded_allele,
                                multiallelic=len(v.ALT) > 1)

    def get_samples(self):
        return self.get_vcf().samples

    def get_number_samples(self):
        return len(self.get_samples())

    def get_number_variants(self):
        raise NotImplementedError("Don't know how to do this using cyvcf2.")
