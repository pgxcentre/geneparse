"""
Reader for VCF files.
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


from .core import Variant, ImputedVariant, Genotypes, GenotypesReader

from cyvcf2 import VCF
import numpy as np


class VCFReader(GenotypesReader):
    def __init__(self, filename, quality_field=None):
        self.get_vcf = lambda: VCF(filename)
        self.quality_field = quality_field

    def iter_genotypes(self):
        for v in self.get_vcf():
            alleles = set(v.REF) | set(v.ALT)

            if self.quality_field:
                variant = ImputedVariant(v.ID, v.CHROM, v.start, alleles,
                                         getattr(v, self.quality_field))
            else:
                variant = Variant(v.ID, v.CHROM, v.start, alleles)

            for coded_allele, g in self._make_genotypes(v.ALT, v.genotypes):
                yield Genotypes(variant, g, v.REF, coded_allele)

    @staticmethod
    def _make_genotypes(alleles, genotypes):
        out = []
        for allele_symbol, allele in enumerate(alleles):
            variant_geno = np.empty(len(genotypes))

            for i, g in enumerate(genotypes):
                variant_geno[i] = 0
                a1, a2, phase = g

                if a1 == allele_symbol + 1:
                    variant_geno[i] += 1

                if a2 == allele_symbol + 1:
                    variant_geno[i] += 1

            out.append((allele, variant_geno))

        return out

    def iter_variants(self):
        for v in self.get_vcf():
            yield Variant(v.ID, v.CHROM, v.start, {v.REF} | set(v.ALT))

    def get_variant_genotypes(self, variant):
        region = self.get_vcf()(
            "{}:{}-{}".format(variant.chrom, variant.pos, variant.pos)
        )
        genotypes = []

        for v in region:
            for coded_allele, g in self._make_genotypes(v.ALT, v.genotypes):
                alleles = {v.REF.upper(), coded_allele.upper()}
                match = (
                    alleles.issubset(set(variant.alleles)) or
                    variant.alleles is None
                )

                if match:
                    genotypes.append(Genotypes(g, v.REF, coded_allele))

        return genotypes

    def get_variants_in_region(self, chrom, start, end):
        region = self.get_vcf()(
            "{}:{}-{}".format(chrom, start, end)
        )
        for v in region:
            for coded_allele, g in self._make_genotypes(v.ALT, v.genotypes):
                variant = Variant(
                    v.ID, v.CHROM, v.start, [v.ref, coded_allele]
                )
                yield Genotypes(variant, g, v.REF, coded_allele)

    def get_samples(self):
        return self.get_vcf().samples

    def get_number_samples(self):
        return len(self.get_samples())

    def get_number_variants(self):
        raise NotImplementedError("Don't know how to do this using cyvcf2.")
