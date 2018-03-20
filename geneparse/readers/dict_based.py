"""Reader implementation based on a Python dict.

This is useful when you have pickled genotypes instances, for example.

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


import pickle


from ..core import GenotypesReader, Variant
from .. import logging


class DictBasedReader(GenotypesReader):
    def __init__(self, name_to_genotypes, samples):
        self._dict = name_to_genotypes
        self.samples = samples

    def iter_genotypes(self):
        return (i.copy() for i in self._dict.values())

    def iter_variants(self):
        for g in self._dict.values():
            yield g.variant.copy()

    def get_variant_by_name(self, name):
        try:
            return [self._dict[name].copy()]
        except KeyError:
            logging.variant_name_not_found(name)
            return []

    def get_variant_genotypes(self, variant):
        out = []
        for g in self._dict.values():
            if g.variant == variant:
                out.append(g.copy())
        return out

    def get_variants_in_region(self, chrom, start, end):
        chrom = Variant._encode_chr(chrom)

        def _in_region(v):
            return (chrom == v.chrom and
                    start <= v.pos <= end)

        return [
            g.copy() for g in self._dict.values() if _in_region(g.variant)
        ]

    def get_samples(self):
        return self.samples

    def get_number_samples(self):
        return len(self.samples)

    def get_number_variants(self):
        return len(self._dict)


class PickleBasedReader(DictBasedReader):
    def __init__(self, filename):
        with open(filename, "rb") as f:
            data = pickle.load(f)

        samples = data.pop("samples")

        geno = {}
        i = 0
        for v in data.values():
            if v.variant.name:
                name = v.variant.name
                if name in geno:
                    name += "_2"
                    logging.warning(
                        "Renaming variant {variant} to {variant}_2 because "
                        "there are duplicate names.".format(variant=v.variant)
                    )
            else:
                name = "Variant{}".format(i)
                i += 1

            geno[name] = v

        super().__init__(geno, samples)
