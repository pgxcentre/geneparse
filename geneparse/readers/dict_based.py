"""
Reader implementation based on a Python dict.

This is useful when you have pickled genotypes instances, for example.

"""


from ..core import GenotypesReader, Variant
from .. import logging


class DictBasedReader(GenotypesReader):
    def __init__(self, name_to_genotypes, samples):
        self._dict = name_to_genotypes
        self.samples = samples

    def iter_genotypes(self):
        return self._dict.values()

    def iter_variants(self):
        for g in self._dict.values():
            yield g.variant

    def get_variant_by_name(self, name):
        try:
            return [self._dict[name]]
        except KeyError:
            logging.variant_name_not_found(name)
            return []

    def get_variant_genotypes(self, variant):
        out = []
        for g in self._dict.values():
            if g.variant == variant:
                out.append(g)
        return out

    def get_variants_in_region(self, chrom, start, end):
        chrom = Variant._encode_chr(chrom)
        def _in_region(v):
            return (chrom == v.chrom and
                    start <= v.pos <= end)

        return [
            g for g in self._dict.values() if _in_region(g.variant)
        ]

    def get_samples(self):
        return self.samples

    def get_number_samples(self):
        return len(self.samples)

    def get_number_variants(self):
        return len(self._dict)
