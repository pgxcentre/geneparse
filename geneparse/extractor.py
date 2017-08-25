"""
Extractor class.
"""


class Extractor(object):
    def __init__(self, parser, names=None, variants=None):
        if names is None and variants is None:
            raise ValueError(
                "You need to provide either a list of variants instances or a "
                "list of names to extract."
            )

        self.parser = parser
        self.names = names
        self.variants = variants

    @property
    def by_names(self):
        return self.names is not None

    @property
    def by_variants(self):
        return self.variants is not None

    def iter_genotypes(self):
        if self.by_names:
            for name in self.names:
                for genotype in self.parser.get_variant_by_name(name):
                    yield genotype

        if self.by_variants:
            for variant in self.variants:
                for genotype in self.parser.get_variant_genotypes(variant):
                    yield genotype

    def __getattr__(self, key):
        return getattr(self.parser, key)
