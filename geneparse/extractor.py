"""
Extractor class.
"""


import logging

from .core import complement_alleles
from . import config

logger = logging.getLogger(__name__)


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

        # Assume i'm looking for A/C
        # I don't find it, but I find T/G with coded G
        if self.by_variants:
            for variant in self.variants:
                hits = _get_variant_silent(self.parser, variant)
                flip_strand = False
                if not hits and not variant.alleles_ambiguous():
                    # Try looking at the other strand.
                    comp_variant = variant.complementary_strand_copy()
                    hits = _get_variant_silent(self.parser, comp_variant)
                    flip_strand = True

                # Still couldn't find a variant.
                if not hits:
                    logger.warning("Could not extract variant {}."
                                   "".format(variant))

                for genotype in hits:
                    if flip_strand:
                        genotype.variant = variant
                        genotype.reference = complement_alleles(
                            genotype.reference
                        )
                        genotype.coded = complement_alleles(
                            genotype.coded
                        )

                    yield genotype

    def __getattr__(self, key):
        return getattr(self.parser, key)


def _get_variant_silent(parser, variant):
    """Gets a variant from the parser while disabling logging."""
    prev_log = config.LOG_NOT_FOUND
    config.LOG_NOT_FOUND = False
    results = parser.get_variant_genotypes(variant)
    config.LOG_NOT_FOUND = prev_log
    return results
