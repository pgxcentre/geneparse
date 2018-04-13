"""Extractor class."""

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


import logging

from ..core import complement_alleles
from .. import config

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
            yield from self.parser.iter_variants_by_names(self.names)

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
