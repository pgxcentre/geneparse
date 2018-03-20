"""Utilities to facilitate unittests."""

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


import random

import numpy as np

from ..core import Variant, Genotypes


def simulate_genotypes_for_variant(variant, coded, coded_freq, n, call_rate=1):
    if variant.alleles is None or len(variant.alleles) != 2:
        raise ValueError(
            "Can only simulate genotypes for biallelic variants (with defined "
            "alleles)."
        )

    # Simulate genotypes.
    g = np.random.binomial(2, coded_freq, size=n).astype(float)

    if call_rate < 1:
        missings = np.random.binomial(1, 1 - call_rate, size=n).astype(bool)
        g[missings] = np.nan

    return Genotypes(
        variant=variant,
        genotypes=g,
        coded=coded,
        reference=(variant.alleles_set - {coded}).pop(),
        multiallelic=False
    )


def simulate_genotypes(coded_freq, n, call_rate=1):
    v = Variant(
        "simulated",
        random.randint(1, 22),
        random.randint(12345, 99999999),
        random.sample("ATGC", 2)
    )

    coded_allele = v.alleles_set.pop()

    return simulate_genotypes_for_variant(v, coded_allele, coded_freq, n,
                                          call_rate)
