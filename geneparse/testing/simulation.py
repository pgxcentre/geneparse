"""
Utilities to facilitate unittests.
"""

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
