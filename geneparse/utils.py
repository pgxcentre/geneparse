"""
Utilities
"""


import numpy as np


def flip_alleles(genotypes):
    """Flip the alleles of an Genotypes instance."""
    genotypes.reference, genotypes.coded = (genotypes.coded,
                                            genotypes.reference)
    genotypes.genotypes = 2 - genotypes.genotypes
    return genotypes


def code_minor(genotypes):
    """Encode the genotypes with respect to the minor allele.

    This confirms that "reference" is the major allele and that "coded" is
    the minor allele.

    In other words, this function can be used to make sure that the genotype
    value is the number of minor alleles for an individual.

    """
    _, minor_coded = maf(genotypes)
    if not minor_coded:
        return flip_alleles(genotypes)

    return genotypes


def maf(genotypes):
    """Computes the MAF and returns a boolean indicating if the minor allele
    is currently the coded allele.

    """
    g = genotypes.genotypes
    maf = np.nansum(g) / (2 * np.sum(~np.isnan(g)))
    if maf > 0.5:
        maf = 1 - maf
        return maf, False

    return maf, True
