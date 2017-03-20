"""
Ground truth for the genotype data available in the data/ directory.

Some conventions are used:

    - Entries in variants that correspond to fully defined variants have the
      rs as the key.
    - 'locus_':  Variants starting with 'locus_' followed by an identifier
      correspond to the variant, but with the alleles unset (alleles=None).
    - 'uk_': Variants starting with 'uk_' followed by an identifier correspond
      to the variants, but without the name.
    - 'synal_n_': Variants starting with 'synal_n_' where n in an ascending
      integer identifier correspond to the variants with a synonymous notation
      for the alleles (e.g. CAGG -> C is equivalent to AGG -> '-').
    - 'subal_n_': Variants starting with 'subal_n_' for allele subsets for
      multiallelic variants.
    - 'codechr_': Variants with a numerically encoded chromosome instead of the
      usual str.

"""


from ..core import Variant, Genotypes

import numpy as np
na = np.nan


def strip_key(k):
    """Sanitize variant identifiers.

    This is necessary because variants with the same chrom, pos and alleles
    hash to the same value.

    """
    if k.startswith("uk_"):
        return k[3:]
    return k


variants = {
    "rs785467": Variant("rs785467", 1, 46521559, ["A", "T"]),
    "locus_rs785467": Variant("rs785467", 1, 46521559, None),
    "uk_rs785467": Variant(None, 1, 46521559, ["A", "T"]),

    "rs146589823": Variant("rs146589823", 2, 74601606, ["CAGG", "C"]),
    "locus_rs146589823": Variant("rs146589823", 2, 74601606, None),
    "uk_rs146589823": Variant(None, 2, 74601605, ["CAGG", "C"]),
    "synal_1_rs146589823": Variant("rs146589823", 2, 74601606, ["AGG", "-"]),

    "rs9628434": Variant("rs9628434", 22, 16615065, ["G", "A", "T"]),
    "locus_rs9628434": Variant("rs9628434", 22, 16615065, None),
    "uk_rs9628434": Variant(None, 22, 16615065, ["G", "A", "T"]),
    "subal_1_rs9628434": Variant("rs9628434", 22, 16615065, ["A", "T"]),
    "subal_2_rs9628434": Variant("rs9628434", 22, 16615065, ["G", "T"]),
    "subal_3_rs9628434": Variant("rs9628434", 22, 16615065, ["G", "A"]),

    "rs140543381": Variant("rs140543381", "X", 89932529, ["A", "T"]),
    "locus_rs140543381": Variant("rs140543381", "X", 89932529, None),
    "uk_rs140543381": Variant(None, "X", 89932529, ["A", "T"]),
    "codechr_rs140543381": Variant("rs140543381", 23, 89932529, ["A", "T"]),
}
variant_to_key = {v: strip_key(k) for k, v in variants.items()}

# Genotypes -> variant, genotype, reference, coded.
genotypes = {
    "rs785467": Genotypes(
        variants["rs785467"], np.array([0, 1, 2, 0, 0]), "A", "T", False
    ),
    "rs146589823": Genotypes(
        variants["rs146589823"], np.array([2, 1, 0, 0, 0]), "CAGG", "C", False
    ),
    "subal_2_rs9628434": Genotypes(
        variants["subal_2_rs9628434"], np.array([1, 1, na, 0, 0]), "G", "T", True
    ),
    "subal_3_rs9628434": Genotypes(
        variants["subal_3_rs9628434"], np.array([1, 0, na, 1, 0]), "G", "A", True
    ),
    "rs140543381": Genotypes(
        variants["rs140543381"], np.array([1, 2, 0, 0, 1]), "A", "T", False
    )
}

samples = ["SAMPLE{}".format(i) for i in range(1, 6)]
