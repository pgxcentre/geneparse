"""
Compare the genotype calls between two files.
"""

import numpy as np

from ..utils import flip_alleles


def compare(reader1, reader2):
    samples1 = reader1.get_samples()
    samples2 = reader2.get_samples()

    common_samples = list(set(samples1) & set(samples2))

    sample_to_index_1 = {sample: i for i, sample in enumerate(samples1)}
    sample_to_index_2 = {sample: i for i, sample in enumerate(samples2)}

    idx1 = np.array((sample_to_index_1[s] for s in common_samples))
    idx2 = np.array((sample_to_index_2[s] for s in common_samples))

    f = open("compare_calls.csv", "w")
    f.write("name,chrom,pos,a1,a2,n_samples,n_match,n_mismatch,n_missing_1,"
            "n_missing_2\n")

    for geno1 in reader1.iter_genotypes():

        # Get the variant in reader2.
        geno2 = reader2.get_variant_genotypes(geno1.variant)

        if len(geno2) == 0:
            # Could not find variant in second container.
            continue

        elif len(geno2) == 1:
            # Single hit, we compare.
            geno2 = geno2.pop()

        else:
            # Select the variant with alleles matching if available.
            match = False
            for g in geno2:
                alleles = {g.reference, g.coded}
                if geno1.reference in alleles and geno1.coded in alleles:
                    geno2 = g
                    match = True
                    break

            if not match:
                continue

        counts = count_match_mismatch(geno1, idx1, geno2, idx2)

        f.write(",".join([str(i) for i in [
            "{} / {}".format(geno1.variant.name, geno2.variant.name),
            geno1.variant.chrom,
            geno1.variant.pos,
            geno1.reference,
            geno1.coded,
            counts["n"],
            counts["match"],
            counts["mismatch"],
            counts["missing_1"],
            counts["missing_2"],
        ]]) + "\n")

    f.close()


def count_match_mismatch(geno1, idx1, geno2, idx2):
    if geno1.reference == geno2.reference and geno1.coded == geno2.coded:
        # Exact same variant.
        pass

    elif geno1.reference == geno2.coded and geno1.coded == geno2.reference:
        # Requires flipping.
        geno2 = flip_alleles(geno2)

    else:
        # Variants do not match.
        raise ValueError("Variants do not match ({} != {}).".format(
            geno1.variant, geno2.variant
        ))

    g1 = geno1.genotypes
    g2 = geno2.genotypes

    # Reorder samples and round values (if dosages)
    g1 = np.round(g1[idx1])
    g2 = np.round(g2[idx2])
    assert g1.shape == g2.shape

    missing_1 = np.isnan(g1)
    missing_2 = np.isnan(g2)
    missing = missing_1 | missing_2

    g1 = g1[~missing]
    g2 = g2[~missing]

    return {
        "n": g1.shape[0],
        "match": np.sum(g1 == g2),
        "mismatch": np.sum(g1 != g2),
        "missing_1": np.sum(missing_1),
        "missing_2": np.sum(missing_2),
    }
