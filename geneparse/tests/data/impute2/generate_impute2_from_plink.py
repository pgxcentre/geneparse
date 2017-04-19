#!/usr/bin/env python


import re
import itertools

from pyplink import PyPlink


def create_probs_from_genotypes(genotype):
    """Creates probabilities from an additive genotype."""
    if genotype == 0:
        return (1, 0, 0)

    if genotype == 1:
        return (0, 1, 0)

    if genotype == 2:
        return (0, 0, 1)

    if genotype == -1:
        # Lower than normal probabilities
        return (0.8, 0.1, 0.1)


with PyPlink("../plink/btest") as bed, \
        open("impute2_test.impute2", "w") as impute2_f, \
        open("impute2_test.sample", "w") as impute2_s:
    # Getting the FAM and the BIM
    fam = bed.get_fam()
    bim = bed.get_bim()

    # Generating the IMPUTE2 file
    for v, genotypes in bed.iter_geno():
        info = bim.loc[v, :]
        assert v == info.name

        r = re.search(r"(:dup[0-9]+)$", v)
        if r:
            v = v.replace(r.group(1), "")

        probs = [create_probs_from_genotypes(g) for g in genotypes]
        print(info.chrom, v, info.pos, info.a2, info.a1,
              *itertools.chain(*probs), sep=" ", file=impute2_f)

    # Generating the SAMPLE file
    print("ID_1", "ID_2", "missing", "father", "mother", "sex", "plink_pheno",
          file=impute2_s)
    print(0, 0, 0, "D", "D", "D", "B", file=impute2_s)
    for _, sample_info in fam.iterrows():
        print(sample_info.fid, sample_info.iid, 0, sample_info.father,
              sample_info.mother, sample_info.gender, sample_info.status,
              file=impute2_s)
