"""
Utilities
"""

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


import urllib
import json
import logging

import numpy as np

from .core import Variant


logger = logging.getLogger(__name__)


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


def rsids_to_variants(li):
    url = "http://grch37.rest.ensembl.org/variation/homo_sapiens"

    req = urllib.request.Request(
        url=url,
        data=json.dumps({"ids": li}).encode("utf-8"),
        headers={
            "Content-type": "application/json",
            "Accept": "application/json",
        },
        method="POST"
    )

    with urllib.request.urlopen(req) as f:
        data = json.loads(f.read().decode("utf-8"))

    out = {}
    for name, info in data.items():
        # Check the mappings.
        found = False
        for mapping in info["mappings"]:
            chrom = mapping.get("seq_region_name")
            pos = mapping.get("start")
            alleles = mapping.get("allele_string").split("/")

            assembly = mapping.get("assembly_name")

            valid = (assembly == "GRCh37" and
                     chrom is not None and
                     pos is not None and
                     len(alleles) >= 2)

            if found and valid:
                logger.warning("Multiple mappings for '{}'.".format(name))
            elif valid:
                found = True
                out[name] = Variant(name, chrom, pos, alleles)

        if not found:
            logger.warning(
                "Could not find mappings for '{}'.".format(name)
            )

    return out
