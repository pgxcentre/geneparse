"""Plink file reader based on PyPlink."""

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


from pyplink import PyPlink
import numpy as np

from ..core import GenotypesReader, Variant, Genotypes
from .. import logging


CHROM_STR_TO_INT = {str(c): c for c in range(1, 23)}
CHROM_STR_TO_INT["X"] = 23
CHROM_STR_TO_INT["Y"] = 24
CHROM_STR_TO_INT["XY"] = 25
CHROM_STR_TO_INT["MT"] = 26
CHROM_STR_TO_INT["Unknown"] = 0  # TODO What is plink chromosome 0?


CHROM_INT_TO_STR = {v: k for k, v in CHROM_STR_TO_INT.items()}


class PlinkReader(GenotypesReader):
    def __init__(self, prefix):
        """Binary plink file reader.
        Args:
            prefix (str): the prefix of the Plink binary files.

        """
        self.bed = PyPlink(prefix)
        self.bim = self.bed.get_bim()
        self.fam = self.bed.get_fam()

        # Identify all multi-allelics.
        self.bim["multiallelic"] = False
        self.bim.loc[
            self.bim.duplicated(["chrom", "pos"], keep=False),
            "multiallelic"
        ] = True

        # We want to set the index for the FAM file
        try:
            self.fam = self.fam.set_index("iid", verify_integrity=True)
        except ValueError:
            logging.info(
                "Setting the index as 'fid_iid' because the individual IDs "
                "are not unique."
            )

            self.fam["fid_iid"] = [
                "{fid}_{iid}".format(fid=fid, iid=iid)
                for fid, iid in zip(self.fam.fid, self.fam.iid)
            ]
            self.fam = self.fam.set_index("fid_iid", verify_integrity=True)

    def close(self):
        self.bed.close()

    def get_variant_genotypes(self, variant):
        """Get the genotypes from a well formed variant instance.

        Args:
            marker (Variant): A Variant instance.

        Returns:
            A list of Genotypes instance containing a pointer to the variant as
            well as a vector of encoded genotypes.

        Note
        ====
            If the sample IDs are not unique, the index is changed to be the
            sample family ID and individual ID (i.e. fid_iid).

        """
        # Find the variant in the bim.
        try:
            plink_chrom = CHROM_STR_TO_INT[variant.chrom.name]
        except KeyError:
            raise ValueError(
                "Invalid chromosome ('{}') for Plink.".format(variant.chrom)
            )

        info = self.bim.loc[
            (self.bim.chrom == plink_chrom) &
            (self.bim.pos == variant.pos), :
        ]

        if info.shape[0] == 0:
            logging.variant_not_found(variant)
            return []

        elif info.shape[0] == 1:
            return self._get_biallelic_variant(variant, info)

        else:
            return self._get_multialleic_variant(variant, info)

    def _get_biallelic_variant(self, variant, info, _check_alleles=True):
        # From 1.3.2 onwards, PyPlink sets unique names.
        info = info.iloc[0, :]
        variant_alleles = variant._encode_alleles([info.a2, info.a1])
        if (_check_alleles and variant_alleles != variant.alleles):
            # Variant with requested alleles is unavailable.
            logging.variant_not_found(variant)
            return []

        geno = self._normalize_missing(self.bed.get_geno_marker(info.name))
        return [Genotypes(variant, geno, info.a2, info.a1, False)]

    def _get_multialleic_variant(self, variant, info):
        # Check if alleles are specified.
        out = []
        if variant.alleles is None:
            # If no alleles are specified, we return all the possible
            # bi-allelic variats.
            for name, row in info.iterrows():
                geno = self.bed.get_geno_marker(name)
                geno = self._normalize_missing(geno)
                out.append(Genotypes(
                    variant, geno, row.a2, row.a1, True
                ))

        else:
            # Find the requested alleles.
            for name, row in info.iterrows():
                row_alleles = set(Variant._encode_alleles((row.a1, row.a2)))
                if row_alleles.issubset(variant.alleles_set):
                    out.extend(self._get_biallelic_variant(
                        variant,
                        info.loc[[name], :],
                        _check_alleles=False
                    ))

        return out

    def iter_genotypes(self):
        """Iterates on available markers.

        Returns:
            Genotypes instances.

        Note
        ====
            If the sample IDs are not unique, the index is changed to be the
            sample family ID and individual ID (i.e. fid_iid).

        """
        # Iterating over all markers
        for i, (_, genotypes) in enumerate(self.bed.iter_geno()):
            info = self.bim.iloc[i, :]

            yield Genotypes(
                Variant(info.name, CHROM_INT_TO_STR[info.chrom],
                        info.pos, [info.a1, info.a2]),
                self._normalize_missing(genotypes),
                reference=info.a2,
                coded=info.a1,
                multiallelic=info.multiallelic
            )

    def iter_variants(self):
        """Iterate over marker information."""
        for idx, row in self.bim.iterrows():
            yield Variant(
                row.name, CHROM_INT_TO_STR[row.chrom], row.pos,
                [row.a1, row.a2]
            )

    def get_variants_in_region(self, chrom, start, end):
        """Iterate over variants in a region."""
        bim = self.bim.loc[
            (self.bim["chrom"] == CHROM_STR_TO_INT[chrom]) &
            (start <= self.bim["pos"]) &
            (self.bim["pos"] <= end)
        ]
        for i, g in enumerate(self.bed.iter_geno_marker(bim.index)):
            info = bim.iloc[i, :]
            name, geno = g
            yield Genotypes(
                Variant(info.name, CHROM_INT_TO_STR[info.chrom],
                        info.pos, [info.a1, info.a2]),
                self._normalize_missing(geno),
                reference=info.a2,
                coded=info.a1,
                multiallelic=info.multiallelic
            )

    def get_variant_by_name(self, name):
        """Get the genotype of a marker using it's name.

        Args:
            name (str): The name of the marker.

        Returns:
            list: A list of Genotypes (only one for PyPlink, see note below).

        Note
        ====
            From PyPlink version 1.3.2 and onwards, each name is unique in the
            dataset. Hence, we can use the 'get_geno_marker' function and be
            sure only one variant is returned.

        """
        # From 1.3.2 onwards, PyPlink sets unique names.
        # Getting the genotypes
        try:
            geno, i = self.bed.get_geno_marker(name, return_index=True)

        except ValueError:
            if name in self.bed.get_duplicated_markers():
                # The variant is a duplicated one, so we go through all the
                # variants with the same name and the :dupx suffix
                return [
                    self.get_variant_by_name(dup_name).pop()
                    for dup_name in self.bed.get_duplicated_markers()[name]
                ]

            else:
                # The variant is not in the BIM file, so we return an empty
                # list
                logging.variant_name_not_found(name)
                return []

        else:
            info = self.bim.iloc[i, :]
            return [Genotypes(
                Variant(info.name, CHROM_INT_TO_STR[info.chrom], info.pos,
                        [info.a1, info.a2]),
                self._normalize_missing(geno),
                reference=info.a2,
                coded=info.a1,
                multiallelic=info.multiallelic,
            )]

    def get_number_samples(self):
        """Returns the number of samples.
        Returns:
            int: The number of samples.
        """
        return self.bed.get_nb_samples()

    def get_number_variants(self):
        """Returns the number of markers.
        Returns:
            int: The number of markers.
        """
        return self.bed.get_nb_markers()

    def get_samples(self):
        return list(self.fam.index)

    @staticmethod
    def _normalize_missing(g):
        """Normalize a plink genotype vector."""
        g = g.astype(float)
        g[g == -1.0] = np.nan
        return g
