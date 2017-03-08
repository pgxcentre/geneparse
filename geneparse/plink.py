"""
Plink file reader based on PyPlink.
"""

import logging

import numpy as np
from pyplink import PyPlink

from .core import GenotypeReader, Variant, Genotypes


logger = logging.getLogger(__name__)


CHROM_MAP = {str(c): c for c in range(1, 23)}
CHROM_MAP["X"] = 23
CHROM_MAP["Y"] = 24


class PlinkReader(GenotypeReader):
    def __init__(self, prefix):
        """Binary plink file reader.
        Args:
            prefix (str): the prefix of the Plink binary files.

        """
        self.bed = PyPlink(prefix)
        self.bim = self.bed.get_bim()
        self.fam = self.bed.get_fam()

        # We want to set the index for the FAM file
        try:
            self.fam = self.fam.set_index("iid", verify_integrity=True)
        except ValueError:
            logger.info(
                "Setting the index as 'fid_iid' because the individual IDs "
                "are not unique."
            )

            self.fam["fid_iid"] = [
                "{fid}_{iid}".format(fid=fid, iid=iid)
                for fid, iid in zip(self.fam.fid, self.fam.iid)
            ]
            self.fam = self.fam.set_index("fid_iid", verify_integrity=True)

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
        try:
            idx = self.bim.loc[[variant.name], :]

            return [

            ]

        except KeyError:
            # Looking up variant in the bim.
            if variant.chrom == "X":
                idx = (
                    self.bim["chrom"].isin({23, 25}) &
                    (self.bim["pos"] == variant.pos)
                )

            else:
                idx = (
                    (self.bim["chrom"] == CHROM_MAP[variant.chrom]) &
                    (self.bim["pos"] == variant.pos)
                )

            if not idx.any():
                return []

            return [
                Genotypes(variant, self.bed.get_geno_marker(m), )
                for m in idx[idx].index
            ]

    def iter_genotypes(self):
        """Iterates on available markers.

        Returns:
            MarkerGenotypes: A named tuple containing the dataframe with the
            encoded genotypes for all samples (the index of the dataframe will
            be the sample IDs), the minor and major alleles.
        Note
        ====
            If the sample IDs are not unique, the index is changed to be the
            sample family ID and individual ID (i.e. fid_iid).
        """
        # Iterating over all markers
        for i, (_, genotypes) in enumerate(self.bed.iter_geno()):
            info = self.bim.iloc[i, :]

            yield Genotypes(
                Variant(info.name, info.chrom, info.pos, [info.a1, info.a2]),
                genotypes,
                reference=info.a2,
                coded=info.a1,
            )

    def iter_marker_info(self):
        """Iterate over marker information."""
        for idx, row in self.bim.iterrows():
            yield Variant(
                row.name, row.chrom, row.pos, [row.a1, row.a2]
            )

    def _create_genotypes(self, variant, genotypes):
        """Creates the genotype dataframe from an binary Plink file.

        Args:
            variant (str): The name of the marker.
            genotypes (numpy.ndarray): The genotypes.

        Returns:
            np.array: The genotypes.
        """
        # Getting and formatting the genotypes
        additive = self.create_geno_df(
            genotypes=genotypes,
            samples=self.fam.index,
        )

        # Checking the format is fine
        additive, minor, major = self.check_genotypes(
            genotypes=additive,
            minor=self.bim.loc[marker, "a1"],
            major=self.bim.loc[marker, "a2"],
        )

        chrom, pos = self.bim.loc[marker, ["chrom", "pos"]].values
        info = MarkerInfo(marker, chrom, pos, a1=minor, a2=major,
                          minor=MarkerInfo.A1)

        # Returning the value as ADDITIVE representation
        if self._representation == Representation.ADDITIVE:
            return MarkerGenotypes(info, additive)

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
