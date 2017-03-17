"""
Plink file reader based on PyPlink.
"""

import logging

from pyplink import PyPlink
import numpy as np

from .core import GenotypesReader, Variant, Genotypes


logger = logging.getLogger(__name__)


CHROM_STR_TO_INT = {str(c): c for c in range(1, 23)}
CHROM_STR_TO_INT["X"] = 23
CHROM_STR_TO_INT["Y"] = 24


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
            logger.info(
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
        plink_chrom = CHROM_STR_TO_INT[variant.chrom]
        info = self.bim.loc[
            (self.bim.chrom == plink_chrom) &
            (self.bim.pos == variant.pos), :
        ]

        if info.shape[0] == 0:
            raise KeyError(variant)

        elif info.shape[0] == 1:
            return self._get_biallelic_variant(variant, info)

        else:
            return self._get_multialleic_variant(variant, info)

    def _get_biallelic_variant(self, variant, info):
        # From 1.3.2 onwards, PyPlink sets unique names.
        info = info.iloc[0, :]
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
                        info.loc[[name], :]
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

    @staticmethod
    def _normalize_missing(g):
        """Normalize a plink genotype vector."""
        g = g.astype(float)
        g[g == -1.0] = np.nan
        return g
