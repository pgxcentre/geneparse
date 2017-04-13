"""
DataFrame file reader.
"""


import logging

from .core import GenotypesReader, Genotypes, Variant


logger = logging.getLogger(__name__)


class DataFrameReader(GenotypesReader):
    def __init__(self, dataframe, map_info):
        """Reads genotypes from a pandas DataFrame.

        Args:
            dataframe (pandas.DataFrame): The data.
            map_info (pandas.DataFrame): The mapping information.

        Note
        ====
            The index of the dataframe should be the sample IDs. The index of
            the map_info should be the variant name, and there should be
            columns named chrom and pos.

        """
        self.df = dataframe
        self.map_info = map_info

    def iter_genotypes(self):
        """Iterates on available markers.

        Returns:
            Genotypes instances.

        """
        # Parsing each column of the dataframe
        for variant in self.df.columns:
            genotypes = self.df.loc[:, variant].values
            info = self.map_info.loc[variant, :]
            print("\n%%%%%%%%%%%%%%%%")
            print(info)

            yield Genotypes(
                Variant(info.name, info.chrom, info.pos, [info.a1, info.a2]),
                genotypes,
                reference=info.a2,
                coded=info.a1,
                multiallelic=False,
            )

    def get_variant_by_name(self, name):
        """Get the genotypes for a given variant (by name).

        Args:
            name (str): The name of the variant to retrieve the genotypes.

        Returns:
            list: A list of Genotypes. This is a list in order to keep the same
            behaviour as the other functions.

        """
        try:
            geno = self.df.loc[:, name].values
            info = self.map_info.loc[name, :]

        except KeyError:
            # The variant is not in the data, so we return an empty
            # list
            logger.warning("Variant {} was not found".format(name))
            return []

        else:
            return [Genotypes(
                Variant(info.name, info.chrom, info.pos, [info.a1, info.a2]),
                geno,
                reference=info.a2,
                coded=info.a1,
                multiallelic=False,
            )]

    def get_samples(self):
        """Get an ordered collection of the samples in the genotype container.
        """
        return self.df.index.tolist()

    def get_number_samples(self):
        """Return the number of samples."""
        return self.df.shape[0]

    def get_number_variants(self):
        """Return the number of variants in the file."""
        return self.df.shape[1]
