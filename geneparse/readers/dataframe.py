"""DataFrame file reader."""

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


from ..core import GenotypesReader, Genotypes, Variant
from .. import logging


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
            logging.variant_name_not_found(name)
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
