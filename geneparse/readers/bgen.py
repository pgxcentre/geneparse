"""
BGEN file reader.
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


from ..core import GenotypesReader


class BGENReader(GenotypesReader):
    def __init__(self, filename, sample_filename=None,
                 probability_threshold=0.9):
        """BGEN file reader.

        Args:
            filename (str): The name of the BGEN file.
            sample_filename (str): The name of the sample file (optional).
            probability_threshold (float): The probability threshold.

        """
        # The BGEN file
        self._bgen_file = None

        # If we have an index, we open it (sqlite3)
        self._bgen_index = None

        # The probability threshold
        self.prob_t = probability_threshold

    def close(self):
        if self._bgen_file:
            self._bgen_file.close()

    def get_variant_genotypes(self, variant):
        """Get the genotypes from a well formed variant instance.

        Args:
            marker (Variant): A Variant instance.

        Returns:
            A list of Genotypes instance containing a pointer to the variant as
            well as a vector of encoded genotypes.

        """
        pass

    def iter_genotypes(self):
        """Iterates on available markers.

        Returns:
            Genotypes instances.

        """
        pass

    def iter_variants(self):
        """Iterate over marker information."""
        pass

    def get_variants_in_region(self, chrom, start, end):
        """Iterate over variants in a region."""
        pass

    def get_variant_by_name(self, name, variant_info=None):
        """Get the genotype of a marker using it's name.

        Args:
            name (str): The name of the marker.
            variant_info (pandas.Series): The marker information (e.g. seek).

        Returns:
            list: A list of Genotypes (only one for PyPlink, see note below).

        Note
        ====
            From PyPlink version 1.3.2 and onwards, each name is unique in the
            dataset. Hence, we can use the 'get_geno_marker' function and be
            sure only one variant is returned.

        """
        pass

    def get_number_samples(self):
        """Returns the number of samples.

        Returns:
            int: The number of samples.

        """
        pass

    def get_number_variants(self):
        """Returns the number of markers.

        Returns:
            int: The number of markers.

        """
        pass

    def get_samples(self):
        pass
