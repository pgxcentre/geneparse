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


import zlib
from struct import unpack

import numpy as np

import zstd

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
        self._bgen_file = open(filename, "rb")
        self._parse_header_block()

        # Getting the sample
        if not self._has_sample:
            if sample_filename is None:
                raise ValueError("No sample information in BGEN file, "
                                 "requires a 'sample_filename'")
            # TODO: Parse sample file here

        else:
            self._parse_sample_identifier_block()

        # If we have an index, we open it (sqlite3)
        self._bgen_index = None

        # The probability threshold
        self.prob_t = probability_threshold

    def close(self):
        if self._bgen_file:
            self._bgen_file.close()

    def _parse_header_block(self):
        """Parses the BGEN file header."""
        # Getting the data offset (the start point of the data
        self._offset = unpack("I", self._bgen_file.read(4))[0]

        # Getting the header size
        self._header_size = unpack("I", self._bgen_file.read(4))[0]

        # Getting the number of samples and variants
        self.nb_variants = unpack("I", self._bgen_file.read(4))[0]
        self.nb_samples = unpack("I", self._bgen_file.read(4))[0]

        # Checking the magic number
        magic = self._bgen_file.read(4)
        if magic != b"bgen":
            # The magic number might be 0, then
            if unpack("I", magic)[0] != 0:
                raise ValueError(
                    "{}: invalid BGEN file.".format(self._bgen_file.name)
                )

        # Passing through the "free data area"
        self._bgen_file.read(self._header_size - 20)

        # Reading the flag
        flag = np.unpackbits(
            np.array([[_] for _ in self._bgen_file.read(4)], dtype=np.uint8),
            axis=1,
        )

        # Getting the compression type from the layout
        compression = self._bits_to_int(flag[0, -2:])
        if compression == 0:
            # No decompression required
            self._decompress = self._no_decompress
        elif compression == 1:
            # ZLIB decompression
            self._decompress = zlib.decompress
        elif compression == 2:
            # ZSTANDARD decompression (needs to be check)
            self._decompress = zstd.ZstdDecompressor().decompress

        # Getting the layout
        layout = self._bits_to_int(flag[0, -6:-2])
        if layout == 0:
            raise ValueError(
                "{}: invalid BGEN file".format(self._bgen_file.name)
            )
        elif layout == 1:
            self._layout = 1
        elif layout == 2:
            self._layout = 2
        else:
            raise ValueError(
                "{}: {} invalid layout type".format(self._bgen_file.name,
                                                    layout)
            )

        # Checking if samples are in the file
        self._has_sample = flag[-1, 0] == 1

    def _parse_sample_identifier_block(self):
        """Parses the sample identifier block."""
        # Getting the block size
        block_size = unpack("I", self._bgen_file.read(4))[0]
        if block_size + self._header_size > self._offset:
            raise ValueError(
                "{}: invalid BGEN file".format(self._bgen_file.name)
            )

        # Checking the number of samples
        n = unpack("I", self._bgen_file.read(4))[0]
        if n != self.nb_samples:
            raise ValueError(
                "{}: invalid BGEN file".format(self._bgen_file.name)
            )

        # Getting the sample information
        samples = []
        for i in range(self.nb_samples):
            size = unpack("H", self._bgen_file.read(2))[0]
            samples.append(self._bgen_file.read(size).decode())
        self.samples = tuple(samples)

    @staticmethod
    def _bits_to_int(bits):
        """Converts bits to int."""
        result = 0
        for bit in bits:
            result = (result << 1) | bit
        return result

    def _no_decompress(data):
        """No compression, so we return the data as is."""
        return data

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
        return self.nb_samples

    def get_number_variants(self):
        """Returns the number of markers.

        Returns:
            int: The number of markers.

        """
        return self.nb_variants

    def get_samples(self):
        return self.samples
