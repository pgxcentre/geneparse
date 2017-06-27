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
import sqlite3
from os import path
from struct import unpack

import numpy as np

from ..core import GenotypesReader, Genotypes, Variant

try:
    import zstd
    HAS_ZSTD = True
except ImportError:
    HAS_ZSTD = False


CHROM_STR_ENCODE = {"0{}".format(chrom): str(chrom) for chrom in range(1, 10)}
CHROM_STR_ENCODE["23"] = "X"
CHROM_STR_ENCODE["24"] = "Y"


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
            # TODO: Parse sample file here
            if sample_filename is None:
                raise ValueError("No sample information in BGEN file, "
                                 "requires a 'sample_filename'")

        else:
            self._parse_sample_identifier_block()

        # If we have an index, we open it (sqlite3)
        if not path.isfile(filename + ".bgi"):
            raise ValueError("{}: no index file".format(filename))
        self._bgen_index = sqlite3.connect(filename + ".bgi")

        # The probability threshold
        self.prob_t = probability_threshold

    def close(self):
        if self._bgen_file:
            self._bgen_file.close()
        if self._bgen_index:
            self._bgen_index.close()

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
        self._is_compressed = False
        if compression == 0:
            # No decompression required
            self._decompress = self._no_decompress

        elif compression == 1:
            # ZLIB decompression
            self._decompress = zlib.decompress
            self._is_compressed = True

        elif compression == 2:
            if not HAS_ZSTD:
                raise ValueError("zstandard module is not installed")

            # ZSTANDARD decompression (needs to be check)
            self._decompress = zstd.ZstdDecompressor().decompress
            self._is_compressed = True

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

    def _seek_to_first_variant(self):
        """Seeks to the first variant of the file."""
        c = self._bgen_index.cursor()
        c.execute(
            "SELECT file_start_position "
            "FROM Variant "
            "ORDER BY file_start_position "
            "LIMIT 1",
        )
        self._bgen_file.seek(c.fetchone()[0])

    def get_variant_genotypes(self, variant):
        """Get the genotypes from a well formed variant instance.

        Args:
            marker (Variant): A Variant instance.

        Returns:
            A list of Genotypes instance containing a pointer to the variant as
            well as a vector of encoded genotypes.

        """
        pass

    def _get_curr_variant_info(self):
        """Gets the current variant's information."""
        if self._layout == 1:
            n = unpack("I", self._bgen_file.read(4))[0]
            if n != self.nb_samples:
                raise ValueError(
                    "{}: invalid BGEN file".format(self._bgen_file.name),
                )

        # Reading the variant id
        var_id = self._bgen_file.read(
            unpack("H", self._bgen_file.read(2))[0]
        ).decode()

        # Reading the variant rsid
        rs_id = self._bgen_file.read(
            unpack("H", self._bgen_file.read(2))[0]
        ).decode()

        # Reading the chromosome
        chrom = self._bgen_file.read(
            unpack("H", self._bgen_file.read(2))[0]
        ).decode()

        # Reading the position
        pos = unpack("I", self._bgen_file.read(4))[0]

        # Getting the number of alleles
        nb_alleles = 2
        if self._layout == 2:
            nb_alleles = unpack("H", self._bgen_file.read(2))[0]

        # Getting the alleles
        alleles = []
        for _ in range(nb_alleles):
            alleles.append(self._bgen_file.read(
                unpack("I", self._bgen_file.read(4))[0]
            ).decode())

        return var_id, rs_id, chrom, pos, tuple(alleles)

    def _get_curr_variant_dosage(self):
        """Gets the current variant's dosage."""
        dosage = None
        if self._layout == 1:
            c = self.nb_samples
            if self._is_compressed:
                c = unpack("I", self._bgen_file.read(4))[0]

            # Getting the probabilities
            probs = np.fromstring(
                self._decompress(self._bgen_file.read(c)),
                dtype="u2",
            ) / 32768
            probs.shape = (self.nb_samples, 3)

            # Computing the dosage
            dosage = self._layout_1_probs_to_dosage(probs)

        else:
            # The total length C of the rest of the data for this variant
            c = unpack("I", self._bgen_file.read(4))[0]

            # The number of bytes to read
            to_read = c

            # D = C if no compression
            d = c
            if self._is_compressed:
                # The total length D of the probability data after
                # decompression
                d = unpack("I", self._bgen_file.read(4))[0]
                to_read = c - 4

            # Reading the data and checking
            data = self._decompress(self._bgen_file.read(to_read))
            if len(data) != d:
                raise ValueError(
                    "{}: invalid BGEN file".format(self._bgen_file.name)
                )

            # Checking the number of samples
            n = unpack("I", data[:4])[0]
            if n != self.nb_samples:
                raise ValueError(
                    "{}: invalid BGEN file".format(self._bgen_file.name)
                )
            data = data[4:]

            # Checking the number of alleles (we only accept 2 alleles)
            nb_alleles = unpack("H", data[:2])[0]
            if nb_alleles != 2:
                raise ValueError(
                    "{}: only two alleles are "
                    "supported".format(self._bgen_file.name)
                )
            data = data[2:]

            # TODO: Check ploidy for sexual chromosomes
            # The minimum and maximum for ploidy (we only accept ploidy of 2)
            min_ploidy, max_ploidy = data[:2]
            if min_ploidy != 2 and max_ploidy != 2:
                raise ValueError(
                    "{}: only accepting ploidy of "
                    "2".format(self._bgen_file.name)
                )
            data = data[2:]

            # Check the list of N bytes for missingness (since we assume only
            # diploid values for each sample)
            ploidy_info = data[:n]
            ploidy_info = np.unpackbits(
                np.array([[_] for _ in ploidy_info], dtype=np.uint8),
                axis=1,
            )
            missing_data = ploidy_info[:, 0] == 1
            data = data[n:]

            # TODO: Permit phased data
            # Is the data phased?
            is_phased = data[0] == 1
            if is_phased:
                raise ValueError(
                    "{}: only accepting unphased "
                    "data".format(self._bgen_file.name)
                )
            data = data[1:]

            # The number of bytes used to encode each probabilities
            b = data[0]
            if b % 8 != 0:
                raise ValueError(
                    "{}: only multuple of 8 encoding is "
                    "accepted".format(self._bgen_file.name)
                )
            data = data[1:]

            # Reading the probabilities (don't forget we allow only for diploid
            # values)
            probs = np.fromstring(
                data, dtype="u{}".format(b // 8)
            ) / (2**b - 1)
            probs.shape = (self.nb_samples, 2)

            # Computing the dosage
            dosage = self._layout_2_probs_to_dosage(probs)

            # Setting the missing to NaN
            dosage[missing_data] = np.nan

        return dosage

    def _layout_1_probs_to_dosage(self, probs):
        """Transforms probability values to dosage (from layout 1)"""
        # Constructing the dosage
        dosage = 2 * probs[:, 2] + probs[:, 1]
        if self.prob_t > 0:
            dosage[~np.any(probs >= self.prob_t, axis=1)] = np.nan

        return dosage

    def _layout_2_probs_to_dosage(self, probs):
        """Transforms probability values to dosage (from layout 2)."""
        # Computing the last genotype's probabilities
        last_probs = 1 - np.sum(probs, axis=1)

        # Constructing the dosage
        dosage = 2 * last_probs + probs[:, 1]

        # Setting low quality to NaN
        if self.prob_t > 0:
            good_probs = (
                np.any(probs >= self.prob_t, axis=1) |
                (last_probs >= self.prob_t)
            )
            dosage[~good_probs] = np.nan

        return dosage

    def iter_genotypes(self):
        """Iterates on available markers.

        Returns:
            Genotypes instances.

        """
        # Seeking at the beginning of the variants
        self._seek_to_first_variant()

        # Checking the number of samples
        for i in range(self.nb_variants):
            # Getting the variant's information
            var_id, rs_id, chrom, pos, alleles = self._get_curr_variant_info()

            # Getting the variant's dosage
            dosage = self._get_curr_variant_dosage()

            # TODO: Check multiallelic status below
            yield Genotypes(
                Variant(rs_id, CHROM_STR_ENCODE.get(chrom, chrom), pos,
                        alleles),
                dosage,
                reference=alleles[0],
                coded=alleles[1],
                multiallelic=False,
            )

    def iter_variants(self):
        """Iterate over marker information."""
        c = self._bgen_index.cursor()
        c.execute(
            "SELECT chromosome, position, rsid, allele1, allele2 "
            "FROM Variant "
        )

        # The array size
        array_size = 10000

        # Fetching the results
        results = c.fetchmany(array_size)
        while results:
            for chrom, pos, rsid, a1, a2 in results:
                yield Variant(rsid, CHROM_STR_ENCODE.get(chrom, chrom),
                              pos, [a1, a2])
            results = c.fetchmany(array_size)

    def get_variants_in_region(self, chrom, start, end):
        """Iterate over variants in a region."""
        pass

    def get_variant_by_name(self, name):
        """Get the genotype of a marker using it's name.

        Args:
            name (str): The name of the marker.

        Returns:
            list: A list of Genotypes.

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
        """Returns the list of samples."""
        return list(self.samples)
