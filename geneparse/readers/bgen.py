"""BGEN file reader based on PyBGEN."""

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


from pybgen import PyBGEN, ParallelPyBGEN

from .. import logging
from ..core import GenotypesReader, Genotypes, Variant


CHROM_STR_ENCODE = {"0{}".format(chrom): str(chrom) for chrom in range(1, 10)}
CHROM_STR_ENCODE["23"] = "X"
CHROM_STR_ENCODE["24"] = "Y"

CHROM_STR_DECODE = {v: k for k, v in CHROM_STR_ENCODE.items()}


class BGENReader(GenotypesReader):
    def __init__(self, filename, sample_filename=None, chromosome=None,
                 probability_threshold=0.9, cpus=1):
        """BGEN file reader.

        Args:
            filename (str): The name of the BGEN file.
            sample_filename (str): The name of the sample file (optional).
            probability_threshold (float): The probability threshold.

        """
        # The BGEN reader (parallel or no)
        if cpus == 1:
            self.is_parallel = False
            self._bgen = PyBGEN(filename, prob_t=probability_threshold)
        else:
            self.is_parallel = True
            self._bgen = ParallelPyBGEN(filename, prob_t=probability_threshold,
                                        cpus=cpus)

        # Getting the samples
        self.samples = self._bgen.samples
        if self.samples is None:
            # The BGEN file didn't contain samples, so we read from file
            if sample_filename is None:
                raise ValueError("No sample information in BGEN file, "
                                 "requires a 'sample_filename'")
            self._parse_sample_file(sample_filename)

        # Does the user ask for a chromosome?
        self.chrom = chromosome

    def close(self):
        self._bgen.close()

    def get_variant_genotypes(self, variant):
        """Get the genotypes from a well formed variant instance.

        Args:
            marker (Variant): A Variant instance.

        Returns:
            A list of Genotypes instance containing a pointer to the variant as
            well as a vector of encoded genotypes.

        """
        # The chromosome to search for (if a general one is set, that's the one
        # we need to search for)
        chrom = variant.chrom.name
        if self.chrom is not None and chrom == self.chrom:
            chrom = "NA"

        # Getting the results
        results = []
        iterator = self._bgen.iter_variants_in_region(
            CHROM_STR_DECODE.get(chrom, chrom), variant.pos, variant.pos,
        )
        for info, dosage in iterator:
            if (variant.alleles is None or
                    variant.iterable_alleles_eq([info.a1, info.a2])):
                results.append(Genotypes(
                    Variant(
                        info.name,
                        CHROM_STR_ENCODE.get(info.chrom, info.chrom),
                        info.pos, [info.a1, info.a2],
                    ),
                    dosage,
                    reference=info.a1,
                    coded=info.a2,
                    multiallelic=True,
                ))

        # If there are no results
        if not results:
            logging.variant_name_not_found(variant)

        return results

    def iter_genotypes(self):
        """Iterates on available markers.

        Returns:
            Genotypes instances.

        """
        for info, dosage in self._bgen.iter_variants():
            yield Genotypes(
                Variant(
                    info.name, CHROM_STR_ENCODE.get(info.chrom, info.chrom),
                    info.pos, [info.a1, info.a2],
                ),
                dosage,
                reference=info.a1,
                coded=info.a2,
                multiallelic=True,
            )

    def iter_variants(self):
        """Iterate over marker information."""
        for variant in self._bgen.iter_variant_info():
            yield Variant(
                variant.name,
                CHROM_STR_ENCODE.get(variant.chrom, variant.chrom),
                variant.pos, [variant.a1, variant.a2],
            )

    def get_variants_in_region(self, chrom, start, end):
        """Iterate over variants in a region."""
        if self.chrom is not None and chrom == self.chrom:
            # We are going to search for 'NA' since the chromosome was set
            chrom = "NA"

        iterator = self._bgen.iter_variants_in_region(
            CHROM_STR_DECODE.get(chrom, chrom), start, end,
        )
        for info, dosage in iterator:
            yield Genotypes(
                Variant(
                    info.name, CHROM_STR_ENCODE.get(info.chrom, info.chrom),
                    info.pos, [info.a1, info.a2],
                ),
                dosage,
                reference=info.a1,
                coded=info.a2,
                multiallelic=True,
            )

    def iter_variants_by_names(self, names):
        """Iterates over the genotypes for variants using a list of names.

        Args:
            names (list): The list of names for variant extraction.

        """
        if not self.is_parallel:
            yield from super().iter_variants_by_names(names)

        else:
            for info, dosage in self._bgen.iter_variants_by_names(names):
                yield Genotypes(
                    Variant(info.name,
                            CHROM_STR_ENCODE.get(info.chrom, info.chrom),
                            info.pos, [info.a1, info.a2]),
                    dosage,
                    reference=info.a1,
                    coded=info.a2,
                    multiallelic=True,
                )

    def get_variant_by_name(self, name):
        """Get the genotype of a marker using it's name.

        Args:
            name (str): The name of the marker.

        Returns:
            list: A list of Genotypes.

        """
        results = []

        try:
            for info, dosage in self._bgen.get_variant(name):
                results.append(Genotypes(
                    Variant(
                        info.name,
                        CHROM_STR_ENCODE.get(info.chrom, info.chrom),
                        info.pos,
                        [info.a1, info.a2],
                    ),
                    dosage,
                    reference=info.a1,
                    coded=info.a2,
                    multiallelic=False,
                ))

        except ValueError:
            logging.variant_name_not_found(name)

        return results

    def get_number_samples(self):
        """Returns the number of samples.

        Returns:
            int: The number of samples.

        """
        return self._bgen.nb_samples

    def get_number_variants(self):
        """Returns the number of markers.

        Returns:
            int: The number of markers.

        """
        return self._bgen.nb_variants

    def get_samples(self):
        """Returns the list of samples."""
        return list(self.samples)

    def _parse_sample_file(self, sample_filename):
        # Reading the sample file
        samples = []
        with open(sample_filename) as f:
            header = None
            for line in f:
                row = line.rstrip("\r\n").split(" ")

                if header is None:
                    header = {name: i for i, name in enumerate(row)}
                    for name in ("ID_1", "ID_2"):
                        if name not in header:
                            raise ValueError("{}: no column named {}".format(
                                name,
                            ))
                    continue

                # The first row is not a sample
                if row[:2] == ["0", "0"]:
                    continue

                samples.append((row[header["ID_1"]], row[header["ID_2"]]))

        # Checking ID_2 is unique
        id_2 = tuple(_[1] for _ in samples)
        if len(set(id_2)) == len(samples):
            self.samples = id_2

        else:
            logging.info(
                "Setting the index as 'fid_iid' because the individual IDs "
                "are not unique."
            )

            self.samples = tuple(
                "{}_{}".format(*_) for _ in samples
            )
