"""IMPUTE2 file reader."""

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


from os import path
from collections import Counter

import numpy as np
import pandas as pd

from ..index.impute2 import get_open_func, get_index

from ..core import GenotypesReader, Variant, Genotypes, VALID_CHROMOSOMES
from .. import logging


CHROM_STR_TO_INT = {str(c): c for c in range(1, 23)}
CHROM_STR_TO_INT["X"] = 23
CHROM_STR_TO_INT["Y"] = 24
CHROM_STR_TO_INT["XY"] = 25
CHROM_STR_TO_INT["MT"] = 26
CHROM_STR_TO_INT["Unknown"] = 0


CHROM_STR_ENCODE = {
    "23": VALID_CHROMOSOMES["X"],
    "24": VALID_CHROMOSOMES["Y"],
    "25": VALID_CHROMOSOMES["XY"],
    "26": VALID_CHROMOSOMES["MT"]
}


class Impute2Reader(GenotypesReader):
    def __init__(self, filename, sample_filename, probability_threshold=0.9):
        """IMPUTE2 file reader.

        Args:
            filename (str): The name of the IMPUTE2 file.
            sample_filename (str): The name of the SAMPLE file.
            probability_threshold (float): The probability threshold.

        Note
        ====
            If the sample IDs are not unique, the index is changed to be the
            sample family ID and individual ID (i.e. fid_iid).

        """
        # Reading the samples
        self.samples = pd.read_csv(sample_filename, sep=" ", skiprows=2,
                                   names=["fid", "iid", "missing", "father",
                                          "mother", "sex", "plink_geno"],
                                   dtype=dict(fid=str, iid=str))

        # We want to set the index for the samples
        try:
            self.samples = self.samples.set_index("iid", verify_integrity=True)

        except ValueError:
            logging.info(
                "Setting the index as 'fid_iid' because the individual IDs "
                "are not unique."
            )

            self.samples["fid_iid"] = [
                "{fid}_{iid}".format(fid=fid, iid=iid)
                for fid, iid in zip(self.samples.fid, self.samples.iid)
            ]
            self.samples = self.samples.set_index(
                "fid_iid", verify_integrity=True,
            )

        # The IMPUTE2 file
        self._impute2_file = get_open_func(filename)(filename, "r")

        # If we have an index, we read it
        self.has_index = path.isfile(filename + ".idx")
        self._impute2_index = None
        self._index_has_location = False
        if self.has_index:
            self._impute2_index = get_index(
                filename,
                cols=[0, 1, 2],
                names=["chrom", "name", "pos"],
                sep=" ",
            )

            # Checking for duplicated marker iD
            try:
                self._impute2_index = self._impute2_index.set_index(
                    "name", verify_integrity=True,
                )
                self._has_duplicated = False

            except ValueError as e:
                self._has_duplicated = True

                # Finding the duplicated markers
                duplicated = self._impute2_index.name.duplicated(keep=False)
                duplicated_markers = self._impute2_index.loc[
                    duplicated, "name"
                ]
                duplicated_marker_counts = duplicated_markers.value_counts()

                # The dictionary that will contain information about the
                # duplicated markers
                self._dup_markers = {
                    m: [] for m in duplicated_marker_counts.index
                }

                # Logging a warning
                logging.found_duplicates(duplicated_marker_counts.iteritems())

                # Renaming the markers
                counter = Counter()
                for i, marker in duplicated_markers.iteritems():
                    counter[marker] += 1
                    new_name = "{}:dup{}".format(marker, counter[marker])
                    self._impute2_index.loc[i, "name"] = new_name

                    # Updating the dictionary containing the duplicated markers
                    self._dup_markers[marker].append(new_name)

                # Resetting the index
                self._impute2_index = self._impute2_index.set_index(
                    "name", verify_integrity=True,
                )

            # Checking if we have chrom/pos in the index
            self._index_has_location = (
                "chrom" in self._impute2_index.columns and
                "pos" in self._impute2_index.columns
            )
            if self._index_has_location:
                # Setting the multiallelic values
                self._impute2_index["multiallelic"] = False
                self._impute2_index.loc[
                    self._impute2_index.duplicated(["chrom", "pos"],
                                                   keep=False),
                    "multiallelic"
                ] = True

        # Saving the probability threshold
        self.prob_t = probability_threshold

    def get_duplicated_markers(self):
        """Returns the duplicated markers, if any.

        Args:
            dict: The map for duplicated marker (might be empty).

        """
        if self._has_duplicated:
            return self._dup_markers
        else:
            return {}

    def close(self):
        if self._impute2_file:
            self._impute2_file.close()

    def get_variant_genotypes(self, variant):
        """Get the genotypes from a well formed variant instance.

        Args:
            marker (Variant): A Variant instance.

        Returns:
            A list of Genotypes instance containing a pointer to the variant as
            well as a vector of encoded genotypes.

        """
        if not self.has_index:
            raise NotImplementedError("Not implemented when IMPUTE2 file is "
                                      "not indexed (see genipe)")

        # Find the variant in the index
        try:
            impute2_chrom = CHROM_STR_TO_INT[variant.chrom.name]
        except KeyError:
            raise ValueError(
                "Invalid chromosome ('{}') for IMPUTE2.".format(variant.chrom)
            )

        variant_info = self._impute2_index[
            (self._impute2_index.chrom == impute2_chrom) &
            (self._impute2_index.pos == variant.pos)
        ]

        if variant_info.shape[0] == 0:
            logging.variant_not_found(variant)
            return []

        elif variant_info.shape[0] == 1:
            return self._get_biallelic_variant(variant, variant_info)

        else:
            return self._get_multialleic_variant(variant, variant_info)

    def _get_biallelic_variant(self, variant, info, _check_alleles=True):
        """Creates a bi-allelic variant."""
        info = info.iloc[0, :]
        assert not info.multiallelic

        # Seeking and parsing the file
        self._impute2_file.seek(info.seek)
        genotypes = self._parse_impute2_line(self._impute2_file.readline())

        variant_alleles = variant._encode_alleles([
            genotypes.reference, genotypes.coded,
        ])
        if (_check_alleles and variant_alleles != variant.alleles):
            # Variant with requested alleles is unavailable.
            logging.variant_not_found(variant)
            return []

        return [genotypes]

    def _get_multialleic_variant(self, variant, info):
        # Check if alleles are specified.
        out = []
        if variant.alleles is None:
            # If no alleles are specified, we return all the possible
            # bi-allelic variants.
            for name, row in info.iterrows():
                assert row.multiallelic

                # Seeking and parsing the file
                self._impute2_file.seek(row.seek)
                genotypes = self._parse_impute2_line(
                    self._impute2_file.readline(),
                )

                # fixing
                self._fix_genotypes_object(genotypes, row)

                out.append(genotypes)

        else:
            # Find the requested alleles.
            for name, row in info.iterrows():
                assert row.multiallelic

                # Seeking and parsing the file
                self._impute2_file.seek(row.seek)
                genotypes = self._parse_impute2_line(
                    self._impute2_file.readline(),
                )

                # Checking the alleles
                row_alleles = set(Variant._encode_alleles(
                    (genotypes.reference, genotypes.coded),
                ))
                if row_alleles.issubset(variant.alleles_set):
                    # Fixing
                    self._fix_genotypes_object(genotypes, row)
                    out.append(genotypes)

        return out

    def iter_genotypes(self):
        """Iterates on available markers.

        Returns:
            Genotypes instances.

        """
        # Seeking at the beginning of the file
        self._impute2_file.seek(0)

        # Parsing each lines of the IMPUTE2 file
        for i, line in enumerate(self._impute2_file):
            genotypes = self._parse_impute2_line(line)

            variant_info = None
            if self.has_index:
                variant_info = self._impute2_index.iloc[i, :]
            self._fix_genotypes_object(genotypes, variant_info)

            yield genotypes

    def iter_variants(self):
        """Iterate over marker information."""
        if not self.has_index:
            raise NotImplementedError("Not implemented when IMPUTE2 file is "
                                      "not indexed (see genipe)")

        for name, row in self._impute2_index.iterrows():
            # Seeking to the right place in the file
            f = self._impute2_file
            f.seek(int(row.seek))
            chrom, name, pos, a1, a2 = f.read(1024).split(" ")[:5]
            pos = int(pos)

            yield Variant(name, CHROM_STR_ENCODE.get(chrom, chrom), pos,
                          [a1, a2])

    def get_variants_in_region(self, chrom, start, end):
        """Iterate over variants in a region."""
        if not self.has_index:
            raise NotImplementedError("Not implemented when IMPUTE2 file is "
                                      "not indexed (see genipe)")

        if not self._index_has_location:
            raise NotImplementedError("Not implemented when index doesn't "
                                      "have location information.")

        # Getting the required variants
        required = self._impute2_index.loc[
            (self._impute2_index.chrom == CHROM_STR_TO_INT[chrom]) &
            (start <= self._impute2_index.pos) &
            (self._impute2_index.pos <= end)
        ]

        for name, variant_info in required.iterrows():
            for genotypes in self.get_variant_by_name(name, variant_info):
                self._fix_genotypes_object(genotypes, variant_info)
                yield genotypes

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
        # From 1.3.2 onwards, PyPlink sets unique names.
        if not self.has_index:
            raise NotImplementedError("Not implemented when IMPUTE2 file is "
                                      "not indexed (see genipe)")

        # Getting the seek position
        if variant_info is None:
            try:
                variant_info = self._impute2_index.loc[name, :]

            except KeyError:
                if name in self.get_duplicated_markers():
                    # The variant is a duplicated one, so we go through all the
                    # variants with the same name and the :dupx suffix
                    return [
                        self.get_variant_by_name(dup_name).pop()
                        for dup_name in self.get_duplicated_markers()[name]
                    ]

                else:
                    # The variant is not in the index
                    logging.variant_name_not_found(name)
                    return []

        # Seeking to the right place in the file
        self._impute2_file.seek(variant_info.seek)

        # Parsing the file
        genotypes = self._parse_impute2_line(self._impute2_file.readline())

        # Fixing the object
        self._fix_genotypes_object(genotypes, variant_info)

        return [genotypes]

    def _fix_genotypes_object(self, genotypes, variant_info):
        """Fixes a genotypes object (variant name, multi-allelic value."""
        # Checking the name (if there were duplications)
        if self.has_index and variant_info.name != genotypes.variant.name:
            if not variant_info.name.startswith(genotypes.variant.name):
                raise ValueError("Index file not synced with IMPUTE2 file")
            genotypes.variant.name = variant_info.name

        # Trying to set multi-allelic information
        if self.has_index and self._index_has_location:
            # Location was in the index, so we can automatically set the
            # multi-allelic state of the genotypes
            genotypes.multiallelic = variant_info.multiallelic

        else:
            # Location was not in the index, so we check one marker before and
            # after the one we found
            logging.warning("Multiallelic variants are not detected on "
                            "unindexed files.")

    def get_number_samples(self):
        """Returns the number of samples.

        Returns:
            int: The number of samples.

        """
        return self.samples.shape[0]

    def get_number_variants(self):
        """Returns the number of markers.

        Returns:
            int: The number of markers.

        """
        if self.has_index:
            return self._impute2_index.shape[0]
        else:
            return None

    def get_samples(self):
        return list(self.samples.index)

    def _parse_impute2_line(self, line):
        """Parses the current IMPUTE2 line (a single variant).

        Args:
            line (str): An IMPUTE2 line.

        Returns:
            Genotypes: The genotype in dosage format.

        Warning
        =======
            By default, the genotypes object has multiallelic set to False.

        """
        # Splitting
        row = line.rstrip("\r\n").split(" ")

        # Constructing the probabilities
        prob = np.array(row[5:], dtype=float)
        prob.shape = (prob.shape[0] // 3, 3)

        # Constructing the dosage
        dosage = 2 * prob[:, 2] + prob[:, 1]
        if self.prob_t > 0:
            dosage[~np.any(prob >= self.prob_t, axis=1)] = np.nan

        return Genotypes(
            Variant(row[1], CHROM_STR_ENCODE.get(row[0], row[0]), int(row[2]),
                    [row[3], row[4]]),
            dosage,
            reference=row[3],
            coded=row[4],
            multiallelic=False,
        )
