"""
Define the API for geneparse.
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


import numpy as np

from .exceptions import InvalidChromosome


VALID_CHROMOSOMES = set(
    [str(i + 1) for i in range(22)] + ["X", "Y", "XY", "MT"]
)


class Chromosome(object):
    __slots__ = ("name")

    def __init__(self, name):
        self.name = str(name)

    def __repr__(self):
        return "{}".format(self.name)

    def __hash__(self):
        return hash(self.name)

    def __eq__(self, other):
        other_name = getattr(other, "name", str(other))
        if self.name == "?" or other_name == "?":
            return False
        return self.name == other_name


VALID_CHROMOSOMES = {k: Chromosome(k) for k in VALID_CHROMOSOMES}
UNKNOWN_CHROMOSOME = Chromosome("?")


class Variant(object):
    # Subclasses should declare a __slots__ containing only the additional
    # slots.
    __slots__ = ("name", "chrom", "pos", "alleles")

    def __init__(self, name, chrom, pos, alleles):
        self.name = str(name)
        self.chrom = Variant._encode_chr(chrom)
        self.pos = int(pos)

        self.alleles = self._encode_alleles(alleles)

    @staticmethod
    def _encode_chr(chrom):
        # Accept instances of Chromosome as is (useful for contigs or non
        # default chromosomes).
        if isinstance(chrom, Chromosome):
            return chrom

        # We have a special value for unknown chromosome.
        elif chrom is None:
            return UNKNOWN_CHROMOSOME

        # See if the Chromosome is already known.
        chrom = str(chrom).upper()

        if chrom.startswith("CHR"):
            chrom = chrom[3:]

        if chrom not in VALID_CHROMOSOMES:
            raise InvalidChromosome(chrom)

        return VALID_CHROMOSOMES[chrom]

    @staticmethod
    def _encode_alleles(iterable):
        if iterable is None:
            return None
        return tuple(sorted(str(s).upper() for s in iterable))

    def copy(self):
        return Variant(self.name, self.chrom, self.pos, self.alleles)

    def __hash__(self):
        # Two variants will have the same hash if they have the same
        # chromosome and position and **exactly the same alleles**.
        # Is this the behaviour we want?
        return hash((self.chrom, self.pos, self.alleles))

    def alleles_ambiguous(self):
        return self.alleles == ("C", "G") or self.alleles == ("A", "T")

    @property
    def alleles_set(self):
        if self.alleles is None:
            return None
        return set(self.alleles)

    def primitive_locus_eq(self, chrom, pos):
        return self.chrom == chrom and self.pos == pos

    def locus_eq(self, other):
        return self.primitive_locus_eq(other.chrom, other.pos)

    def iterable_alleles_eq(self, alleles):
        return self.alleles == self._encode_alleles(alleles)

    def alleles_eq(self, other):
        return self.iterable_alleles_eq(other.alleles)

    def __eq__(self, other):
        """Tests for the equality between two variants.

        If any variant has undefined alleles, we return the locus equality.

        Else, we return True if at least two alleles are the same in both
        variants.

        """
        locus_match = self.locus_eq(other)
        if self.alleles is None or other.alleles is None:
            return locus_match

        overlap = len(self.alleles_set & other.alleles_set) >= 2
        return locus_match and overlap

    def __repr__(self):
        return "<Variant chr{}:{}_{}>".format(self.chrom, self.pos,
                                              self.alleles)


class ImputedVariant(Variant):
    __slots__ = ("quality", )

    def __init__(self, name, chrom, pos, alleles, quality):
        super().__init__(name, chrom, pos, alleles)
        self.quality = float(quality)

        if 0 <= self.quality <= 1:
            raise ValueError(
                "The 'quality' field for ImputedVariant instances is expected "
                "to be a float value between 0 and 1."
            )


class Genotypes(object):
    __slots__ = ("variant", "genotypes", "reference", "coded", "multiallelic")

    def __init__(self, variant, genotypes, reference, coded, multiallelic):
        """Class holding information on a variant as well as a vector of
        genotypes.

        The "reference" allele corresponds to 0 and the "coded" allele
        corresponds to 1.

        """
        self.variant = variant
        self.genotypes = genotypes

        self.reference = str(reference).upper()

        self.multiallelic = multiallelic

        if variant.alleles and (self.reference not in variant.alleles):
            raise ValueError(
                "reference allele not in the known alleles for the variant "
                "({} not in {}).".format(self.reference, variant.alleles)
            )

        self.coded = str(coded).upper()

        if variant.alleles and (self.coded not in variant.alleles):
            raise ValueError(
                "coded allele not in the known alleles for the variant "
                "({} not in {}).".format(self.coded, variant.alleles)
            )

    def flip(self):
        """Flips the reference and coded alleles of this instance."""
        self.genotypes = 2 - self.genotypes
        self.reference, self.coded = self.coded, self.reference

    def maf(self):
        freq = self.coded_freq()
        if freq > 0.5:
            return 1 - freq
        else:
            return freq

    def coded_freq(self):
        """Gets the frequency of the coded allele."""
        return np.nanmean(self.genotypes) / 2

    def __eq__(self, other):
        # If not the same locus, not equals.
        if not self.variant.locus_eq(other.variant):
            return False

        # If same alleles, return the genotype comparison as-is.
        alleles_eq = (self.reference == other.reference and
                      self.coded == other.coded)
        if alleles_eq:
            return _np_eq(self.genotypes, other.genotypes)

        # Check if it's the same alleles but they are flipped.
        if self.reference == other.coded and self.coded == other.reference:
            return _np_eq(self.genotypes, (2 - other.genotypes))

        raise RuntimeError("Failed equality check between genotypes.")

    def __repr__(self):
        return (
            "<Genotypes for {} Reference:{} Coded:{}, {}>"
            "".format(self.variant, self.reference, self.coded, self.genotypes)
        )


class SplitChromosomeReader(object):
    def __init__(self, chrom_to_reader):
        """Reader to handle genotype access using files split by chromosome.

        A dict mapping chromosomes to instances of GenotypesReader should be
        passed.

        """
        self.chrom_to_reader = chrom_to_reader

        samples = None
        self.n_vars = 0
        for chrom, reader in self.chrom_to_reader.items():
            # Keep track of the total number of variants.
            self.n_vars += reader.get_number_variants()

            # Check that the sample order is the same.
            cur_samples = reader.get_samples()
            if samples is None:
                samples = cur_samples
            else:
                if samples != cur_samples:
                    raise ValueError(
                        "Not all sub-readers have the same sample order."
                    )
                samples = cur_samples

        self.samples = samples

    @staticmethod
    def _unknown_chrom_message(chrom):
        return (
            "Unable to find a reader instance for chromosome '{}'."
            "".format(chrom)
        )

    def iter_variants(self):
        for chrom, reader in self.chrom_to_reader.items():
            for v in reader.iter_variants():
                yield v

    def iter_genotypes(self):
        for chrom, reader in self.chrom_to_reader.items():
            for g in reader.iter_genotypes():
                yield g

    def get_variant_genotypes(self, variant):
        try:
            return self.chrom_to_reader[
                variant.chrom
            ].get_variant_genotypes(variant)
        except KeyError:
            raise ValueError(self._unknown_chrom_message(variant.chrom))

    def get_variant_by_name(self, name):
        out = []
        for chrom, reader in self.chrom_to_reader.items():
            out.extend(reader.get_variant_by_name(name))
        return out

    def get_variants_in_region(self, chrom, start, end):
        try:
            return self.chrom_to_reader[
                chrom
            ].get_variants_in_region(chrom, start, end)
        except KeyError:
            raise ValueError(self._unknown_chrom_message(chrom))

    def get_samples(self):
        return self.samples

    def get_number_samples(self):
        return len(self.samples)

    def get_number_variants(self):
        return self.n_vars


class GenotypesReader(object):
    def __init__(self):
        """Abstract class to read genotypes data."""
        raise NotImplementedError()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    def close(self):
        pass

    def __repr__(self):
        return "<{} {:,d} samples; {:,d} variants>".format(
            self.__class__.__name__,
            self.get_number_samples(),
            self.get_number_variants(),
        )

    # API methods
    def iter_variants(self):
        """Iterate over variants without reading the actual genotypes.

        This is a generator of Variant instances. Also not that subclasses
        can define their own Variant subclasses to represent additional
        fields.

        To improve consistency, the ImputedVariant class is provided. It
        defines a single additional field ("quality") containing a float
        between 0 and 1.

        """
        raise NotImplementedError()

    def iter_genotypes(self):
        """Iterate over variants and read the genotypes.

        This method yields instances of Genotypes.

        """
        raise NotImplementedError()

    def get_variant_genotypes(self, variant):
        """Get the genotypes for a given variant.

        Args:
            variant (Variant): A variant for which to retrieve genotypes.

        Returns:
            list: A list of Genotypes. This is a list because for
            multi-allelics the representation needs multiple entries.

        """
        raise NotImplementedError()

    def get_variant_by_name(self, name):
        """Get the genotypes for a given variant (by name).

        Args:
            name (str): The name of the variant to retrieve the genotypes.

        Returns:
            list: A list of Genotypes. This is a list in order to keep the same
            behaviour as the other functions.

        """
        raise NotImplementedError()

    def get_variants_in_region(self, chrom, start, end):
        """Get the variants in a region.

        Args:
            chrom (str): The chromosome (e.g. 'X' or '3').
            start (int): The start position for the region.
            end (int): The end position for the region.

        """
        raise NotImplementedError()

    def get_samples(self):
        """Get an ordered collection of the samples in the genotype container.
        """
        raise NotImplementedError()

    def get_number_samples(self):
        """Return the number of samples."""
        raise NotImplementedError()

    def get_number_variants(self):
        """Return the number of variants in the file."""
        raise NotImplementedError()


def _np_eq(a, b):
    nan_a = np.isnan(a)
    nan_b = np.isnan(b)
    if not np.all(nan_a == nan_b):
        return False

    return np.all(a[~nan_a] == b[~nan_b])


# Exceptions.
class UnknownMinorAllele(Exception):
    pass
