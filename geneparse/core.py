"""
Define the API for geneparse.
"""


class Variant(object):
    # Subclasses should declare a __slots__ containing only the additional
    # slots.
    __slots__ = ("name", "chrom", "pos", "alleles")

    def __init__(self, name, chrom, pos, alleles):
        self.name = str(name)
        self.chrom = Variant._encode_chr(chrom)
        self.pos = int(pos)

        if alleles is None:
            self.alleles = None
        else:
            self.alleles = tuple(sorted(str(s).upper() for s in alleles))

    @staticmethod
    def _encode_chr(chrom):
        chrom = str(chrom).upper()
        if chrom.startswith("CHR"):
            chrom = chrom[3:]

        return chrom

    # Hashing and comparison.
    def __hash__(self):
        # Two variants will have the same hash if they have the same
        # chromosome and position and **exactly the same alleles**.
        # Is this the behaviour we want?
        return hash((self.chrom, self.pos, self.alleles))

    def locus_eq(self, other):
        return (self.chrom == other.chrom and
                self.pos == other.pos)

    def alleles_subset(self, other):
        if self.alleles is None or other.alleles is None:
            return False

        return len(set(self.alleles) & set(other.alleles)) >= 1

    def __eq__(self, other):
        return self.locus_eq(other) and self.alleles_subset(other)

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
    __slots__ = ("variant", "genotypes", "reference", "coded")

    def __init__(self, variant, genotypes, reference, coded):
        """Class holding information on a variant as well as a vector of
        genotypes.

        The "reference" allele corresponds to 0 and the "coded" allele
        corresponds to 1.

        """
        self.variant = variant
        self.genotypes = genotypes

        self.reference = str(reference).upper()

        if self.reference not in variant.alleles:
            raise ValueError(
                "reference allele not in the known alleles for the variant."
            )

        self.coded = str(coded).upper()

        if self.coded not in variant.alleles:
            raise ValueError(
                "coded allele not in the known alleles for the variant."
            )


class GenotypesReader(object):
    def __init__(self):
        """Abstract class to read genotypes data."""
        raise NotImplementedError()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    def close():
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


# Exceptions.
class UnknownMinorAllele(Exception):
    pass
