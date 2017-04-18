import glob
import re

from . import plink, impute2
from .core import Genotypes, Variant, ImputedVariant, SplitChromosomeReader


# TODO:
# 1. Warn and show last exception if no reader correctly initialized.
# 2. Could also make it async to load faster.
class _SplitChromosomeReaderFactory(object):
    def __init__(self, reader_class):
        self.reader_class = reader_class

    def __call__(self, pattern, *args, **kwargs):
        if "{chrom}" not in pattern:
            raise ValueError("Expected '{chrom}' as a placeholder in the "
                             "pattern.")

        # Explode the path for every possible chromosome.
        chrom_to_reader = {}
        for chrom in list(range(1, 23)) + ["X", "Y", "XY", "MT"]:
            chrom = str(chrom)
            cur = re.sub("{chrom}", chrom, pattern)
            try:
                # Instantiate the reader.
                chrom_to_reader[chrom] = self.reader_class(
                    cur, *args, **kwargs
                )
            except:
                pass

        return SplitChromosomeReader(chrom_to_reader)


parsers = {
    "plink": plink.PlinkReader,
    "chrom-split-plink": _SplitChromosomeReaderFactory(plink.PlinkReader),
    "impute2": impute2.Impute2Reader,
    "chrom-split-impute2": _SplitChromosomeReaderFactory(
        impute2.Impute2Reader
    ),
}
