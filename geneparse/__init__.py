from . import plink
from .core import Genotypes, Variant, ImputedVariant

parsers = {
    "plink": plink.PlinkReader,
}
