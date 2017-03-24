from . import plink, impute2
from .core import Genotypes, Variant, ImputedVariant

parsers = {
    "plink":   plink.PlinkReader,
    "impute2": impute2.Impute2Reader,
}
