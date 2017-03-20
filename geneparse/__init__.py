from . import plink
from .core import Genotypes, Variant

parsers = {
    "plink": plink.PlinkReader,
}
