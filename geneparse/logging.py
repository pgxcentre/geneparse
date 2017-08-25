import logging

from .core import Variant
from . import config

_logger = logging.getLogger("geneparse")


def info(message):
    _logger.info(message)


def warning(message):
    _logger.warning(message)


# Standardized messages.
def found_duplicates(counts):
    """Log that duplicates were found.

    :param counts: A list of duplicate marker names along with their number
                   of occurences.
    :type counts: list

    """
    _logger.warning("Duplicated markers found")
    for marker, count in counts:
        _logger.warning(" - {}: {:,d} times".format(marker, count))
    _logger.warning("Appending ':dupX' to the duplicated markers according "
                    "to their location in the file.")


def variant_not_found(v):
    if not config.LOG_NOT_FOUND:
        return

    if not isinstance(v, Variant):
        raise ValueError("Expected a variant.")

    _logger.warning("Variant {} was not found.".format(v))


def variant_name_not_found(name):
    if not config.LOG_NOT_FOUND:
        return

    _logger.warning("Variant named {} was not found.".format(name))
