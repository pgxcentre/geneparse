"""Facilitates logging for geneparse."""

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
