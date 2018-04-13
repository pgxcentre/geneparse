[![PyPI version](https://badge.fury.io/py/geneparse.svg)](http://badge.fury.io/py/geneparse)
[![Build Status](https://travis-ci.org/pgxcentre/geneparse.svg?branch=master)](https://travis-ci.org/pgxcentre/geneparse)


# geneparse - Module to parse genetics data file

`geneparse` is a module that helps developers to parse multiple genetics file
format (*e.g.* Plink binary files, IMPUTE2 files, BGEN and VCF).


## Dependencies

The tool requires a standard [Python](http://python.org/) installation (3.4 or
higher are supported) with the following modules:

1. [numpy](http://www.numpy.org/)
2. [pandas](http://pandas.pydata.org/)
3. [pyplink](https://github.com/lemieuxl/pyplink)
4. [pybgen](https://github.com/lemieuxl/pybgen)
5. [cyvcf2](https://github.com/brentp/cyvcf2)
6. [biopython](https://github.com/biopython/biopython)

The tool has been tested on *Linux* only, but should work on *MacOS* operating
systems as well.


## Installation

You can install or update `geneparse` using `pip`:

```bash
pip install -U geneparse
```


## Testing

To test the module, just perform the following command:

```console
$ python -m geneparse.tests
.sssss.........s.sssss.........ssssssss...ss.ss...s...................
.....................................................s....ss.....
----------------------------------------------------------------------
Ran 135 tests in 1.064s

OK (skipped=27)
```


## Indexing

Some genotype data require indexing for fast access. This can be done using
geneparse.

```console
$ python -m geneparse.index --help
usage: geneparse-indexer [-h] [--impute2 IMPUTE2 [IMPUTE2 ...]]
                         [--bgen BGEN [BGEN ...]] [--legacy]

Genotype file indexer.

optional arguments:
  -h, --help            show this help message and exit

IMPUTE2 index:
  --impute2 IMPUTE2 [IMPUTE2 ...]
                        Index an IMPUTE2 genotype file format. The file can be
                        plain text or bgzipped.

BGEN index:
  --bgen BGEN [BGEN ...]
                        Index a BGEN genotype file. This requires 'bgenix' to
                        be in the PATH.
  --legacy              Index the file using the '-with-rowid' option. This
                        flag enables compatibility with SQLITE prior to
                        version 3.8.2. See
                        https://bitbucket.org/gavinband/bgen/wiki/bgenix for
                        more information.
```


## Extraction

We provide a simple tool to extract genotypes from different format to either
VCF or Binary plink files.


```console
$ python -m geneparse.extract --help
usage: geneparse-extractor [-h] -f FORMAT [-e FILE] [-k FILE] -o FILE
                           [--output-format FORMAT]
                           PARSER_ARGS [PARSER_ARGS ...]

Genotype file extractor. This tool will extract markers according to names or
to genomic locations.

optional arguments:
  -h, --help            show this help message and exit

Input Options:
  -f FORMAT, --format FORMAT
                        The input file format.
  PARSER_ARGS           The arguments that will be passed to the genotype
                        parsers.

Extract Options:
  -e FILE, --extract FILE
                        The list of markers to extract (one per line, no
                        header).
  -k FILE, --keep FILE  The list of samples to keep (one per line, no header).

Output Options:
  -o FILE, --output FILE
                        The output file (can be '-' for STDOUT when using VCF
                        as output).
  --output-format FORMAT
                        The output file format. Note that the extension will
                        be added if absent.

The parser arguments (PARSER_ARGS) are the same as the one in the API. For
example, the arguments for the Plink parser is 'prefix:PREFIX' (where PREFIX
is the prefix of the BED/BIM/FAM files).
```
