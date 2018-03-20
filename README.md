[![PyPI version](https://badge.fury.io/py/geneparse.svg)](http://badge.fury.io/py/geneparse)
[![Build Status](https://travis-ci.org/pgxcentre/geneparse.svg?branch=master)](https://travis-ci.org/pgxcentre/geneparse)


# geneparse - Module to parse genetics data file

`geneparse` is a module that helps developers to parse multiple genetics file
format (*e.g.* Plink binary files, IMPUTE2 files and VCF).


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

Using `pip`:

```bash
pip install geneparse
```


### Updating

To update the module using `pip`:

```bash
pip install -U geneparse
```


## Testing

To test the module, just perform the following command:

```console
$ python -m geneparse.tests
.sssss.........s.sssss.........ssssssss...ss.ss...s...................
..........................................
----------------------------------------------------------------------
Ran 112 tests in 0.713s

OK (skipped=24)
```
