> (chemfp for short) is a set of formats and related tools for the storage, exchange, and search of cheminformatics fingerprint data sets.

Download [chemfp-1.1p1.tar.gz](http://code.google.com/p/chem-fingerprints/downloads/detail?name=chemfp-1.1p1.tar.gz). Note: this is a source distribution. Precompiled installers look complicated to do, and will wait until a paying customer asks for them.

The main web site is http://chemfp.com/ . This Google Code site hosts the project information and the public source code and distribution. Releases past 1.1 use a delayed release system. [Paying customers](http://chemfp.com/license/) get access to the latest releases, and users of the no-cost version get access to an older version.

# Cheminformatics Fingerprint #

The chem-fingerprints (chemfp) project goals are to define and promote
common file formats for storing and exchanging cheminformatics
fingerprint data sets, and to develop tools which work with that
format.

The [FPS](FPS.md) format is used for dense binary fingerprints of size less
than about 10,000 bits, and usually only about 1,024 bits. The
UseCases page describes a few examples of how people might use the FPS
format.

People mostly use tools, not file formats. The chemfp distribution
includes five command-line programs for working with FPS files:

  * [sdf2fps](sdf2fps.md) - extract fingerprint data from SD tags
  * [ob2fps](ob2fps.md) - use OpenBabel to generate fingerprints from structures
  * [oe2fps](oe2fps.md) - use OEChem/OEGraphSim to generate fingerprints from structures
  * [rdkit2fps](rdkit2fps.md) - use RDKit to generate fingerprints from structures
  * [simsearch](simsearch.md) - do Tanimoto similarity searches between two FPS files

## Documentation ##

The [documentation](https://readthedocs.org/docs/chemfp/en/latest/) includes descriptions of how to
[use the chemfp command-line tools](https://readthedocs.org/docs/chemfp/en/latest/using-tools.html) to
[extract fingerprints from PubChem data](https://chemfp.readthedocs.org/en/latest/using-tools.html#generating-fingerprint-files-from-pubchem-sd-files)
or
[generate them from ChEBI data](https://chemfp.readthedocs.org/en/latest/using-tools.html#using-a-toolkit-to-process-the-chebi-dataset), and how to carry out
[theshold-based](https://chemfp.readthedocs.org/en/latest/using-tools.html#threshold-search)
and
[k-nearest neighbor](https://chemfp.readthedocs.org/en/latest/using-tools.html#k-nearest-neighbor-search) Tanimoto searches.


For the programmers out there, the command-line tools are built on top
of the "chemfp" Python library, and portions of the library are
available for public use and
[documented](https://chemfp.readthedocs.org/en/latest/api.html). You can calculate
[population counts](https://chemfp.readthedocs.org/en/latest/using-api.html#byte-and-hex-fingerprints) directly, but it's better to use the
[built-in Tanimoto search routines](https://chemfp.readthedocs.org/en/latest/using-api.html#how-to-use-query-fingerprints-to-search-for-similar-target-fingerprints) and
[compute a distance matrix](https://chemfp.readthedocs.org/en/latest/using-api.html#computing-a-distance-matrix-for-clustering) or implement the
[Butina clustering algorithm](https://chemfp.readthedocs.org/en/latest/using-api.html#taylor-butina-clustering).


## Status ##

The project started in early 2010. The first preview release was in May 2011, 1.0 was released on 20 September 2011 and 1.1 was released on 5 February 2013. The current version is [1.1](http://code.google.com/p/chem-fingerprints/downloads/detail?name=chemfp-1.1.tar.gz). The FPS format is stable, the algorithms are well tested, and the public APIs documented. It's ready for you to use.

There are many more things it can do. You can help out. If you had
ideas, comment, or code contribution then join the
[mailing list](http://eight.pairlist.net/mailman/listinfo/chemfp) or
send email to [me](mailto:dalke@dalkescientific.com) directly.

The planned features for chemfp-1.2 are:

  * support searching a fingerprint arena larger than 4GB. (This limits chemfp to 32 million 1024 bit fingerprints)
  * develop a binary format for fast loads of memory-mapped data files
  * improve the fingerprint type system

Some of the other ideas for the future are in the [documentation](https://chemfp.readthedocs.org/en/latest/#future).

## Advertising ##
This program was developed by Andrew Dalke of Andrew Dalke Scientific, AB. It is distributed free of charge under the “MIT” license, shown below.

Further chemfp development depends on funding from people like you. Asking for voluntary contributions almost never works. Instead, starting with chemfp-1.1, the source code is distributed under an incentive program. You can pay for the commerical distribution, or use the no-cost version.

If you pay for the commercial distribution then you will get the most recent version of chemfp, free upgrades for one year, support, and a discount on renewing participation in the incentive program.

If you use the no-cost distribution then you will get the 1.1 version of chemfp, limited support, and minor bug fixes and improvements.

The current plan is that older versions of the commercial distribution will be released under the no-cost program. However, the no-cost version will always be at least one, and more likely two to three years behind the version available to those who fund chemfp development.


If you have questions about or with to purchase the commercial distribution, send an email to sales@dalkescientific.com .