#!/usr/bin/env python

from distutils.core import setup, Extension

DESCRIPTION  = """\
chemfp is a set of command-lines tools for generating cheminformatics
fingerprints and searching those fingerprints by Tanimoto similarity,
as well as a Python library which can be used to build new tools.

These algorithms are designed for the dense, 100-10,000 bit
fingerprints which occur in small-molecule/pharmaceutical
chemisty. The Tanimoto search algorithms are implemented in C for
performance and support both threshold and k-nearest searches.

Fingerprint generation can be done either by extracting existing
fingerprint data from an SD file or by using an existing chemistry
toolkit. chemfp supports the Python libraries from Open Babel,
OpenEye, and RDKit toolkits.
"""

import sys
if sys.platform[:3] != "win":
    extra_compile_args = ["-O3"]
else:
    extra_compile_args = []

setup(name = "chemfp",
      version = "1.1a2",
      description = DESCRIPTION,
      author = "Andrew Dalke",
      author_email = 'dalke@dalkescientific.com',
      url = "http://code.google.com/p/chem-fingerprints/",
      license = "MIT",
      classifiers = ["Development Status :: 5 - Production/Stable",
                     "Environment :: Console",
                     "License :: OSI Approved :: MIT License",
                     "Operating System :: OS Independent",
                     "Programming Language :: Python",
                     "Topic :: Scientific/Engineering :: Chemistry",
                     "Topic :: Software Development :: Libraries :: Python Modules"],
      
      packages = ["chemfp", "chemfp.commandline", "chemfp.futures", "chemfp.progressbar"],
      package_data = {"chemfp": ["rdmaccs.patterns", "substruct.patterns"]},
      scripts = ["ob2fps", "oe2fps", "rdkit2fps", "sdf2fps", "simsearch"],

      ext_modules = [Extension("_chemfp",
                               ["src/bitops.c", "src/chemfp.c",
                                "src/heapq.c", "src/fps.c",
                                "src/searches.c",
                                "src/select_popcount.c",
                                "src/python_api.c"],
                               extra_compile_args = extra_compile_args,
                               )],
     )
