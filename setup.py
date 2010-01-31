#!/usr/bin/env python

from distutils.core import setup

setup(name = "chemfp",
      version = "0.5",
      description = "Cheminformatics fingerprint tools",
      author = "Andrew Dalke",
      author_email = 'dalke@dalkescientific.com',
      url = "http://code.google.com/p/chem-fingerprints/",
      packages = ["chemfp"],
      scripts = ["ob2fps", "oe2fps", "rdkit2fps", "sdf2fps"],
     )

