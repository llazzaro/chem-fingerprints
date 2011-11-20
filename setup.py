#!/usr/bin/env python

from distutils.core import setup, Extension
from distutils.command.build_ext import build_ext

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

USE_OPENMP = False

# chemfp has experimental support for OpenMP, but only for the counts.
# Preliminary tests 
def OMP(*args):
    if USE_OPENMP:
        return list(args)
    return []


# Compiler-specific configuration settings due to:
#  http://stackoverflow.com/questions/724664/python-distutils-how-to-get-a-compiler-that-is-going-to-be-used/5192738#5192738
copt =  {
    "msvc": OMP("/openmp") + ["/Ox"],
    "mingw32" : OMP("-fopenmp") + ["-O3", "-ffast-math", "-march=native"],

    "gcc-4.1": ["-O3"], # Doesn't support OpenMP
    "gcc": OMP("-fopenmp") + ["-O3"],
    }

lopt =  {
    "mingw32" : OMP("-fopenmp"),

    "gcc-4.1": ["-O3"], # Doesn't support OpenMP
    "gcc": OMP("-fopenmp") + ["-O3"],
    }


def _is_gcc(compiler):
    return "gcc" in compiler or "g++" in compiler

class build_ext_subclass( build_ext ):
    def build_extensions(self):
        c = self.compiler.compiler_type
        if c == "unix":
            c = self.compiler.compiler[0]
            if _is_gcc(c):
                names = [c, "gcc"]
            else:
                names = [c]
        else:
            names = [c]

        for c in names:
            if copt.has_key(c):
                for e in self.extensions:
                    e.extra_compile_args = copt[ c ]
                break

        for c in names:
            if lopt.has_key(c):
                for e in self.extensions:
                    e.extra_link_args = lopt[ c ]
                break
        
        build_ext.build_extensions(self)




setup(name = "chemfp",
      version = "1.1a3",
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
                                "src/lauradoux.c", "src/lut16.c", "src/gillies.c",
                                "src/python_api.c"],
                               )],
      cmdclass = {"build_ext": build_ext_subclass},
     )
