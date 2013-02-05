.. _intro:

====================
chemfp 1.1 documentation
====================

`chemfp <http://code.google.com/p/chem-fingerprints/>`_ is a set of
tools for working with cheminformatics fingerprints in the FPS format.

Most people will use the command-line programs to generate and search
fingerprint files. :ref:`ob2fps <ob2fps>`, :ref:`oe2fps <oe2fps>`, and
:ref:`rdkit2fps <rdkit2fps>` use respectively the `Open Babel
<http://openbabel.org/>`_, `OpenEye <http://www.eyesopen.com/>`_, and
`RDKit <http://www.rdkit.org/>`_ chemistry toolkits to convert
structure files into fingerprint files. :ref:`sdf2fps <sdf2fps>`
extracts fingerprints encoded in SD tags to make the fingerprint
file. :ref:`simsearch <simsearch>` finds targets in a fingerprint file
which are sufficiently similar to the queries.

The programs are built using the :ref:`chemfp Python library API <chemfp-api>`,
which in turn uses a C extension for the performance
critical sections. The parts of the library API documented here are
meant for public use, along with some examples.


Remember: chemfp cannot generate fingerprints from a structure file
without a third-party chemistry toolkit.

Chemfp is regularly tested on a Mac using multiple versions of
OEChem/OEGraphSim, Open Babel, and RDKit as well as Python 2.5, 2.6,
and 2.7.


.. toctree::
   :maxdepth: 2

   installing
   using-tools
   tool-help
   using-api
   api
   

License and advertisement
=========================

This program was developed by Andrew Dalke of Andrew Dalke Scientific,
AB. It is distributed free of charge under the "MIT" license, shown
below.

Further chemfp development depends on funding from people like
you. Asking for voluntary contributions almost never works. Instead,
starting with chemfp-1.1, the source code is distributed under an
incentive program. You can pay for the commerical distribution, or use
the no-cost version.

If you pay for the commercial distribution then you will get the most
recent version of chemfp, free upgrades for one year, support, and a
discount on renewing participation in the incentive program.

If you use the no-cost distribution then you will get the 1.1 version
of chemfp, limited support, and minor bug fixes and improvements.

The current plan is that older versions of the commercial distribution
will be released under the no-cost program. However, the no-cost
version will always be at least one, and more likely two to three
years behind the version available to those who fund chemfp development.

If you have questions about or with to purchase the commercial
distribution, send an email to dalke@dalkescientific.com .


.. highlight:: none

::

  Copyright (c) 2010-2012 Andrew Dalke Scientific, AB (Gothenburg, Sweden)
  
  Permission is hereby granted, free of charge, to any person obtaining
  a copy of this software and associated documentation files (the
  "Software"), to deal in the Software without restriction, including
  without limitation the rights to use, copy, modify, merge, publish,
  distribute, sublicense, and/or sell copies of the Software, and to
  permit persons to whom the Software is furnished to do so, subject to
  the following conditions:
  
  The above copyright notice and this permission notice shall be
  included in all copies or substantial portions of the Software.
  
  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
  MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
  LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
  OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
  WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

Copyright to portions of the code are held by other people or
organizations, and may be under a different license. See the specific
code for details. These are:

 - OpenMP, cpuid, POPCNT, and Lauradoux implementations by Kim
   Walisch, <kim.walisch@gmail.com>, under the MIT license
 - SSSE3.2 popcount implementation by Stanford Univeristy (written by
   Imran S. Haque <ihaque@cs.stanford.edu>) under the BSD license
 - heapq by the Python Software Foundation under the Python license
 - chemfp/progressbar/ by Nilton Volpato under the LGPL 2.1 and/or BSD license
 - chemfp/futures/ by Brian Quinlan under the Python license
 - chemfp/argparse.py by Steven J. Bethard under the Apache License 2.0


Future
======

The chemfp code base is solid and in use at several companies. It has
great support for fingerprint generation and similarity search and
is very fast, but there's plenty left to do in future.

The most likely near-term chemfp improvements are to fix the
limitations that prevent chemfp from searching arenas which are larger
than 2GB and add the FPB binary file format to improve the fingerprint
loading times for large data sets.

The fingerprint type internals are confusing. I have not documented
them as part of the public API because it will be rewritten. Once
done, you will have access to toolkit-specific functions like: list
the supported fingerprint types, create a fingerprint when the
structure is a string, and test if two fingerprint types use the same
toolkit.

I do not have a Microsoft Windows computer. I looked into how I might
provide an installer for that OS. It's a lot of work, if only because
I need to support Python 2.5, 2.6, and 2.7 and learn how installers
work on Windows. I've decided to put it off until a member of the
incentive program specifically wants it.


The goals beyond that will depend in part on user feedback. The
following are possible ideas, to help inspire you. `Let me know
<mailto:dalke@dalkescientific.com>`_ if you need something like one of
these.

The threshold and k-nearest arena search results store hits using
compressed sparse rows. These work well for sparse results, but when
you want the entire similarity matrix (ie, with a minimum threshold of
0.0) of a large arena, then time and space to maintain the sparse data
structure becomes noticable. It's likely in that case that you want
to store the scores in a 2D NumPy matrix.

I'm really interested in using chemfp to handle different sorts of
clustering. Let me know if there are things I can add to the API which
would help you do that.

The Tanimoto the most common similarity method in cheminformatics, but
not the only one. I could add support for Tversky and other
methods. Some of these can also take benefit from chemfp's sublinear
search method.

I think that an awk-like, or LINQ-like command-line selection tool
would be useful. It could be used to select fingerprint records by
position, id, or popcount, or just pick records at random. It's
tempting to thinking of a full-blown language, but that's what Python
is for.

There are several internal APIs, like the SDF reader code, which
should be refined and made part of the public API. This includes
validating the reader against a larger number of real-world SD files.

If you are not a Python programmer then you might prefer that the core
search routines be made accessible through a C API. That's possible,
in that the software was designed with that in mind, but it needs more
development and testing.

ChemFP 1.1 supports OpenMP. That's great for shared-memory
machines. Are you interested in supporting a distributed computing
version?

There are any number of higher-level tools which can be built on the
chemfp components. For example, what about a wsgi component which
implements a web-based search API for your local network?

There's a paper on using locality-sensitive hashing to find highly
similar fingerprints. Are there cases where it's more useful than
chemfp?

I will support Python 3, but so far no one has asked for it. While the
toolkit interfaces will have to wait until the respective toolkits
support Python 3, the similarity code does not depend on any toolkit
and could be ported on its own.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

