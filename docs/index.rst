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

Chemfp is regularly tested on a Mac using OEChem 1.7.4, OEChem 1.7.6,
OpenBabel 2.2.3, OpenBabel 2.3.0, OpenBabel 2.3+svn, RDKit 201012,
RDKit 201103, RDKit 201106 and RDKit 201112 as well as Python 2.5,
2.6, and 2.7.


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

Further chemfp development depends on fund from people like you. You
can purchase support or consulting hours from me, and you can pay my
to develop new features. Additionally, I make my living developing
custom software for computational chemistry and giving training
courses for computational chemists who want to get up to speed on
Python programming for their field.  If you are interested in hiring
my services, contact me at dalke@dalkescientific.com .


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
 - chemfp/Watcher.py by Allen Downey under the Python license
 - chemfp/argparse.py by Steven J. Bethard under the Apache License 2.0


Future
======

The chemfp code base is solid and in use at several companies. It has
great support for fingerprint generation and similarity search, but
there's plenty left to do in future.

The biggest question is how to fund the future. Long experience in
this field says that few people are willing to pay a few hundred
dollars for 'free' software, even when they are willing to have
in-house developers work for a week to develop code which is not as
powerful. My current plan is that 1.2 and later will have two release
schedules; one for paying customers, and another a few months later
for those who don't want pay. The software in both cases will be under
the current MIT license.

The most likely near-term chemfp improvements are to revamp the
internal fingerprint types API and make it public, and fix limitations
that prevent chemfp from searching arenas which are larger than 2GB.

There are several internal APIs, like the SDF reader code, which
should be refined and made part of the public API. This includes
validating the reader against a larger number of real-world SD files.

If you are not a Python programmer then you might prefer that the core
search routines be made accessible through a C API. That's possible,
in that the software was designed with that in mind, but it's never
been tested.

I can think of a number of possibly useful command-line tools, like
being able to merge a number of FPS files into one, sorted by popcount
and able to handle all of PubChem.  Or selecting N fingerprints at
random from a set of fingerprints.

The 1.1 release supports OpenMP for multi-core support. I'm looking
for people who will fund me to develop a distributed computing version
to run on your compute cluster.

I'm really interested in using chemfp to handle different sorts of
clustering. Let me know if there are things I can add to the API which
would help you do that.

Are you working with large data sets? Is load time a problem for you?
I have ideas on how to develop a binary format and would like to
develop it futher, but it doesn't make much sense if no one is
interested in it.

There are any number of higher-level tools which can be built on the
chemfp components. For example, what about a wsgi component which
implements a web-based search API for your local network?

There's a paper on using locality-sensitive hashing to find highly
similar fingerprints. Are there cases where it's more useful than
chemfp?

But fundamentally, future work will be guided by what people want, and
by funding. `Let me know <mailto:dalke@dalkescientific.com>`_ in
either case.



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

