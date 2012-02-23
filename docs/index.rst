.. _intro:

====================
chemfp documentation
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

  Copyright (c) 2010-2011 Andrew Dalke Scientific, AB (Gothenburg, Sweden)
  
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


Future
======

Chemfp will progress based on a combination of my interests and your
funding. (Funny how giving me money affects my interests.) Here are
some of the possibilities I've thought of:

The core API documented here is unlikely to change in the near future,
but there are many ways to improve it. Internal routines to read an SD
file and manage fingerprint types should be made public, along with
the API for reading SD files. There will also be a way to read
structure input from a string rather than only from a file.

If you are not a Python programmer then you might prefer that the core
search routines be made accessible through a C API. That's possible,
in that the software was designed with that in mind, but it's never
been tested.

The code makes good use of multiple processors on a shared memory
machine. If you want to cluster a large data set with a high
similarity threshold then you should be able to make effective use of
a distributed compute cluster. That sounds like an excellent project
that you can fund.

If you work on the command-line then you would probably like
command-line tools which can merge a number of FPS files into one,
sort the fingerprints by popcount, or select N fingerprints at random.

The dense binary fingerprints which chemfp focuses on are only a
subset of the cheminformatics fingerprint types. I would really like
to support sparse fingerprints and count fingerprints, along with fast
search code and efficient memory use. I estimate that will take
several months of development, so I really needs external funding for
it.

Another advanced project is to add locality-sensitive hashing, which
would turn the O(N**2) nearest-neighbor algorithm into an O(N)
probabilistic algorithm, for a high enough similarity threshold.

There are any number of higher-level tools which can be built on the
chemfp components. For example, what about a wsgi component which
implements a web-based search API for your local network?

But fundamentally, future work will be guided by what people want, and
by funding. `Let me know <mailto:dalke@dalkescientific.com>`_ in
either case.



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

