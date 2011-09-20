.. _intro:

====================
chemfp documentation
====================

chemfp is a set of tools for working with cheminformatics fingerprints
in the FPS format.

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
AB. I make my living doing custom software development for
cheminformatics and related fields, and giving training courses for
computational chemists who want to get up to speed on Python
programming for their field.

If you are interested in hiring my services, contact me at
dalke@dalkescientific.com .

Despite that I like having money, this software is available under
what's often called "the MIT license."

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

There's still a lot to do before this package can be called "mature."
For one, the popcount implementation is an 8-bit lookup table. A
16-bit table is faster on most machines, and if you have an Intel
machine then there's a couple of assembly-based techniques which
should almost triple the performance.

The SDF reader code should be made part of the public API, which means
validating against a larger number of real-world SD files.

There should be a way to generate structure fingerprints when the
input is a string, rather than a file.

There should be direct support for building the upper-triangle
distance matrix used for clustering.

I can think of a number of possibly useful command-line tools, like
being able to merge a number of FPS files into one, sorted by popcount
and able to handle all of PubChem.  Or selecting N fingerprints at
random from a set of fingerprints.

I've discussed and proposed a first draft of a binary format, which
would improve load time. I would like to develop it further, but right
I don't know who is bound by load time performance.

I want to develop multi-threaded versions of the search functions. The
core search algorithms are thread-safe, and the Python code releases
the global interpreter lock before going into them, so this shouldn't
be all that hard. (Famous last words.)

What about a wsgi component which implements a web-based search API
for your local network?

Do you want the underlying C code which does the searches to be
available as a C API so you can call it from non-Python programs?

But fundamentally, future work will be guided by what people want, and
by funding. `Let me know <mailto:dalke@dalkescientific.com>`_ in
either case.



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

