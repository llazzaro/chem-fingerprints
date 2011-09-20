Installing
==========

Download the chemfp source code from the
`project home page <http://code.google.com/p/chem-fingerprints/`_.

The chemfp tools depends on a working Python installation.  You can
download Python 2.7 from http://www.python.org/download/. (Note for
OEChem users: OpenEye doesn't yet support Python 2.7 so you will need
to install Python 2.6 from
http://www.python.org/download/releases/2.6.7/ .)

The core chemfp functionality does not depend on a third-party library
but you will need a chemistry toolkit in order to generate new
fingerprints from structure files. chemfp supports the free Open Babel
and RDKit toolkits and the proprietary OEChem toolkit. Make sure you
install the Python libraries for the toolkit(s) you select.

If you have a source version of chemfp then you will need a C compiler
in order to compile it. This uses Python's standard "setup.py" and you
can see http://docs.python.org/install/index.html for details of how
to use it. The short version is that on Unix systems using sudo (that
is, Mac OS X and most Linux-based OSes) you can do::


  sudo python setup.py install

while for Windows you can do::

   python setup.py install

