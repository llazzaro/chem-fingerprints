Installing
==========

The no-cost version of chemfp is available from `http://pypi.python.org/pypi/chemfp/ PyPI`_
so if you use the third-party tools
`http://www.pip-installer.org/en/latest/index.html pip`_
or 
`http://peak.telecommunity.com/DevCenter/EasyInstall easy_install`_
then you can install it with one of:

  pip install chemfp

  easy_install chemfp

The chemfp tools depends on a working Python installation.  You can
download Python 2.7 from http://www.python.org/download/.  Note that
older versions of OEChem only support up to Python 2.6. At some point
chemfp will support Python 3.3 or later.

The core chemfp functionality does not depend on a third-party library
but you will need a chemistry toolkit in order to generate new
fingerprints from structure files. chemfp supports the free Open Babel
and RDKit toolkits and the proprietary OEChem toolkit. Make sure you
install the Python libraries for the toolkit(s) you select.

If you have a source version of chemfp then you will need a C compiler
in order to compile it. This uses Python's standard "setup.py". Read
http://docs.python.org/install/index.html for details of how to use
it. The short version is that on Unix systems using sudo (that is, Mac
OS X and most Linux-based OSes) you can do::


  sudo python setup.py install

while for Windows you can do::

   python setup.py install

Bear in mind Python 2.7 for Windows was built with Visual
C++ 2008. The setup.py step does not work with Visual C++ 2010 or
later unless you patch your local version of Python. The bug report is
at http://bugs.python.org/issue13210 .
