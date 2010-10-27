# I don't understand the combination of unittest and doctest well
# enough. This feels like a kludge, and it doesn't support
# nosetest.

import unittest

import test_docstrings
import test_chemfp
import test_sdf_reader

suite = unittest.TestSuite()
suite.addTest(test_docstrings.suite())
suite.addTest(unittest.defaultTestLoader.loadTestsFromModule(test_chemfp))
suite.addTest(unittest.defaultTestLoader.loadTestsFromModule(test_sdf_reader))

class MyLoader(object):
    def loadTestsFromModule(self, name):
        return suite

unittest.main(testLoader = MyLoader())
