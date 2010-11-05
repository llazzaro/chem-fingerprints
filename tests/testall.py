# I don't understand the combination of unittest and doctest well
# enough. This feels like a kludge, and it doesn't support
# nosetest.

import unittest

import test_docstrings

suite = unittest.TestSuite()
suite.addTest(test_docstrings.suite())

for name in ("test_chemfp", "test_sdf_reader", "test_sdf2fps", "test_oe2fps",
             "test_rdkit2fps", "test_rdkit", "test_ob2fps"):
    m = __import__(name)
    suite.addTest(unittest.defaultTestLoader.loadTestsFromModule(m))


class MyLoader(object):
    def loadTestsFromModule(self, name):
        return suite

unittest.main(testLoader = MyLoader())
