import sys
import unittest
import doctest

# XXX I do not know how to make doctest work with nosetests
# Even this code here feels like a hack - I want the '-v' support under unittest

class MyLoader(object):
    def loadTestsFromModule(self, name):
        return suite()

def suite():
    return doctest.DocTestSuite("chemfp.decoders")

if __name__ == "__main__":
    unittest.main(testLoader = MyLoader())
