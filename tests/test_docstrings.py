import sys
import unittest2
import doctest

# XXX I do not know how to make doctest work with nosetests
# Even this code here feels like a hack - I want the '-v' support under unittest

#class MyLoader(object):
#    def loadTestsFromModule(self, name):
#        return suite()

def load_tests(loader, tests, ignore):
    tests.addTests(doctest.DocTestSuite("chemfp.decoders"))
    return tests

if __name__ == "__main__":
    unittest2.main()
