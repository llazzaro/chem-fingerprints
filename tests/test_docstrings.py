import sys
import unittest2
import doctest

def load_tests(loader, tests, ignore):
    tests.addTests(doctest.DocTestSuite("chemfp.decoders"))
    return tests

if __name__ == "__main__":
    unittest2.main()
