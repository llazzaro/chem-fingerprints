import sys
import unittest
import doctest

# I don't know the right way to do this.
# Faking it for now. It works under unittest but not nosetest.
class TestEncodings(unittest.TestCase):
    pass


for x in doctest.DocTestSuite("chemfp.decoders")._tests:
    name = "test_" + x._dt_test.name.replace(".", "__")
    setattr(TestEncodings, name, x)

if __name__ == "__main__":
    unittest.main()
