import sys
import unittest
import doctest

# I don't know the right way to do this.
# Faking it for now.
class TestEncodings(unittest.TestCase):
    def test_spam(self):
        pass
    pass


for x in doctest.DocTestSuite("chemfp.decoders")._tests:
    name = "test_" + x._dt_test.name.replace(".", "__")
    setattr(TestEncodings, name, x)

print TestEncodings

if __name__ == "__main__":
    unittest.main()
