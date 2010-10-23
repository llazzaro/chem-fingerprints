import unittest
import re

import chemfp
import _chemfp

version_pattern = re.compile(r"\d+\.\d+(\.\d)?((a|b|pre)\d+)?$")

class SystemTests(unittest.TestCase):
    def test_version(self):
        m = version_pattern.match(_chemfp.version())
        self.assertNotEqual(m, None, "bad version: %s" % (_chemfp.version(),))

if __name__ == "__main__":
    unittest.main()
