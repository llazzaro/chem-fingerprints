import sys
import unittest

import support

from chemfp.commandline import rdkit2fps

runner = support.Runner(rdkit2fps.main)
        
class TestMACCS(unittest.TestCase):
    def test_bitorder(self):
        result = runner.run_fps("--maccs166", 7, "maccs.smi")
        # The fingerprints are constructed to test the first few bytes.
        self.assertEquals(result[0][:6], support.set_bit(2))
        self.assertEquals(result[1][:6], support.set_bit(3))
        self.assertEquals(result[2][:6], support.set_bit(4))
        self.assertEquals(result[3][:6], support.set_bit(5))
        self.assertEquals(result[4][:6], support.set_bit(9))
        self.assertEquals(result[5][:6], support.set_bit(10))
        self.assertEquals(result[6][:6], support.set_bit(16))

if __name__ == "__main__":
    unittest.main()
