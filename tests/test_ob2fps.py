import unittest
import sys

from chemfp.commandline import ob2fps

import support

real_stdout = sys.stdout

runner = support.Runner(ob2fps.main)
run = runner.run
run_fps = runner.run_fps
run_split = runner.run_split

class TestFingerprintTypes(unittest.TestCase):
    def test_FP2(self):
        headers, fps = run_split("--FP2", 19)
        self.assertEquals(headers["#type"], "OpenBabel-FP2/1")
    def test_FP3(self):
        headers, fps = run_split("--FP3", 19)
        self.assertEquals(headers["#type"], "OpenBabel-FP3/1")
    def test_FP4(self):
        headers, fps = run_split("--FP4", 19)
        self.assertEquals(headers["#type"], "OpenBabel-FP4/1")
    def test_MACCS(self):
        headers, fps = run_split("--MACCS", 19)
        self.assertEquals(headers["#type"], "OpenBabel-MACCS/1")

class TestIO(unittest.TestCase):
    pass

class TestMACCS(unittest.TestCase):
    def test_bitorder(self):
        result = runner.run_fps("--MACCS", 7, "maccs.smi")
        # The fingerprints are constructed to test the first few bytes.
        self.assertEquals(result[0][:6], support.set_bit(2))
        self.assertEquals(result[1][:6], support.set_bit(3))
        self.assertEquals(result[2][:6], support.set_bit(4))
        self.assertEquals(result[3][:6], support.set_bit(5))
        self.assertEquals(result[4][:6], support.set_bit(9))
        ## This appears to be a bug in the OpenBabel MACCS definition
        #self.assertEquals(result[5][:6], support.set_bit(10))
        # This is WRONG, since OB has an off-by-one error in the ring sizes
        # Once this is fixed you must update the MACCS keys version number
        self.assertEquals(result[5][:6], "000020")
        self.assertEquals(result[6][:6], support.set_bit(16))

if __name__ == "__main__":
    unittest.main()
