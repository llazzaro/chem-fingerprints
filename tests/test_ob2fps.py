import unittest
import sys

from chemfp.commandline import ob2fps

import support

real_stdout = sys.stdout

runner = support.Runner(ob2fps.main)
run = runner.run
run_fps = runner.run_fps

class TestFingerprintTypes(unittest.TestCase):
    def test_FP2(self):
        result = run_fps("--FP2")
        self.assertEquals(len(result), 19)

if __name__ == "__main__":
    unittest.main()
