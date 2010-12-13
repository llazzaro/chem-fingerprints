import os
import shutil
import sys
import tempfile
import unittest2

from chemfp.commandline import ob2fps

import support

real_stdout = sys.stdout

runner = support.Runner(ob2fps.main)
run = runner.run
run_fps = runner.run_fps
run_split = runner.run_split

class TestFingerprintTypes(unittest2.TestCase):
    def test_FP2(self):
        headers, fps = run_split("--FP2", 19)
        self.assertEquals(headers["#type"], "OpenBabel-FP2/1")
        self.assertEquals(fps[0], "200206000000000402800e00040140010100014008206200000200c0082200080200500201c9804100270538000000402000a2040080c1240001c2c2004600200c200c04020800200410a0001490000200a803c018005400c80c00000000810100840000880064a0124010000000080102060142400110200a00000004800000 9425004")
    def test_FP3(self):
        headers, fps = run_split("--FP3", 19)
        self.assertEquals(headers["#type"], "OpenBabel-FP3/1")
        self.assertEquals(fps[0], "0400000000b001 9425004")
    def test_FP4(self):
        headers, fps = run_split("--FP4", 19)
        self.assertEquals(headers["#type"], "OpenBabel-FP4/1")
        self.assertEquals(fps[0], "1100000000000000000080000000000000010000000c9800000000000000000000000640407800 9425004")
    def test_MACCS(self):
        headers, fps = run_split("--MACCS", 19)
        self.assertEquals(headers["#type"], "OpenBabel-MACCS/1")
        self.assertEquals(fps[0], "800400000002080019cc40eacdec980baea378ef1b 9425004")

class TestIO(unittest2.TestCase):
    def test_compressed_auto(self):
        header, fps = run_split("--FP3", 19, support.PUBCHEM_SDF_GZ)
        self.assertEquals(fps[0], "0400000000b001 9425004")
    def test_compressed_specified(self):
        header, fps = run_split("--FP3 --in sdf.gz", 19, support.PUBCHEM_SDF_GZ)
        self.assertEquals(fps[0], "0400000000b001 9425004")
    def test_format_specified(self):
        header, fps = run_split("--FP3 --in sdf", 19, support.PUBCHEM_ANOTHER_EXT)
        self.assertEquals(fps[0], "0400000000b001 9425004")
    def test_output(self):
        dirname = tempfile.mkdtemp(prefix="test_ob2fps")
        output_filename = os.path.join(dirname, "blah.fps")
        assert len(output_filename.split()) == 1 # ensure no whitespace
        try:
            output = runner.run("--FP3 -o " + output_filename)
            assert len(output) == 0
            with open(output_filename) as f:
                result = f.readlines()
        finally:
            shutil.rmtree(dirname)
        self.assertEquals(result[0], "#FPS1\n")
        fps = [line for line in result if not line.startswith("#")]
        self.assertEquals(len(fps), 19)
        self.assertEquals(fps[0], "0400000000b001 9425004\n")
        

class TestMACCS(unittest2.TestCase):
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
    unittest2.main()
