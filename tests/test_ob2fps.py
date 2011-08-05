from __future__ import with_statement
import os
import shutil
import tempfile
import unittest2

import support

try:
    import openbabel
    has_openbabel = True
    skip_openbabel = False
except ImportError:
    has_openbabel = False
    skip_openbabel = True

    if not support.can_skip("ob"):
        skip_openbabel = False
        import openbabel
        
if has_openbabel:
    import chemfp.openbabel
    from chemfp.commandline import ob2fps
    VERSION = chemfp.openbabel._ob_version

    runner = support.Runner(ob2fps.main)
    run = runner.run
    run_fps = runner.run_fps
    run_split = runner.run_split
    run_exit = runner.run_exit

    HAS_MACCS = chemfp.openbabel.HAS_MACCS
else:
    HAS_MACCS = False

class TestFingerprintTypes(unittest2.TestCase):
    def test_unspecified(self):
        # Should give the same results as FP2
        headers, fps = run_split("", 19)
        self.assertEquals(headers["#type"], "OpenBabel-FP2/1")
        self.assertEquals(fps[0], "200206000000000402800e00040140010100014008206200000200c0082200080200500201c9804100270538000000402000a2040080c1240001c2c2004600200c200c04020800200410a0001490000200a803c018005400c80c00000000810100840000880064a0124010000000080102060142400110200a00000004800000 9425004")
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
        # Sigh. OpenBabel post 2.3 added stereo support. This structure is
        #   Clc1c(/C=C/C(=O)NNC(=O)Cn2nc(cc2C)C)c(F)ccc1\t9425004\n
        # FP4 bits 289 and 290 are "on" for post 2.3.0 releases
        # Those bits are:
        #  289 Cis_double_bond: */[D2]=[D2]\*
        #  290 Trans_double_bond: */[D2]=[D2]/*
        if (VERSION.startswith("2.2") or
            VERSION == "2.3.0"):
            self.assertEquals(fps[0], "1100000000000000000080000000000000010000000c9800000000000000000000000640407800 9425004")  # old (without stereo)
        else:
            self.assertEquals(fps[0], "1100000000000000000080000000000000010000000c9800000000000000000000000640437800 9425004")  # new (with stereo)

    @unittest2.skipUnless(HAS_MACCS, "Missing MACCS support")
    def test_MACCS(self):
        headers, fps = run_split("--MACCS", 19)
        if headers["#type"] == "OpenBabel-MACCS/1":
            # Running on a buggy 2.3.0 release or earlier
            self.assertEquals(fps[0], "800400000002080019cc40eacdec980baea378ef1b 9425004")
        else:
            self.assertEquals(headers["#type"], "OpenBabel-MACCS/2")
            # Running on a corrected post-2.3.0 release
            self.assertEquals(fps[0], "000000000002080019c444eacd6c980baea178ef1f 9425004")
            
    @unittest2.skipIf(HAS_MACCS, "check for missing MACCS support")
    def test_MACCS_does_not_exist(self):
        run_exit("--MACCS")

TestFingerprintTypes = unittest2.skipIf(skip_openbabel, "OpenBabel not installed")(
    TestFingerprintTypes)


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
        
    def test_missing_filename(self):
        errmsg = run_exit("--FP2", "does_not_exist.smi")
        self.assertIn("Cannot read structures", errmsg)
        self.assertIn("No such file", errmsg)
        self.assertIn("does_not_exist.smi", errmsg)
        
    def test_bad_extension(self):
        errmsg = run_exit("--FP2 --in xyzzy")
        self.assertIn("Unknown structure format", errmsg)
        self.assertIn("xyzzy", errmsg)
        
TestIO = unittest2.skipIf(skip_openbabel, "OpenBabel not installed")(TestIO)

class TestMACCS(unittest2.TestCase):
    @unittest2.skipIf(not HAS_MACCS, "need MACCS support")
    def test_bitorder(self):
        result = runner.run_fps("--MACCS", 7, support.fullpath("maccs.smi"))
        # The fingerprints are constructed to test the first few bytes.
        self.assertEquals(result[0][:6], support.set_bit(2))
        self.assertEquals(result[1][:6], support.set_bit(3))
        self.assertEquals(result[2][:6], support.set_bit(4))
        self.assertEquals(result[3][:6], support.set_bit(5))
        self.assertEquals(result[4][:6], support.set_bit(9))
        ## This appears to be a bug in the OpenBabel MACCS definition
        if VERSION in ("2.2.3", "2.3.0"):
            # This is WRONG, since OB has an off-by-one error in the ring sizes
            self.assertEquals(result[5][:6], "000020")
        else:
            # which is fixed in the SVN version
            self.assertEquals(result[5][:6], support.set_bit(10))
        self.assertEquals(result[6][:6], support.set_bit(16))

TestMACCS = unittest2.skipIf(skip_openbabel, "OpenBabel not installed")(TestMACCS)

if __name__ == "__main__":
    unittest2.main()
