from __future__ import with_statement
import sys
import unittest2
import tempfile
import shutil
import os

import support

try:
    import rdkit
    has_rdkit = True
    skip_rdkit = False
except ImportError:
    has_rdkit = False
    skip_rdkit = True

    if not support.can_skip("rdkit"):
        skip_rdkit = False
        import rdkit

if has_rdkit:
    from chemfp.commandline import rdkit2fps

    runner = support.Runner(rdkit2fps.main)
else:
    runner = None

MACCS_SMI = support.fullpath("maccs.smi")

class TestMACCS(unittest2.TestCase):
    def test_bitorder(self):
        result = runner.run_fps("--maccs166", 7, MACCS_SMI)
        # The fingerprints are constructed to test the first few bytes.
        self.assertEquals(result[0][:6], support.set_bit(2))
        self.assertEquals(result[1][:6], support.set_bit(3))
        self.assertEquals(result[2][:6], support.set_bit(4))
        self.assertEquals(result[3][:6], support.set_bit(5))
        self.assertEquals(result[4][:6], support.set_bit(9))
        self.assertEquals(result[5][:6], support.set_bit(10))
        self.assertEquals(result[6][:6], support.set_bit(16))
    def test_type(self):
        for line in runner.run("--maccs166", MACCS_SMI):
            if line.startswith("#type="):
                self.assertEquals(line, "#type=RDKit-MACCS166/1")
                return
        self.assertEquals("could not find", "#type line")

TestMACCS = unittest2.skipIf(skip_rdkit, "RDKit not installed")(TestMACCS)

_fp1 = "32bb93145be9598dc6f22cbd1c781196e1733f7a53ed6f09e9e55e22bd3d3ac9e3be17f187fbcaefea8d2982ba7dab47ae1a3fd1aca52b48c70f540f964f79cd79afd9dc9871717341eaf7d7abe6febbc9bee9a971855ec7d960ecb2dacdbbb9b9b6d05f8ce9b7f4bc57fa7fa4573e95fe5a7dc918883f7fd9a3a825ef8e2fb2df944b94a2fb36c023cef883e967d9cf698fbb927cfe4fcbbaff71f7ada5ced97d5d679764bba6be8ff7d762f98d26bfbb3cb003647e1180966bc7eaffdad9a2ce47c6169bf679639e67e1bf50bd8bf30d3438dc877e67ba4e786fedfb831e56f34abc27bdfdce02c7aa57b36f761deb9d9bd5b2579df169ab0eae547515d2a7"

assert len(_fp1) == 2048 // 4

def get_field_and_first(cmdline, field):
    result = runner.run(cmdline)
    field_value = None
    first = None
    for line in result:
        if line.startswith(field):
            field_value = line
        if not line.startswith("#"):
            first = line
            break
    return (field_value, first)

class TestRDKFingerprints(unittest2.TestCase):
    def assertIn(self, substr, str):
        self.assertEquals(substr in str, True, str)
        
    def test_is_default(self):
        result = runner.run_fps("", 19)
        self.assertEquals(result[0], _fp1 + "\t9425004")
        self.assertNotEquals(result[1].split()[0], _fp1)
        # All must have the same length (since the fp lengths and ids lengths are the same
        self.assertEquals(len(set(map(len, result))), 1, set(map(len, result)))

    def test_as_rdk(self):
        result = runner.run_fps("--RDK", 19)
        self.assertEquals(result[0], _fp1 + "\t9425004")
        self.assertNotEquals(result[1].split()[0], _fp1)
        # All must have the same length (since the fp lengths and ids lengths are the same
        self.assertEquals(len(set(map(len, result))), 1, set(map(len, result)))

    def test_num_bits_default(self):
        result = runner.run_fps("--fpSize 2048", 19)
        self.assertEquals(result[0], _fp1 + "\t9425004")
        self.assertNotEquals(result[1].split()[0], _fp1)

    def test_num_bits_16(self):
        field, first = get_field_and_first("--fpSize 16", "#num_bits=")
        self.assertEquals(field, "#num_bits=16")
        self.assertEquals(first, "ffff\t9425004")

    def test_num_bits_1(self):
        field, first = get_field_and_first("--fpSize 1", "#num_bits=")
        self.assertEquals(field, "#num_bits=1")
        self.assertEquals(first, "01\t9425004")

    def test_num_bits_2(self):
        field, first = get_field_and_first("--fpSize 2", "#num_bits=")
        self.assertEquals(field, "#num_bits=2")
        self.assertEquals(first, "03\t9425004")

    def test_num_bits_too_small(self):
        result = runner.run_exit("--fpSize 0")
        self.assertIn("fpSize must be 1 or greater", result)

    def test_bits_per_hash_default(self):
        field, first = get_field_and_first("--nBitsPerHash 4", "#type=")
        self.assertEquals(field,
  "#type=RDKit-Fingerprint/1 minPath=1 maxPath=7 fpSize=2048 nBitsPerHash=4 useHs=1")
        self.assertEquals(first.split()[0], _fp1)

    def test_bits_per_hash(self):
        field, first = get_field_and_first("--nBitsPerHash 1", "#type")
        self.assertEquals(field,
  "#type=RDKit-Fingerprint/1 minPath=1 maxPath=7 fpSize=2048 nBitsPerHash=1 useHs=1")
        self.assertNotEquals(first.split()[0], _fp1)

    def test_bits_per_hash_too_small(self):
        result = runner.run_exit("--nBitsPerHash 0")
        self.assertIn("nBitsPerHash must be 1 or greater", result)

    def test_min_path_default(self):
        field, first = get_field_and_first("--minPath 1", "#type")
        self.assertEquals(field,
  "#type=RDKit-Fingerprint/1 minPath=1 maxPath=7 fpSize=2048 nBitsPerHash=4 useHs=1")
        self.assertEquals(first.split()[0], _fp1)

    def test_min_path_2(self):
        field, first = get_field_and_first("--minPath 2", "#type")
        self.assertEquals(field,
  "#type=RDKit-Fingerprint/1 minPath=2 maxPath=7 fpSize=2048 nBitsPerHash=4 useHs=1")
        self.assertNotEquals(first.split()[0], _fp1)

    def test_min_path_too_small(self):
        result = runner.run_exit("--minPath 0")
        self.assertIn("minPath must be 1 or greater", result)

    def test_min_path_too_large(self):
        result = runner.run_exit("--minPath 5 --maxPath 4")
        self.assertIn("--minPath must not be greater than --maxPath", result)

    def test_max_path_default(self):
        field, first = get_field_and_first("--maxPath 7", "#type")
        self.assertEquals(field,
  "#type=RDKit-Fingerprint/1 minPath=1 maxPath=7 fpSize=2048 nBitsPerHash=4 useHs=1")
        self.assertEquals(first.split()[0], _fp1)

    def test_max_path_6(self):
        field, first = get_field_and_first("--maxPath 6", "#type")
        self.assertEquals(field,
  "#type=RDKit-Fingerprint/1 minPath=1 maxPath=6 fpSize=2048 nBitsPerHash=4 useHs=1")
        self.assertNotEquals(first.split()[0], _fp1)

#    def test_ignore_Hs(self):
#  I don't have a good test case for this... XXX

TestRDKFingerprints = unittest2.skipIf(skip_rdkit, "RDKit not installed")(TestRDKFingerprints)


_morgan1 = "00000080000200000010010000000040000000400000000000000000000000000000000200000000800000400040000400000000000000004200000000080000000000000000020000000000000000000004000000004008000000000002000000000000800000000800000100080800000000048000000000400000000000000081000002000000010000000000000001000020000000000000000000020000000000000100000020040800100000000000000000000000000000000000000000000000000000040000800000000000000008000000000408004000000000000000000000000100000002000000002000010000100000000000000000000000"

_morgan_radius3 = "00000080000200000110010000000040000100400000000000000000000000004000000201000000800000400040000401000000000000004200000000080000000000000000020000000000000000000004400000004008000000000002000000000000800000000800000100080800000000048000000000408000000000000081000002000000010020000000000001000020000000000002000000020080000000000100000020040800100000000000000000000000000000000000200000000040000000040000800000000000000008000000040408004000000000000000000000000100000002000000002800010000100000000000000000000000"

class TestRDKMorgan(unittest2.TestCase):
    def test_as_morgan(self):
        result = runner.run_fps("--morgan", 19)
        self.assertEquals(result[0], _morgan1 + "\t9425004")
        self.assertNotEquals(result[1].split()[0], _morgan1)
        # All must have the same length (since the fp lengths and ids lengths are the same
        self.assertEquals(len(set(map(len, result))), 1, set(map(len, result)))
        
    def test_num_bits_default(self):
        result = runner.run_fps("--morgan --fpSize 2048", 19)
        self.assertEquals(result[0], _morgan1 + "\t9425004")
        self.assertNotEquals(result[1].split()[0], _morgan1)

    def test_num_bits_16(self):
        field, first = get_field_and_first("--morgan --fpSize 16", "#num_bits=")
        self.assertEquals(field, "#num_bits=16")
        self.assertEquals(first, "fbff\t9425004")

    def test_num_bits_1(self):
        field, first = get_field_and_first("--morgan --fpSize 1", "#num_bits=")
        self.assertEquals(field, "#num_bits=1")
        self.assertEquals(first, "01\t9425004")

    def test_num_bits_2(self):
        field, first = get_field_and_first("--morgan --fpSize 2", "#num_bits=")
        self.assertEquals(field, "#num_bits=2")
        self.assertEquals(first, "03\t9425004")

    def test_num_bits_too_small(self):
        result = runner.run_exit("--morgan --fpSize 0")
        self.assertIn("fpSize must be 1 or greater", result)

    def test_radius_default(self):
        result = runner.run_fps("--morgan --radius 2", 19)
        self.assertEquals(result[0], _morgan1 + "\t9425004")
        self.assertNotEquals(result[1].split()[0], _morgan1)

    def test_radius_3(self):
        result = runner.run_fps("--morgan --radius 3", 19)
        self.assertEquals(result[0], _morgan_radius3 + "\t9425004")
        self.assertNotEquals(result[1].split()[0], _morgan1)

    def test_radius_too_small(self):
        result = runner.run_exit("--morgan --radius -1")
        self.assertIn("radius must be 0 or greater", result)

    def test_default_use_options(self):
        field, first = get_field_and_first("--morgan --useFeatures 0 --useChirality 0 --useBondTypes 1",
                                      "#type")
        self.assertEquals(field,
                          "#type=RDKit-Morgan/1 radius=2 fpSize=2048 useFeatures=0 useChirality=0 useBondTypes=1")
        self.assertEquals(first, _morgan1 + "\t9425004")

    # This isn't a complete test of the different options. I don't think it's worth the effort
    def test_useChirality(self):
        field, first = get_field_and_first("--morgan --useFeatures 1 --useChirality 1 --useBondTypes 0",
                                           "#type=")
        self.assertEquals(field,
                          "#type=RDKit-Morgan/1 radius=2 fpSize=2048 useFeatures=1 useChirality=1 useBondTypes=0")
        self.assertNotEquals(first, _morgan1 + "\t9425004")
        
TestRDKMorgan = unittest2.skipIf(skip_rdkit, "RDKit not installed")(TestRDKMorgan)

class TestIO(unittest2.TestCase, support.TestIdAndErrors):
    _runner = runner
    def test_input_format(self):
        def without_source_header(cmdline, source):
            return [line for line in runner.run(cmdline, source)
                        if not line.startswith("#source=") and
                           not line.startswith("#date=")]
        result1 = without_source_header("", support.PUBCHEM_SDF)
        result2 = without_source_header("", support.PUBCHEM_SDF_GZ)
        self.assertEquals(result1, result2)

        result3 = without_source_header("--in sdf.gz", support.PUBCHEM_SDF_GZ)
        self.assertEquals(result1, result3)
        
        result4 = without_source_header("--in sdf", support.PUBCHEM_ANOTHER_EXT)
        self.assertEquals(result1, result4)

    def test_output(self):
        dirname = tempfile.mkdtemp(prefix="test_rdkit2fps")
        output_filename = os.path.join(dirname, "blah.fps")
        assert len(output_filename.split()) == 1
        try:
            result = runner.run("-o " + output_filename)
            assert len(result) == 0
            with open(output_filename) as f:
                result = f.readlines()
        finally:
            shutil.rmtree(dirname)
        self.assertEquals(result[0], "#FPS1\n")
        while result and result[0].startswith("#"):
            del result[0]
        self.assertEquals(len(result), 19)
        self.assertEquals(result[0], _fp1 + "\t9425004\n")

    def test_bad_format(self):
        result = runner.run_exit("--in spam")
        self.assertIn("Unsupported format specifier: 'spam'", result)

TestIO = unittest2.skipIf(skip_rdkit, "RDKit not installed")(TestIO)

if __name__ == "__main__":
    unittest2.main()
