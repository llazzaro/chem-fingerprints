from __future__ import with_statement
import sys
import unittest2
import tempfile
import shutil
import os
from cStringIO import StringIO
import tempfile

import support

try:
    import chemfp.rdkit
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
TRP = open(support.fullpath("tryptophan.sdf")).read()

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

_atom_pair_fingerprints = {
None: None,
"1": {"2048": "100000100008007000045020000008b0080220a4420084c000054010800300e02040000088c080010000800101404023000000100000000020000000004000a00002060400040800000000000002c00108000000000030009d00100002001001900080003081058010000400200209004000000000050b00800008084042060801000800000000010200030000000000040000000080000000000400000000021000708000100000600010008080200008000c8020000004040008000000600000000008000000100000000402000000400080300000600090004020000020000002008100000800000100020000000000000000000008100000000002000000",
       "128": "fd42febdfaddfff5ff05df3fe3c3fffb",
       "minLength": "dd02bebd328cbff5be055f3e6242ff32",
       "maxLength": "3042c000e8d1e141f101d02181c160cb",
       },

"2": {"2048": "0100070010000000101100000013010000000000001100000010000001001703100000000007000011000301000310001000000010000000003000110100731000310001300000000000101000033110100000010000001000037311000000370313033003010000000101070000130030000010330000000000170031001077000013301000000300003133030000300030133003000131011100100f000010000013010300000000030310310000000300101030000011033010077100100000300003000000000011000000000110000010000301000037300000000101000001303000000000003000000010000010000001100000073001100100101010",
       "128": "77f7fff7ff17017f7fffff7fff3fffff",
       "minLength": "71777f777317003377ffff37733fff77",
       "maxLength": "073033737f0001370100337f7f101077",
       },
}
if not (skip_rdkit or chemfp.rdkit.ATOM_PAIR_VERSION is None):
    _atom_pair_fps = _atom_pair_fingerprints[chemfp.rdkit.ATOM_PAIR_VERSION]
    PAIR_TYPE = "RDKit-AtomPair/" + chemfp.rdkit.ATOM_PAIR_VERSION + " "

class TestAtomPairFingerprinter(unittest2.TestCase):
    def test_pair_defaults(self):
        header, output = runner.run_split("--pair", 19)
        self.assertEqual(header["#type"], PAIR_TYPE + "fpSize=2048 minLength=1 maxLength=30")
        self.assertEqual(output[0], _atom_pair_fps["2048"] + "\t9425004")
    def test_pair_explicit_defaults(self):
        header, output = runner.run_split("--pair --fpSize 2048 --minLength 1 --maxLength 30", 19)
        self.assertEqual(header["#type"], PAIR_TYPE + "fpSize=2048 minLength=1 maxLength=30")
        self.assertEqual(output[0], _atom_pair_fps["2048"] + "\t9425004")
    def test_num_bits_128(self):
        header, output = runner.run_split("--pair --fpSize 128", 19)
        self.assertEqual(header["#type"], PAIR_TYPE + "fpSize=128 minLength=1 maxLength=30")
        self.assertEqual(output[0], _atom_pair_fps["128"] + "\t9425004")
    def test_num_bits_error(self):
        errmsg = runner.run_exit("--pair --fpSize 0")
        self.assertIn("fpSize must be 1 or greater", errmsg)
        errmsg = runner.run_exit("--pair --fpSize 2.3")
        self.assertIn("fpSize must be 1 or greater", errmsg)
    def test_min_length(self):
        header, output = runner.run_split("--pair --fpSize 128 --minLength 4", 19)
        self.assertEqual(header["#type"], PAIR_TYPE + "fpSize=128 minLength=4 maxLength=30")
        self.assertEqual(output[0], _atom_pair_fps["minLength"] + "\t9425004")
    def test_max_length(self):
        header, output = runner.run_split("--pair --fpSize 128 --maxLength 3", 19)
        self.assertEqual(header["#type"], PAIR_TYPE + "fpSize=128 minLength=1 maxLength=3")
        self.assertEqual(output[0], _atom_pair_fps["maxLength"] + "\t9425004")
    def test_min_length_error(self):
        errmsg = runner.run_exit("--pair --minLength spam")
        self.assertIn("minLength must be 0 or greater", errmsg)
    def test_max_length_error(self):
        errmsg = runner.run_exit("--pair --maxLength -3")
        self.assertIn("maxLength must be 0 or greater", errmsg)
        errmsg = runner.run_exit("--pair --maxLength spam")
        self.assertIn("maxLength must be 0 or greater", errmsg)
    def test_invalid_min_max_lengths(self):
        errmsg = runner.run_exit("--pair --maxLength 0") # default minLength is 1
        self.assertIn("--minLength must not be greater than --maxLength", errmsg)
        errmsg = runner.run_exit("--pair --minLength 4 --maxLength 3")
        self.assertIn("--minLength must not be greater than --maxLength", errmsg)
    def test_valid_min_max_lengths(self):
        header, output = runner.run_split("--pair --minLength 0")
        self.assertEqual(header["#type"], PAIR_TYPE + "fpSize=2048 minLength=0 maxLength=30")
        header, output = runner.run_split("--pair --minLength 0 --maxLength 0")
        self.assertEqual(header["#type"], PAIR_TYPE + "fpSize=2048 minLength=0 maxLength=0")
        header, output = runner.run_split("--pair --minLength 5 --maxLength 5")
        self.assertEqual(header["#type"], PAIR_TYPE + "fpSize=2048 minLength=5 maxLength=5")
        header, output = runner.run_split("--pair --minLength 6 --maxLength 8")
        self.assertEqual(header["#type"], PAIR_TYPE + "fpSize=2048 minLength=6 maxLength=8")


if skip_rdkit:
    TestAtomPairFingerprinter = unittest2.skipIf(skip_rdkit, "RDKit not installed")(TestAtomPairFingerprinter)
else:
    TestAtomPairFingerprinter = unittest2.skipIf(chemfp.rdkit.ATOM_PAIR_VERSION is None,
   "This version of RDKit has a broken GetHashedAtomPairFingerprintAsBitVect")(TestAtomPairFingerprinter)


_torsion_fingerprints = {
"1": {
    "2048": "000000010100000000000040800008000000000000000000000000000000000040004000000000000c0000000000000040000000003010000000000000000800000000000000000000000000000000000000000000100000000000000000000000000080000000000000000000000000000000000000000000000000400000000000000000000000000000000000000000000002000000000000000000000000000000000000002000000000000000000000000000000020200080000000000000100000000000000000000000000000800000000000008400000000000000000200000000000000000000020000000000000000010008000000000200000000",
    "128": "c2104083013018a42c008042c0000800",
    "targetSize": "1491150001c0f010648000086245052c",
    },
"2": {
    "2048": "00000000000003000000000001000000010000000000000000000000000000000000000000000000000000000000000700000010000000000010000000000000000000100000000000000000100000000000000000000100000000000001003000000000000000000000000000000000000000000000000000000000000000001000000100000000000000010000001000000000000000000000000000001000001100000000000000000000100000000000000000000000000000000000000000000001000000000000000000000000010000000000330000100100000010000000100000000000000000101000000000000001000000000030000000000000",
    "128": "13111033000037000070011131013037",
    "targetSize": "33037307030103730303131100331100",
      }
}
if not skip_rdkit:
    _torsion_fps = _torsion_fingerprints[chemfp.rdkit.TORSION_VERSION]
    TORSION_TYPE = "RDKit-Torsion/" + chemfp.rdkit.TORSION_VERSION + " "

class TestTorsionFingerprinter(unittest2.TestCase):
    def test_torsion_defaults(self):
        header, output = runner.run_split("--torsion", 19)
        self.assertEqual(header["#type"], "RDKit-Torsion/1 fpSize=2048 targetSize=4")
        self.assertEqual(output[0], _torsion_fps["2048"] + "\t9425004")
    def test_torsion_explicit_defaults(self):
        header, output = runner.run_split("--torsion --fpSize 2048 --targetSize 4", 19)
        self.assertEqual(header["#type"], "RDKit-Torsion/1 fpSize=2048 targetSize=4")
        self.assertEqual(output[0], _torsion_fps["2048"] + "\t9425004")
    def test_num_bits_128(self):        
        header, output = runner.run_split("--torsion --fpSize 128 --targetSize 4", 19)
        self.assertEqual(header["#type"], "RDKit-Torsion/1 fpSize=128 targetSize=4")
        self.assertEqual(output[0], _torsion_fps["128"] + "\t9425004")
    def test_num_bits_error(self):
        errmsg = runner.run_exit("--torsion --fpSize 0")
        self.assertIn("fpSize must be 1 or greater", errmsg)
        errmsg = runner.run_exit("--torsion --fpSize 2.3")
        self.assertIn("fpSize must be 1 or greater", errmsg)
    def test_target_size(self):
        header, output = runner.run_split("--torsion --fpSize 128 --targetSize 5", 19)
        self.assertEqual(header["#type"], "RDKit-Torsion/1 fpSize=128 targetSize=5")
        self.assertEqual(output[0], _torsion_fps["targetSize"] + "\t9425004")
    def test_target_size_error(self):
        errmsg = runner.run_exit("--torsion --fpSize 128 --targetSize -1")
        self.assertIn("targetSize must be 1 or greater", errmsg)
        errmsg = runner.run_exit("--torsion --fpSize 128 --targetSize spam")
        self.assertIn("targetSize must be 1 or greater", errmsg)
        errmsg = runner.run_exit("--torsion --fpSize 128 --targetSize 0")
        self.assertIn("targetSize must be 1 or greater", errmsg)

TestTorsionFingerprinter = unittest2.skipIf(skip_rdkit, "RDKit not installed")(TestTorsionFingerprinter)

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


class TestBadStructureFiles(unittest2.TestCase):
    def setUp(self):
        self.dirname = tempfile.mkdtemp()
    def tearDown(self):
        shutil.rmtree(self.dirname)

    def _make_datafile(self, text, ext):
        filename = os.path.join(self.dirname, "input."+ext)
        with open(filename, "w") as f:
            f.write(text)
        return filename
    
    def test_blank_line_in_smiles(self):
        filename = self._make_datafile("C methane\n\nO water\n", "smi")
        msg = runner.run_exit([filename])
        self.assertIn("Unexpected blank line at line 2 of", msg)

    def test_bad_smiles(self):
        filename = self._make_datafile("C methane\nQ Q-ane\nO water\n", "smi")
        msg = runner.run_exit([filename])
        self.assertIn("Cannot parse the SMILES 'Q' at line 2", msg)

    def test_smiles_without_title(self):
        filename = self._make_datafile("C methane\nO water\n[235U]\n", "smi")
        msg = runner.run_exit([filename])
        self.assertIn("Missing SMILES name (second column) at line 3", msg)

    def test_sdf_with_bad_record(self):
        # Three records, second one is bad
        input = TRP + TRP.replace("32 28", "40 21") + TRP
        filename = self._make_datafile(input, "sdf")
        msg = runner.run_exit([filename])
        self.assertIn("Could not parse molecule block at line 70", msg)
        self.assertIn("input.sdf", msg)

    def test_sdf_with_bad_record_checking_id(self):
        # This tests a different code path than the previous
        input = TRP + TRP.replace("32 28", "40 21") + TRP
        filename = self._make_datafile(input, "sdf")
        msg = runner.run_exit(["--id-tag", "COMPND", filename])
        self.assertIn("Could not parse molecule block at line 70", msg)
        self.assertIn("input.sdf", msg)

    def test_sdf_with_missing_id(self):
        filename = self._make_datafile(TRP, "sdf")
        msg = runner.run_exit(["--id-tag", "SPAM", filename])
        self.assertIn("Missing id tag 'SPAM' for record #1 at line 1", msg)
        self.assertIn("input.sdf", msg)

    def test_ignore_errors(self):
        input = TRP + TRP.replace("32 28", "40 21") + TRP.replace("COMPND", "BLAH")
        filename = self._make_datafile(input, "sdf")
        header, output = runner.run_split(["--errors", "ignore", "--id-tag", "BLAH"], source=filename)
        self.assertEqual(len(output), 1)

    def test_unsupported_format(self):
        filename = self._make_datafile("Unknown", "xyzzy")
        result = runner.run_exit([filename])
        self.assertIn("Unknown structure filename extension", result)
        self.assertIn("input.xyzzy", result)
        
TestBadStructureFiles = unittest2.skipIf(skip_rdkit, "RDKit not installed")(TestBadStructureFiles)
        

# Some code to test the internal interface
class TestInternals(unittest2.TestCase):
    def test_make_rdk_fingerprinter(self):
        # Make sure that I can call with the defaults
        chemfp.rdkit.make_rdk_fingerprinter()
            
    def test_make_rdk_fingerprinter_bad_fpSize(self):
        with self.assertRaisesRegexp(ValueError, "fpSize must be positive"):
            chemfp.rdkit.make_rdk_fingerprinter(fpSize=0)
        with self.assertRaisesRegexp(ValueError, "fpSize must be positive"):
            chemfp.rdkit.make_rdk_fingerprinter(fpSize=-10)

    def test_make_rdk_fingerprinter_min_path(self):
        with self.assertRaisesRegexp(ValueError, "minPath must be positive"):
            chemfp.rdkit.make_rdk_fingerprinter(minPath=0)
        with self.assertRaisesRegexp(ValueError, "minPath must be positive"):
            chemfp.rdkit.make_rdk_fingerprinter(monPath=-3)

    def test_make_rdk_fingerprinter_max_path(self):
        chemfp.rdkit.make_rdk_fingerprinter(minPath=2, maxPath=2)
        with self.assertRaisesRegexp(ValueError, "maxPath must not be smaller than minPath"):
            chemfp.rdkit.make_rdk_fingerprinter(minPath=3, maxPath=2)

    def test_make_rdk_fingerprinter_min_path(self):
        with self.assertRaisesRegexp(ValueError, "nBitsPerHash must be positive"):
            chemfp.rdkit.make_rdk_fingerprinter(nBitsPerHash=0)
        with self.assertRaisesRegexp(ValueError, "nBitsPerHash must be positive"):
            chemfp.rdkit.make_rdk_fingerprinter(nBitsPerHash=-1)


    def test_make_morgan_fingerprinter(self):
        chemfp.rdkit.make_morgan_fingerprinter()
        
    def test_make_morgan_fingerprinter_bad_fpSize(self):
        with self.assertRaisesRegexp(ValueError, "fpSize must be positive"):
            chemfp.rdkit.make_morgan_fingerprinter(fpSize=0)
        with self.assertRaisesRegexp(ValueError, "fpSize must be positive"):
            chemfp.rdkit.make_morgan_fingerprinter(fpSize=-10)

    def test_make_morgan_fingerprinter_bad_radius(self):
        with self.assertRaisesRegexp(ValueError, "radius must be positive or zero"):
            chemfp.rdkit.make_morgan_fingerprinter(radius=-1)
        with self.assertRaisesRegexp(ValueError, "radius must be positive or zero"):
            chemfp.rdkit.make_morgan_fingerprinter(radius=-10)

    
TestInternals = unittest2.skipIf(skip_rdkit, "RDKit not installed")(TestInternals)
        
if __name__ == "__main__":
    unittest2.main()
