import unittest2
import sys
import os
from cStringIO import StringIO as SIO

import support

try:
    from openeye import oechem  # These tests require OEChem
    has_oechem = True
    skip_oechem = False
except ImportError:
    has_oechem = False
    skip_oechem = True
    if not support.can_skip("oe"):
        skip_oechem = False
        from openeye import oechem

if has_oechem:
    from chemfp.commandline import oe2fps
    import chemfp.openeye
    chemfp.openeye._USE_SELECT = False # Grrr. Needed to automate testing.

    real_stdout = sys.stdout
    real_stderr = sys.stderr

    PUBCHEM_SDF = support.fullpath("pubchem.sdf")
    PUBCHEM_SDF_GZ = support.fullpath("pubchem.sdf.gz")
    PUBCHEM_ANOTHER_EXT = support.fullpath("pubchem.should_be_sdf_but_is_not")


    oeerrs = oechem.oeosstream()
    oechem.OEThrow.SetOutputStream(oeerrs)
    
def _check_for_oe_errors():
    lines = oeerrs.str().splitlines()
    for line in lines:
        if line.startswith("Warning: Stereochemistry corrected on atom number"):
            continue
        if line.startswith("Warning: Unknown file format set in input stream"):
            # There's a bug in OEChem where it generates this warning on unknown
            # file extensions even after SetFormat has been called
            continue
        raise AssertionError("Unexpected message from OEChem: %r" % (line,))


# I build the fingerprints using bit offsets to ensure that the test
# data matches the actual bit results from OEChem. While I could
# reproduce the method in chemfp.openeye.get_maccs_fingerprinter, that
# would be cheating. I could test against ToHexString() but then I
# would have the nagging feeling that I got the ordering
# backwards. Instead, I do it from scratch using the bit offset.

def _construct_test_values():
    from openeye.oechem import oemolistream
    from openeye.oegraphsim import OEFingerPrint, OEMakePathFP
    fp = OEFingerPrint()
    ifs = oemolistream()
    assert ifs.open(PUBCHEM_SDF)
    hex_data = []

    def _convert_to_chemfp_order(s):
        # The FPS format allows either case but I prefer lowercase
        s = s.lower()
        # OpenEye orders hex values on nibbles. Chemfp orders on bytes.
        return "".join( (s[i+1]+s[i]) for i in range(0, len(s), 2))

    for mol in ifs.GetOEGraphMols():
        OEMakePathFP(fp, mol)
        # Set the byte values given the bit offsets
        bytes = [0] * (4096//8)
        i = fp.FirstBit()
        while i >= 0:
            bytes[i//8] |= 1<<(i%8)
            i = fp.NextBit(i)
        as_hex = "".join("%02x" % i for i in bytes)
        assert len(as_hex) == 2*(4096//8), len(as_hex)
        # Double-check that it matches the (reordered) ToHexString()
        oe_hex = fp.ToHexString()[:-1]
        assert as_hex == _convert_to_chemfp_order(oe_hex), (
            as_hex, _convert_to_chemfp_order(oe_hex))
        
        hex_data.append("%s %s" % (as_hex, mol.GetTitle()))
    return hex_data

# I have this to flag any obvious changes in the OEChem algorithm and
# to help with figuring out how to build a test case.

_fp1 = "00001002200200000000000000000000000008400020000300801002300000000200000840000000000080000000000000204008000000000010000c10000000400000010100000210800002000000009400000000020020088000000000010000918000200000580400002000010020002440000008001001404000000200010000a8c00020400200002000004084000000030100820000000000000002000000510001800000010001000081100110000800480000100400000c00004c000800000808000100000022000228800020004000000200182100000100000000101000010004004808000000800000000001010010201000000090400000100000020010000010201000000040300100000000580000000000000200000000401000000000000008040004000000008002080820000310280200004040a000000010000080005000004010010018000000800000020008208040000400000200000000000000000800080050000008400100004000000200ac001000000000800100200060900010002000000040200000000000040808000048400040000000020000001001000000000302002008200000a044000180800000100000000200000049004080080000100022a00084000400280480000000402400080400404100000000040000020c10000000000c000100002000080010002080100002000600"

if has_oechem:
    hex_test_values = _construct_test_values()
    assert hex_test_values[0].startswith(_fp1)

class OERunner(support.Runner):
    def pre_run(self):
        oeerrs.clear()
    def post_run(self):
        _check_for_oe_errors()

    def run_stdin(self, cmdline, source):
        fd = os.open(source, os.O_RDONLY)
        oechem.oein.openfd(fd, 0)
        try:
            return self.run(cmdline, None)
        finally:
            oechem.oein.openfd(0, 0)
            os.close(fd)

if has_oechem:
    runner = OERunner(oe2fps.main)
    run = runner.run
    run_stdin = runner.run_stdin
    run_fps = runner.run_fps
    run_exit = runner.run_exit

def headers(lines):
    assert lines[0] == "FPS1"
    del lines[0]
    return [line for line in lines if line.startswith("#")]


class TestMACCS(unittest2.TestCase):
    def test_bitorder(self):
        result = run_fps("--maccs166", 7, support.fullpath("maccs.smi"))
        # The fingerprints are constructed to test the first few bytes.
        self.assertEquals(result[0][:6], support.set_bit(2))
        self.assertEquals(result[1][:6], support.set_bit(3))
        self.assertEquals(result[2][:6], support.set_bit(4))
        self.assertEquals(result[3][:6], support.set_bit(5))
        self.assertEquals(result[4][:6], support.set_bit(9))
        self.assertEquals(result[5][:6], support.set_bit(10))
        self.assertEquals(result[6][:6], support.set_bit(16))

TestMACCS = unittest2.skipIf(skip_oechem, "OEChem not installed")(TestMACCS)

class TestPath(unittest2.TestCase):
    def test_default(self):
        result = run_fps("", 19)
        hexfp, id = result[0].split()
        self.assertEquals(len(hexfp), 4096//4)
        self.assertEquals(result[0], hex_test_values[0])
        self.assertEquals(result, hex_test_values)

    def test_path_option(self):
        result = run_fps("--path", 19)
        self.assertEquals(result, hex_test_values)

    def test_num_bits(self):
        result = run_fps("--numbits 16", 19)
        self.assertEquals(result[0][:5], "ff1f ")

    def test_min_bonds_default(self):
        result = run_fps("--minbonds 0", 19)
        self.assertEquals(result, hex_test_values)

    def test_min_bonds_changed(self):
        result = run_fps("--minbonds 1", 19)
        self.assertNotEquals(result, hex_test_values)

    def test_max_bonds_default(self):
        result = run_fps("--maxbonds 5", 19)
        self.assertEquals(result, hex_test_values)
        
    def test_max_bonds_changed(self):
        result = run_fps("--minbonds 4", 19)
        self.assertNotEquals(result, hex_test_values)

    def test_atype_default_named(self):
        result = run_fps("--atype DefaultAtom", 19)
        self.assertEquals(result, hex_test_values)

    def test_atype_default_flags(self):
        result = run_fps("--atype Aromaticity|AtomicNumber|Chiral|EqHalogen|"
                         "FormalCharge|HvyDegree|Hybridization", 19)
        self.assertEquals(result, hex_test_values)


    def test_atype_default_flags_with_duplicates(self):
        result = run_fps("--atype Aromaticity|Chiral|AtomicNumber|AtomicNumber|EqHalogen|"
                         "HvyDegree|FormalCharge|Hybridization", 19)
        self.assertEquals(result, hex_test_values)

    # Make sure that each of the flags returns some other answer
    def test_atype_different_1(self):
        result = run_fps("--atype AtomicNumber|Chiral|EqHalogen|"
                         "FormalCharge|HvyDegree|Hybridization", 19)
        self.assertNotEquals(result, hex_test_values)

    def test_atype_different_2(self):
        result = run_fps("--atype Aromaticity|Chiral|EqHalogen|"
                         "FormalCharge|HvyDegree|Hybridization", 19)
        self.assertNotEquals(result, hex_test_values)

    def test_atype_different_3(self):
        result = run_fps("--atype Aromaticity|AtomicNumber|EqHalogen|"
                         "FormalCharge|HvyDegree|Hybridization", 19)
        self.assertNotEquals(result, hex_test_values)

    def test_atype_different_4(self):
        result = run_fps("--atype Aromaticity|AtomicNumber|Chiral|"
                         "FormalCharge|HvyDegree|Hybridization", 19)
        self.assertNotEquals(result, hex_test_values)

    def test_atype_different_5(self):
        result = run_fps("--atype Aromaticity|AtomicNumber|Chiral|EqHalogen|"
                         "HvyDegree|Hybridization", 19)
        self.assertNotEquals(result, hex_test_values)

    def test_atype_different_6(self):
        result = run_fps("--atype Aromaticity|AtomicNumber|Chiral|EqHalogen|"
                         "FormalCharge|Hybridization", 19)
        self.assertNotEquals(result, hex_test_values)

    def test_atype_different_7(self):
        result = run_fps("--atype Aromaticity|AtomicNumber|Chiral|EqHalogen|"
                         "FormalCharge|HvyDegree", 19)
        self.assertNotEquals(result, hex_test_values)


    def test_btype_default_named(self):
        result = run_fps("--btype DefaultBond", 19)
        self.assertEquals(result, hex_test_values)

    def test_btype_default_flags(self):
        result = run_fps("--btype BondOrder|Chiral", 19)
        self.assertEquals(result, hex_test_values)

    def test_btype_different_1(self):
        result = run_fps("--btype BondOrder", 19)
        self.assertNotEquals(result, hex_test_values)

    def test_btype_different_2(self):
        result = run_fps("--btype Chiral", 19)
        self.assertNotEquals(result, hex_test_values)

TestPath = unittest2.skipIf(skip_oechem, "OEChem not installed")(TestPath)

class TestIO(unittest2.TestCase):
    def test_compressed_input(self):
        result = run_fps("", source=PUBCHEM_SDF_GZ)
    def test_unknown_extension(self):
        # OEChem's default assumes SMILES. This will parse some of the
        # SD file lines as SMILES and skip the ones it doesn't know.
        # The error output will have a bunch of warnings, starting
        # with the "Unknown file format ... " warning, and then this
        # string about a SMILES parse error.
        try:
            run("", source=PUBCHEM_ANOTHER_EXT)
        except AssertionError, x:
            self.assertEquals("Problem parsing SMILES" in str(x), True, str(x))
            
    def test_specify_input_format(self):
        result = run_fps("--in sdf", source=PUBCHEM_ANOTHER_EXT)


    def test_from_stdin(self):
        run_stdin("--in sdf", source=PUBCHEM_SDF)

    def test_from_gziped_stdin(self):
        run_stdin("--in sdf.gz", source=PUBCHEM_SDF_GZ)

    def test_unknown_format(self):
        msg = run_exit("--in blah")
        self.assertEquals("Unknown format 'blah'" in msg, True, msg)

    def test_file_does_not_exist(self):
        msg = run_exit("", source="/asdfaserwe.does.not.exist")
        self.assertEquals("No such file or directory" in msg, True, repr(msg))
        self.assertEquals("asdfaserwe.does.not.exist" in msg, True, repr(msg))

# XXX how to test that this generates a warning?
#    def test_specify_input_format_with_dot(self):
#        result = run_fps("--in .sdf", source=PUBCHEM_ANOTHER_EXT)

TestIO = unittest2.skipIf(skip_oechem, "OEChem not installed")(TestIO)

class TestArgErrors(unittest2.TestCase):
    def _run(self, cmd, expect):
        msg = run_exit(cmd)
        self.assertEquals(expect in msg, True, msg)

    def test_two_fp_types(self):
        self._run("--maccs166 --path", "Cannot specify both --maccs166 and --path")

    def test_num_bits_too_small(self):
        self._run("--numbits 0", "between 16 and 65536 bits")
        self._run("--numbits 1", "between 16 and 65536 bits")
        self._run("--numbits 15", "between 16 and 65536 bits")

    def test_num_bits_too_large(self):
        self._run("--numbits 65537", "between 16 and 65536 bits")
        # Check for overflow, even though I know it won't happen in Python
        self._run("--numbits %(big)s"%dict(big=2**32+32), "between 16 and 65536 bits")

    def test_min_bonds_too_small(self):
        self._run("--minbonds=-1", "0 or greater")

    def test_min_bonds_larger_than_default_max_bonds(self):
        self._run("--minbonds=6", "--maxbonds must not be smaller than --minbonds")

    def test_min_bonds_too_large(self):
        self._run("--minbonds=4 --maxbonds=3",
                  "--maxbonds must not be smaller than --minbonds")

    def test_bad_atype(self):
        self._run("--atype spam", "Unknown atom type 'spam'")

    def test_bad_atype2(self):
        self._run("--atype DefaultAtom|spam", "Unknown atom type 'spam'")

    def test_bad_atype3(self):
        self._run("--atype DefaultAtom|", "Missing atom flag")

    def test_bad_btype(self):
        self._run("--btype eggs", "Unknown bond type 'eggs'")

    def test_bad_btype2(self):
        self._run("--btype DefaultBond|eggs", "Unknown bond type 'eggs'")

    def test_bad_btype3(self):
        self._run("--btype DefaultBond|", "Missing bond flag")

TestArgErrors = unittest2.skipIf(skip_oechem, "OEChem not installed")(TestArgErrors)

class TestHeaderOutput(unittest2.TestCase):
    def _field(self, s, field):
        result = run(s)
        filtered = [line for line in result if line.startswith(field)]
        self.assertEquals(len(filtered), 1, result)
        return filtered[0]

    def test_software(self):
        result = self._field("", "#software")
        self.assertEquals("#software=OEGraphSim/" in result, True, result)
        self.assertEquals("(" in result, True, result)
        self.assertEquals(")" in result, True, result)
        result = self._field("--maccs166", "#software")
        self.assertEquals("#software=OEGraphSim/" in result, True, result)

    def test_type(self):
        result = self._field("", "#type")
        self.assertEquals(result,
  "#type=OpenEye-Path/1 numbits=4096 minbonds=0 maxbonds=5 atype=DefaultAtom btype=DefaultBond")

    def test_default_atom_and_bond(self):
        result = self._field(
            "--atype=Aromaticity|AtomicNumber|Chiral|EqHalogen|FormalCharge|HvyDegree|Hybridization "
            "--btype=BondOrder|Chiral", "#type")
        self.assertEquals(result,
  "#type=OpenEye-Path/1 numbits=4096 minbonds=0 maxbonds=5 atype=DefaultAtom btype=DefaultBond")

        
    # different flags. All flags? and order
    def test_num_bits(self):
        result = self._field("--numbits 38", "#num_bits")
        self.assertEquals(result, "#num_bits=38")
        
    def test_atype_flags(self):
        result = self._field("--atype FormalCharge|FormalCharge", "#type") + " "
        self.assertEquals(" atype=FormalCharge " in result, True, result)
    
    def test_btype_flags(self):
        result = self._field("--btype Chiral|BondOrder", "#type") + " "
        self.assertEquals(" btype=DefaultBond " in result, True, result)
        result = self._field("--btype BondOrder|Chiral", "#type") + " "
        self.assertEquals(" btype=DefaultBond " in result, True, result)
    
    def test_pipe_or_comma(self):
        result = self._field("--atype HvyDegree,FormalCharge --btype Chiral,BondOrder",
                             "#type") + " "
        self.assertEquals(" atype=FormalCharge|HvyDegree " in result, True, result)
        self.assertEquals(" btype=DefaultBond " in result, True, result)
        
    
    def test_maccs_header(self):
        result = self._field("--maccs166", "#type")
        self.assertEquals(result, "#type=OpenEye-MACCS166/1")
    
TestHeaderOutput = unittest2.skipIf(skip_oechem, "OEChem not installed")(TestHeaderOutput)
        
if __name__ == "__main__":
    unittest2.main()
