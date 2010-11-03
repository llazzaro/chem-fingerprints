import unittest
import sys
import os
from cStringIO import StringIO as SIO

from openeye import oechem

from chemfp.commandline import oe2fps
import chemfp.openeye
chemfp.openeye._USE_SELECT = False # Grrr. Needed to automate testing.

real_stdout = sys.stdout
real_stderr = sys.stderr

PUBCHEM_SDF = "pubchem.sdf"
PUBCHEM_SDF_GZ = "pubchem.sdf.gz"
PUBCHEM_ANOTHER_EXT = "pubchem.should_be_sdf_but_is_not"


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
        raise AssertionError("Unexpected message from OEChem: {msg!r}".format(msg=line))


# I build the fingerprints using bit offsetsto ensure that the test
# data matches the actual bit results from OEChem. While I could
# reproduce the method in chemfp.openeye.get_maccs_fingerprinter, that
# would be cheating. I could test against ToHexString() but then I
# would have the nagging feeling that I got the ordering
# backwards. Instead, I do it from scratch using the bit offset.

def _construct_test_values():
    from openeye.oechem import oemolistream
    from openeye.oegraphsim import OEFingerPrint, OEMakePathFP
    fp = OEFingerPrint()
    ifs = oemolistream("pubchem.sdf")
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

hex_test_values = _construct_test_values()

# I have this to flag any obvious changes in the OEChem algorithm and
# to help with figuring out how to build a test case.

_fp1 = "00001002200200000000000000000000000008400020000300801002300000000200000840000000000080000000000000204008000000000010000c10000000400000010100000210800002000000009400000000020020088000000000010000918000200000580400002000010020002440000008001001404000000200010000a8c00020400200002000004084000000030100820000000000000002000000510001800000010001000081100110000800480000100400000c00004c000800000808000100000022000228800020004000000200182100000100000000101000010004004808000000800000000001010010201000000090400000100000020010000010201000000040300100000000580000000000000200000000401000000000000008040004000000008002080820000310280200004040a000000010000080005000004010010018000000800000020008208040000400000200000000000000000800080050000008400100004000000200ac001000000000800100200060900010002000000040200000000000040808000048400040000000020000001001000000000302002008200000a044000180800000100000000200000049004080080000100022a00084000400280480000000402400080400404100000000040000020c10000000000c000100002000080010002080100002000600"

assert hex_test_values[0].startswith(_fp1)
        

def run(s, source=PUBCHEM_SDF):
    args = s.split()
    if source is not None:
        args = args + [source]
    oeerrs.clear()
    try:
        errstr = oeerrs
        sys.stdout = stdout = SIO()
        oe2fps.main(args)
    finally:
        sys.stdout = real_stdout
    _check_for_oe_errors()
    result = stdout.getvalue().splitlines()
    if result:
        assert result[0] == "#FPS1"
    return result

def run_stdin(s, source):
    fd = os.open(source, os.O_RDONLY)
    oechem.oein.openfd(fd, 0)
    try:
        return run(s, None)
    finally:
        oechem.oein.openfd(0, 0)
        os.close(fd)

def run_fps(s, expect_length=None, source=PUBCHEM_SDF):
    result = run(s, source)
    while result[0].startswith("#"):
        del result[0]
    if expect_length is not None:
        assert len(result) == expect_length, (len(result), expect_length)
    return result


def run_exit(s, source=PUBCHEM_SDF):
    sys.stderr = stderr = SIO()
    try:
        try:
            run(s, source)
        except SystemExit:
            pass
        else:
            raise AssertionError("should have exited: %r" % (s,))
    finally:
        sys.stderr = real_stderr
    return stderr.getvalue()

def headers(lines):
    assert lines[0] == "FPS1"
    del lines[0]
    return [line for line in lines if line.startswith("#")]

def _set_bit(n):
    assert n <= 16
    bytes = [0, 0, 0]
    bytes[n//8] = 1<<(n%8)
    return "%02x%02x%02x" % tuple(bytes)

class TestMACCS(unittest.TestCase):
    def test_bitorder(self):
        result = run_fps("--maccs166", 7, "maccs.smi")
        # The fingerprints are constructed to test the first few bytes.
        self.assertEquals(result[0][:6], _set_bit(2))
        self.assertEquals(result[1][:6], _set_bit(3))
        self.assertEquals(result[2][:6], _set_bit(4))
        self.assertEquals(result[3][:6], _set_bit(5))
        self.assertEquals(result[4][:6], _set_bit(9))
        self.assertEquals(result[5][:6], _set_bit(10))
        self.assertEquals(result[6][:6], _set_bit(16))

class TestPath(unittest.TestCase):
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
        result = run_fps("--num-bits 16", 19)
        self.assertEquals(result[0][:5], "ff1f ")

    def test_min_bonds_default(self):
        result = run_fps("--min-bonds 0", 19)
        self.assertEquals(result, hex_test_values)

    def test_min_bonds_changed(self):
        result = run_fps("--min-bonds 1", 19)
        self.assertNotEquals(result, hex_test_values)

    def test_max_bonds_default(self):
        result = run_fps("--max-bonds 5", 19)
        self.assertEquals(result, hex_test_values)
        
    def test_max_bonds_changed(self):
        result = run_fps("--min-bonds 4", 19)
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

class TestIO(unittest.TestCase):
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

class TestArgErrors(unittest.TestCase):
    def _run(self, cmd, expect):
        msg = run_exit(cmd)
        self.assertEquals(expect in msg, True, msg)

    def test_two_fp_types(self):
        self._run("--maccs166 --path", "Cannot specify both --maccs166 and --path")

    def test_num_bits_too_small(self):
        self._run("--num-bits 0", "between 16 and 65536 bits")
        self._run("--num-bits 1", "between 16 and 65536 bits")
        self._run("--num-bits 15", "between 16 and 65536 bits")

    def test_num_bits_too_large(self):
        self._run("--num-bits 65537", "between 16 and 65536 bits")
        # Check for overflow, even though I know it won't happen in Python
        self._run("--num-bits {big}".format(big=2**32+32), "between 16 and 65536 bits")

    def test_min_bonds_too_small(self):
        self._run("--min-bonds=-1", "0 or greater")

    def test_min_bonds_larger_than_default_max_bonds(self):
        self._run("--min-bonds=6", "--max-bonds must not be smaller than --min-bonds")

    def test_min_bonds_too_large(self):
        self._run("--min-bonds=4 --max-bonds=3",
                  "--max-bonds must not be smaller than --min-bonds")

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

class TestHeaderOutput(unittest.TestCase):
    def _field(self, s, field):
        result = run(s)
        filtered = [line for line in result if line.startswith(field)]
        self.assertEquals(len(filtered), 1, result)
        return filtered[0]

    def test_software(self):
        result = self._field("", "#software")
        self.assertEquals("#software=OEGraphSim/1.0.0 (20100809)" in result, True, result)
        result = self._field("--maccs166", "#software")
        self.assertEquals("#software=OEGraphSim/1.0.0 (20100809)" in result, True, result)

    def test_type(self):
        result = self._field("", "#type")
        self.assertEquals(result,
  "#type=OpenEye-Path/1 min_bonds=0 max_bonds=5 "
  "atype=Aromaticity|AtomicNumber|Chiral|EqHalogen|FormalCharge|HvyDegree|Hybridization "
  "btype=BondOrder|Chiral")

        
    # different flags. All flags? and order
    def test_num_bits(self):
        result = self._field("--num-bits 38", "#num_bits")
        self.assertEquals(result, "#num_bits=38")
        
    def test_atype_flags(self):
        result = self._field("--atype FormalCharge|FormalCharge", "#type") + " "
        self.assertEquals(" atype=FormalCharge " in result, True, result)
        
    def test_btype_flags(self):
        result = self._field("--btype Chiral|BondOrder", "#type") + " "
        self.assertEquals(" btype=BondOrder|Chiral " in result, True, result)
    
    def test_pipe_or_comma(self):
        result = self._field("--atype HvyDegree,FormalCharge --btype Chiral,BondOrder",
                             "#type") + " "
        self.assertEquals(" atype=FormalCharge|HvyDegree " in result, True, result)
        self.assertEquals(" btype=BondOrder|Chiral " in result, True, result)
        
    
    def test_maccs_header(self):
        result = self._field("--maccs166", "#type")
        self.assertEquals(result, "#type=OpenEye-MACCS166/1")
    
        
if __name__ == "__main__":
    unittest.main()
