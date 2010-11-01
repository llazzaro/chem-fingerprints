import unittest
import sys
from cStringIO import StringIO as SIO

from openeye import oechem

from chemfp.commandline import oe2fps

real_stdout = sys.stdout

PUBCHEM_SDF = "pubchem.sdf" # Must not contain a space!

oeerrs = oechem.oeosstream()
oechem.OEThrow.SetOutputStream(oeerrs)
def _check_for_oe_errors():
    lines = oeerrs.str().splitlines()
    for line in lines:
        if line.startswith("Warning: Stereochemistry corrected on atom number"):
            continue
        raise AssertionError("Unexpected message from OEChem: {msg!r}".format(msg=line))

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

def run_fps(s, expect_length=None, source=PUBCHEM_SDF):
    result = run(s, source)
    while result[0].startswith("#"):
        del result[0]
    if expect_length is not None:
        assert len(result) == expect_length, (len(result), expect_length)
    return result

def _set_bit(n):
    assert n <= 16
    bytes = [0, 0, 0]
    bytes[n//8] = 1<<(n%8)
    return "%02x%02x%02x" % tuple(bytes)

class TestMACCS(unittest.TestCase):
    def test_FP2(self):
        result = run_fps("--maccs", 7, "maccs.smi")
        # The fingerprints are constructed to test the first few bytes.
        self.assertEquals(result[0][:6], _set_bit(2))
        self.assertEquals(result[1][:6], _set_bit(3))
        self.assertEquals(result[2][:6], _set_bit(4))
        self.assertEquals(result[3][:6], _set_bit(5))
        self.assertEquals(result[4][:6], _set_bit(9))
        self.assertEquals(result[5][:6], _set_bit(10))
        self.assertEquals(result[6][:6], _set_bit(16))

if __name__ == "__main__":
    unittest.main()
