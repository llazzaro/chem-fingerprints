import unittest
import sys

from chemfp.commandline import ob2fps

real_stdout = sys.stdout

PUBCHEM_SDF = "pubchem.sdf"

def run(s):
    args = s.split()
    try:
        sys.stdin = open(PUBCHEM_SDF)
        sys.stdout = stdout = SIO()
        sdf2fps.main(args)
    finally:
        sys.stdout = real_stdout
        sys.stdin = real_stdin
    result = stdout.getvalue().splitlines()
    if result:
        assert result[0] == "#FPS1"
    return result

def run_fps(s, expect_length=None):
    result = run(s)
    while result[0].startswith("#"):
        del result[0]
    if expect_length is not None:
        assert len(result) == expect_length
    return result


class TestFingerprintTypes(unittest.TestCase):
    def test_FP2(self):
        result = run("--FP2")
        self.assertEquals(len(result), 11)

if __name__ == "__main__":
    unittest.main()
