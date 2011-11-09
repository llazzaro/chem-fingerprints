import sys
import unittest2
from cStringIO import StringIO as SIO

from chemfp.commandline import sdf2fps

import support

real_stdin = sys.stdin
real_stdout = sys.stdout
real_stderr = sys.stderr

DECODER_SDF = support.fullpath("decoder.sdf")

def run(s):
    args = s.split()
    try:
        sys.stdin = open(DECODER_SDF)
        sys.stdout = stdout = SIO()
        sdf2fps.main(args)
    finally:
        sys.stdout = real_stdout
        sys.stdin = real_stdin
    result = stdout.getvalue().splitlines()
    if result:
        assert result[0] == "#FPS1"
    return result

def run_failure(s):
    sys.stderr = stderr = SIO()
    try:
        try:
            run(s)
        except SystemExit:
            pass
        else:
            raise AssertionError("should have exited: %r" % (s,))
    finally:
        sys.stderr = real_stderr
    return stderr.getvalue()

def run_warning(s):
    sys.stderr = stderr = SIO()
    try:
        try:
            run(s)
        finally:
            sys.stderr = real_stderr
    except SystemExit, e:
        raise AssertionError("unexpected SystemExit: %r and %r" % (e, stderr.getvalue()))
    return stderr.getvalue()

def run_fps(s, expect_length=None):
    result = run(s)
    while result[0].startswith("#"):
        del result[0]
    if expect_length is not None:
        assert len(result) == expect_length
    return result


class TestDecoderFlags(unittest2.TestCase):
    def test_cactvs(self):
        result = run_fps("--cactvs --fp-tag PUBCHEM_CACTVS_SUBSKEYS")
        self.assertEquals(result, ["07de8d002000000000000000000000000080060000000c000000000000000080030000f8401800000030508379344c014956000055c0a44e2a0049200084e140581f041d661b10064483cb0f2925100619001393e10001007000000000008000000000000000400000000000000000\t9425004",
                                   "07de0d000000000000000000000000000080460300000c0000000000000000800f0000780038000000301083f920cc09695e0800d5c0e44e6e00492190844145dc1f841d261911164d039b8f29251026b9401313e0ec01007000000000000000000000000000000000000000000000\t9425009"])

    def test_binary40(self):
        result = run_fps("--binary --fp-tag binary40", 2)
        self.assertEquals(result[0], "000500c000\t9425004")
        self.assertEquals(result[1], "00fab75300\t9425009")

    def test_binary_msb40(self):
        result = run_fps("--binary-msb --fp-tag binary40", 2)
        self.assertEquals(result[0], "000300a000\t9425004")
        self.assertEquals(result[1], "00caed5f00\t9425009")

    def test_binary3(self):
        result = run_fps("--binary --fp-tag binary3", 2)
        self.assertEquals(result[0], "04\t9425004")
        self.assertEquals(result[1], "03\t9425009")

    def test_binary_msb3(self):
        result = run_fps("--binary-msb --fp-tag binary3", 2)
        self.assertEquals(result[0], "01\t9425004")
        self.assertEquals(result[1], "06\t9425009")

    def test_binary8(self):
        result = run_fps("--binary --fp-tag binary8", 2)
        self.assertEquals(result[0], "76\t9425004")
        self.assertEquals(result[1], "bc\t9425009")

    def test_binary_msb8(self):
        result = run_fps("--binary-msb --fp-tag binary8", 2)
        self.assertEquals(result[0], "6e\t9425004")
        self.assertEquals(result[1], "3d\t9425009")


    def test_binary17(self):
        result = run_fps("--binary --fp-tag binary17", 2)
        self.assertEquals(result[0], "38b701\t9425004")
        self.assertEquals(result[1], "489d01\t9425009")

    def test_binary_msb17(self):
        result = run_fps("--binary-msb --fp-tag binary17", 2)
        self.assertEquals(result[0], "db3900\t9425004")
        self.assertEquals(result[1], "732500\t9425009")

    def test_hex2(self):
        result = run_fps("--hex --fp-tag hex2", 2)
        self.assertEquals(result[0], "ab\t9425004")
        self.assertEquals(result[1], "01\t9425009")
        
    def test_hex_lsb2(self):
        # 0xab == 0b10101011
        # 10101011 with LSB first is 5 d => "d5"
        # 0x01 == 0b00000001 => 80 when in LSB first
        result = run_fps("--hex-lsb --fp-tag hex2", 2)
        self.assertEquals(result[0], "d5\t9425004")
        self.assertEquals(result[1], "80\t9425009")

    def test_hex_msb2(self):
        # With 2 nibbles the result is the same as hex
        result = run_fps("--hex-msb --fp-tag hex2", 2)
        self.assertEquals(result[0], "ab\t9425004")
        self.assertEquals(result[1], "01\t9425009")

    def test_hex16(self):
        result = run_fps("--hex --fp-tag hex16", 2)
        self.assertEquals(result[0], "0123456789abcdef\t9425004")
        self.assertEquals(result[1], "abcdef0123456789\t9425009")
        
    def test_hex_lsb16(self):
        result = run_fps("--hex-lsb --fp-tag hex16", 2)
        # 0123456789abcdef in LSB form => 
        # 084c2a6e195d3b7f when nibbles bits are in MSB form but nibbles are LSB
        # 80 c4 a2 e6 91 d5 b3 f7 when byte bits are in MSB and bytes are LSB
        self.assertEquals(result[0], "80c4a2e691d5b3f7\t9425004")
        # abcdef0123456789 in LSB form =>
        # 5d3b7f084c2a6e19 =>
        # d5 b3 f7 80 c4 a2 e6 91
        self.assertEquals(result[1], "d5b3f780c4a2e691\t9425009")

    def test_hex_msb16(self):
        # Just a bit of reordering
        result = run_fps("--hex-msb --fp-tag hex16", 2)
        self.assertEquals(result[0], "efcdab8967452301\t9425004")
        self.assertEquals(result[1], "8967452301efcdab\t9425009")
        
    def test_base64_16(self):
        result = run_fps("--base64 --fp-tag base64_16", 2)
        self.assertEquals(result[0], "Greetings, human".encode("hex") + "\t9425004")
        self.assertEquals(result[1], "blahblahspamblah".encode("hex") + "\t9425009")


    def test_bad_decoding(self):
        msg = run_warning("--base64 --fp-tag binary17 --errors report")
        self.assertIn("Could not base64 decode 'binary17' value", msg)
        self.assertIn("Skipping.", msg)

class TestBitSizes(unittest2.TestCase):
    def test_exact_fingerprint_bits(self):
        result = run("--binary --fp-tag binary3")
        self.assertIn("#num_bits=3", result)
        
    def test_user_bits_match_fingerprint_bits(self):
        result = run("--binary --fp-tag binary3 --num-bits 3")
        self.assertIn("#num_bits=3", result)
        self.assertIn("04\t9425004", result)
        self.assertIn("03\t9425009", result)

    def test_user_bits_disagree_with_fingerprint_bits(self):
        errmsg = run_failure("--binary --fp-tag binary3 --num-bits 2")
        self.assertIn("has 3 bits", errmsg)
        self.assertIn(" 2", errmsg)

    def test_implied_from_fingerprint_bytes(self):
        result = run("--hex --fp-tag hex2")
        self.assertIn("#num_bits=8", result)

    def test_user_bits_matches_fingerprint_bytes(self):
        result = run("--hex --fp-tag hex2 --num-bits 8")
        self.assertIn("#num_bits=8", result)

    def test_user_bits_too_large_for_bytes(self):
        result = run_failure("--hex --fp-tag hex2 --num-bits 9")
        self.assertIn("1 <= num-bits <= 8, not 9", result)

    def test_user_bits_acceptably_smaller_than_bytes(self):
        result = run("--hex --fp-tag hex2 --num-bits 6")
        self.assertIn("#num_bits=6", result)

    def test_user_bits_too_much_smaller_than_bytes(self):
        result = run_failure("--hex --fp-tag hex16 --num-bits 56")
        self.assertIn("57 <= num-bits <= 64, not 56", result)

class TestTitleProcessing(unittest2.TestCase):
    def test_title_from_title_tag(self):
        result = run("--hex --fp-tag hex2 --id-tag binary3")
        self.assertIn("ab\t001", result)

    def test_missing_title_from_title_line(self):
        warning = run_warning("--hex --fp-tag hex2 --id-tag FAKE_TITLE --errors report")
        self.assertIn("Missing id tag 'FAKE_TITLE' in the record starting at line 146 of ", warning)
        self.assertIn("decoder.sdf", warning)
        self.assertIn("title='9425009'", warning)
        self.assertIn("Skipping.", warning)

    def test_missing_all_titles(self):
        warning = run_warning("--hex --fp-tag hex2 --id-tag DOES_NOT_EXIST --errors report")
        self.assertIn("Missing id tag 'DOES_NOT_EXIST'", warning)
        self.assertIn("line 1 of", warning)
        self.assertIn("line 146 of", warning)

class TestShortcuts(unittest2.TestCase):
    def test_pubchem(self):
        result = run("--pubchem")
        self.assertIn("#num_bits=881", result)
        self.assertIn("#software=CACTVS/unknown", result)
        self.assertIn("#type=CACTVS-E_SCREEN/1.0 extended=2", result)
        self.assertIn("07de8d002000000000000000000000000080060000000c000000000000000080030000f8401800000030508379344c014956000055c0a44e2a0049200084e140581f041d661b10064483cb0f2925100619001393e10001007000000000008000000000000000400000000000000000\t9425004", result)
        self.assertIn("07de0d000000000000000000000000000080460300000c0000000000000000800f0000780038000000301083f920cc09695e0800d5c0e44e6e00492190844145dc1f841d261911164d039b8f29251026b9401313e0ec01007000000000000000000000000000000000000000000000\t9425009", result)

class TestBadArgs(unittest2.TestCase):
    def test_missing_fp_tag(self):
        msg = run_failure("")
        self.assertIn("argument --fp-tag is required", msg)

    def test_num_bits_positive(self):
        msg = run_failure("--fp-tag SPAM --num-bits 0")
        self.assertIn("--num-bits must be a positive integer", msg)
        msg = run_failure("--fp-tag SPAM --num-bits -1")
        self.assertIn("--num-bits must be a positive integer", msg)

    def test_bad_char(self):
        msg = run_failure("--fp-tag SPAM --software this\bthat")
        self.assertIn("--software", msg)
        self.assertIn("'\\x08'", msg)
        
        
if __name__ == "__main__":
    unittest2.main()
