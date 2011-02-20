import unittest2
import sys

from chemfp.commandline import simsearch

from cStringIO import StringIO
import support

class SimsearchRunner(support.Runner):
    def verify_result(self, result):
        assert result[0] == "#Simsearch/1", result[0]

runner = SimsearchRunner(simsearch.main)
run = runner.run
run_split = runner.run_split


def run_split_stdin(input, cmdline, expect_length=None, source="simple.fps"):
    old_stdin = sys.stdin
    sys.stdin = StringIO(input)
    try:
        return run_split(cmdline, expect_length, source)
    finally:
        sys.stdin = old_stdin


# The values I get using gmpy are:
#    [(1.0, 'deadbeef'),
#     (0.95999999999999996, 'Deaf_Beef'),
#     (0.83999999999999997, 'DEADdead'),
#     (0.23999999999999999, 'several'),
#     (0.041666666666666664, 'bit1'),
#     (0.040000000000000001, 'two_bits'),
#     (0.0, 'zeros')]

class TestOptions(unittest2.TestCase):
    def test_default(self):
        header, lines = run_split("--hex-query deadbeef", 1, "simple.fps")
        self.assertEquals(header,
                          {"#num_bits": "32",
                           "#software": "chemfp/0.9.1",
                           "#type": "Tanimoto k=3 threshold=0.0"})
        self.assertEquals(lines,
                          ["3 Query1 1.000 deadbeef 0.960 Deaf_Beef 0.840 DEADdead"])

    def test_k_3(self):
        header, lines = run_split("--hex-query deadbeef -k 3", 1, "simple.fps")
        self.assertEquals(header,
                          {"#num_bits": "32",
                           "#software": "chemfp/0.9.1",
                           "#type": "Tanimoto k=3 threshold=0.0"})
        self.assertEquals(lines,
                          ["3 Query1 1.000 deadbeef 0.960 Deaf_Beef 0.840 DEADdead"])

    def test_k_2(self):
        header, lines = run_split("--hex-query deadbeef -k 2", 1, "simple.fps")
        self.assertEquals(header,
                          {"#num_bits": "32",
                           "#software": "chemfp/0.9.1",
                           "#type": "Tanimoto k=2 threshold=0.0"})
        self.assertEquals(lines,
                          ["2 Query1 1.000 deadbeef 0.960 Deaf_Beef"])

    def test_k_1(self):
        header, lines = run_split("--hex-query deadbeef -k 1", 1, "simple.fps")
        self.assertEquals(header,
                          {"#num_bits": "32",
                           "#software": "chemfp/0.9.1",
                           "#type": "Tanimoto k=1 threshold=0.0"})
        self.assertEquals(lines,
                          ["1 Query1 1.000 deadbeef"])

    def test_knearest_1(self):
        header, lines = run_split("--hex-query deadbeef --k-nearest 1", 1, "simple.fps")
        self.assertEquals(header,
                          {"#num_bits": "32",
                           "#software": "chemfp/0.9.1",
                           "#type": "Tanimoto k=1 threshold=0.0"})
        self.assertEquals(lines,
                          ["1 Query1 1.000 deadbeef"])

    def test_k_10(self):
        # Asked for 10 but only 7 are available
        header, lines = run_split("--hex-query deadbeef -k 10", 1, "simple.fps")
        self.assertEquals(header,
                          {"#num_bits": "32",
                           "#software": "chemfp/0.9.1",
                           "#type": "Tanimoto k=10 threshold=0.0"})
        self.assertEquals(lines,
                          ["7 Query1 1.000 deadbeef 0.960 Deaf_Beef 0.840 DEADdead "
                           "0.240 several 0.042 bit1 0.040 two_bits 0.000 zeros"])

    def test_threshold(self):
        header, lines = run_split("--hex-query deadbeef --threshold 0.9", 1, "simple.fps")
        self.assertEquals(header,
                          {"#num_bits": "32",
                           "#software": "chemfp/0.9.1",
                           "#type": "Tanimoto k=3 threshold=0.9"})
        self.assertEquals(lines,
                          ["2 Query1 1.000 deadbeef 0.960 Deaf_Beef"])

    def test_threshold_and_k(self):
        header, lines = run_split("--hex-query deadbeef -t 0.9 -k 1", 1, "simple.fps")
        self.assertEquals(header,
                          {"#num_bits": "32",
                           "#software": "chemfp/0.9.1",
                           "#type": "Tanimoto k=1 threshold=0.9"})
        self.assertEquals(lines,
                          ["1 Query1 0.960 Deaf_Beef"])
    
    def test_stdin(self):
        header, lines = run_split_stdin("deadbeef spam\n", "", 1, "simple.fps")
        self.assertEquals(header,
                          {"#num_bits": "32",
                           "#software": "chemfp/0.9.1",
                           "#type": "Tanimoto k=3 threshold=0.0"})
        self.assertEquals(lines,
                          ["3 spam 1.000 deadbeef 0.960 Deaf_Beef 0.840 DEADdead"])

    def test_stdin2(self):
        header, lines = run_split_stdin("deadbeef spam\nDEADBEEF eggs\n",
                                        "", 2, "simple.fps")
        self.assertEquals(header,
                          {"#num_bits": "32",
                           "#software": "chemfp/0.9.1",
                           "#type": "Tanimoto k=3 threshold=0.0"})
        self.assertEquals(lines,
                          ["3 spam 1.000 deadbeef 0.960 Deaf_Beef 0.840 DEADdead",
                           "3 eggs 1.000 deadbeef 0.960 Deaf_Beef 0.840 DEADdead"])

    def test_stdin3(self):
        header, lines = run_split_stdin("deadbeef spam\n87654321 countdown\n",
                                        "--threshold 0.9", 2, "simple.fps")
        self.assertEquals(header,
                          {"#num_bits": "32",
                           "#software": "chemfp/0.9.1",
                           "#type": "Tanimoto k=3 threshold=0.9"})
        self.assertEquals(lines,
                          ["2 spam 1.000 deadbeef 0.960 Deaf_Beef",
                           "0 countdown"])

class _AgainstSelf:
    def test_with_threshold(self):
        header, lines = run_split("--queries simple.fps --threshold 0.8" + self.extra_arg,
                                  7, "simple.fps")
        self.assertEquals(lines,
                          ["0 zeros",
                           "1 bit1 1.000 bit1",
                           "1 two_bits 1.000 two_bits",
                           "1 several 1.000 several",
                           "3 deadbeef 1.000 deadbeef 0.960 Deaf_Beef 0.840 DEADdead",
                           "3 DEADdead 1.000 DEADdead 0.840 deadbeef 0.808 Deaf_Beef",
                           "3 Deaf_Beef 1.000 Deaf_Beef 0.960 deadbeef 0.808 DEADdead"])

    def test_with_count_and_threshold(self):
        header, lines = run_split("--queries simple.fps --threshold 0.8 --count" + self.extra_arg,
                                  7, "simple.fps")
        self.assertEquals(lines,
                          ["0 zeros",
                           "1 bit1",
                           "1 two_bits",
                           "1 several",
                           "3 deadbeef",
                           "3 DEADdead",
                           "3 Deaf_Beef"])

class TestAgainstSelf(unittest2.TestCase, _AgainstSelf):
    extra_arg = ""

class TestAgainstSelfInMemory(unittest2.TestCase, _AgainstSelf):
    extra_arg = " --in-memory"

class TestAgainstSelfFileScan(unittest2.TestCase, _AgainstSelf):
    extra_arg = " --file-scan"


                          
        
if __name__ == "__main__":
    unittest2.main()
