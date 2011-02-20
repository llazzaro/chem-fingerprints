import unittest2

from chemfp.commandline import simsearch

import support

class SimsearchRunner(support.Runner):
    def verify_result(self, result):
        assert result[0] == "#Simsearch/1", result[0]

runner = SimsearchRunner(simsearch.main)
run = runner.run
run_split = runner.run_split

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
        header, lines = run_split("--hex-query deadbeef", None, "simple.fps")
        self.assertEquals(header,
                          {"#num_bits": "32",
                           "#software": "chemfp/0.9.1",
                           "#type": "Tanimoto k=3 threshold=0.0"})
        self.assertEquals(lines,
                          ["3 Query1 1.000 deadbeef 0.960 Deaf_Beef 0.840 DEADdead"])

    def test_k_3(self):
        header, lines = run_split("--hex-query deadbeef -k 3", None, "simple.fps")
        self.assertEquals(header,
                          {"#num_bits": "32",
                           "#software": "chemfp/0.9.1",
                           "#type": "Tanimoto k=3 threshold=0.0"})
        self.assertEquals(lines,
                          ["3 Query1 1.000 deadbeef 0.960 Deaf_Beef 0.840 DEADdead"])

    def test_k_2(self):
        header, lines = run_split("--hex-query deadbeef -k 2", None, "simple.fps")
        self.assertEquals(header,
                          {"#num_bits": "32",
                           "#software": "chemfp/0.9.1",
                           "#type": "Tanimoto k=2 threshold=0.0"})
        self.assertEquals(lines,
                          ["2 Query1 1.000 deadbeef 0.960 Deaf_Beef"])

    def test_k_1(self):
        header, lines = run_split("--hex-query deadbeef -k 1", None, "simple.fps")
        self.assertEquals(header,
                          {"#num_bits": "32",
                           "#software": "chemfp/0.9.1",
                           "#type": "Tanimoto k=1 threshold=0.0"})
        self.assertEquals(lines,
                          ["1 Query1 1.000 deadbeef"])

    def test_knearest_1(self):
        header, lines = run_split("--hex-query deadbeef --k-nearest 1", None, "simple.fps")
        self.assertEquals(header,
                          {"#num_bits": "32",
                           "#software": "chemfp/0.9.1",
                           "#type": "Tanimoto k=1 threshold=0.0"})
        self.assertEquals(lines,
                          ["1 Query1 1.000 deadbeef"])

    def test_k_10(self):
        # Asked for 10 but only 7 are available
        header, lines = run_split("--hex-query deadbeef -k 10", None, "simple.fps")
        self.assertEquals(header,
                          {"#num_bits": "32",
                           "#software": "chemfp/0.9.1",
                           "#type": "Tanimoto k=10 threshold=0.0"})
        self.assertEquals(lines,
                          ["7 Query1 1.000 deadbeef 0.960 Deaf_Beef 0.840 DEADdead "
                           "0.240 several 0.042 bit1 0.040 two_bits 0.000 zeros"])

    def test_threshold(self):
        header, lines = run_split("--hex-query deadbeef --threshold 0.9", None, "simple.fps")
        self.assertEquals(header,
                          {"#num_bits": "32",
                           "#software": "chemfp/0.9.1",
                           "#type": "Tanimoto k=3 threshold=0.9"})
        self.assertEquals(lines,
                          ["2 Query1 1.000 deadbeef 0.960 Deaf_Beef"])
    
if __name__ == "__main__":
    unittest2.main()
