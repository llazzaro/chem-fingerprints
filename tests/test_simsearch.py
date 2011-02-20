import unittest2

from chemfp.commandline import simsearch

import support

class SimsearchRunner(support.Runner):
    def verify_result(self, result):
        assert result[0] == "#Simsearch/1", result[0]

runner = SimsearchRunner(simsearch.main)
run = runner.run
run_split = runner.run_split

class TestSimple(unittest2.TestCase):
    def test_deadbeef(self):
        header, lines = run_split("--hex-query deadbeef", None, "simple.fps")
        self.assertEquals(header,
                          {"#num_bits": "32",
                           "#software": "chemfp/0.9.1",
                           "#type": "Tanimoto k=3 threshold=0.0"})
        self.assertEquals(lines,
                          ["3 Query1 1.000 deadbeef 0.960 Deaf_Beef 0.840 DEADdead"])
        
if __name__ == "__main__":
    unittest2.main()
