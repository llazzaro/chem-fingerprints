import unittest2
import sys
import gzip
from cStringIO import StringIO

import chemfp
from chemfp.commandline import simsearch

SOFTWARE = "chemfp/" + chemfp.__version__

import support

class SimsearchRunner(support.Runner):
    def verify_result(self, result):
        assert result[0] == "#Simsearch/1", result[0]
class CountRunner(support.Runner):
    def verify_result(self, result):
        assert result[0] == "#Count/1", result[0]

runner = SimsearchRunner(simsearch.main)
run = runner.run
run_split = runner.run_split
run_exit = runner.run_exit

count_runner = CountRunner(simsearch.main)
count_run_split = count_runner.run_split
count_run_exit = count_runner.run_exit

SIMPLE_FPS = support.fullpath("simple.fps")
SIMPLE_FPS_GZ = support.fullpath("simple.fps.gz")

def run_split_stdin(input, cmdline, expect_length=None, source=SIMPLE_FPS):
    old_stdin = sys.stdin
    sys.stdin = StringIO(input)
    try:
        return run_split(cmdline, expect_length, source)
    finally:
        sys.stdin = old_stdin

def gzip_compress(s):
    f = StringIO()
    g = gzip.GzipFile(fileobj=f, mode="w")
    g.write(s)
    g.close()
    return f.getvalue()


# The values I get using gmpy are:
#    [(1.0, 'deadbeef'),
#     (0.95999999999999996, 'Deaf Beef'),
#     (0.83999999999999997, 'DEADdead'),
#     (0.23999999999999999, 'several'),
#     (0.041666666666666664, 'bit1'),
#     (0.040000000000000001, 'two_bits'),
#     (0.0, 'zeros')]

class TestOptions(unittest2.TestCase):
    def test_default(self):
        header, lines = run_split("--hex-query deadbeef -t 0.1", 1, SIMPLE_FPS)
        self.assertIn("simple.fps", header.pop("#targets"))
        self.assertEquals(header,
                          {"#num_bits": "32",
                           "#software": SOFTWARE,
                           "#type": "Tanimoto k=all threshold=0.1"})
        self.assertEquals(len(lines), 1, lines)
        fields = lines[0].split("\t")
        self.assertEquals(fields[:2], ["4", "Query1"])
        hits = zip(fields[2::2], fields[3::2])
        hits.sort()
        self.assertEquals(hits, [("DEADdead", "0.840"), ("Deaf Beef", "0.960"),
                                 ("deadbeef", "1.000"), ('several', '0.240')])

    def test_k_3(self):
        header, lines = run_split("--hex-query deadbeef -k 3 --threshold 0.8", 1, SIMPLE_FPS)
        self.assertIn("simple.fps", header.pop("#targets"))
        self.assertEquals(header,
                          {"#num_bits": "32",
                           "#software": SOFTWARE,
                           "#type": "Tanimoto k=3 threshold=0.8"})
        self.assertEquals(lines,
                          ["3\tQuery1\tdeadbeef\t1.000\tDeaf Beef\t0.960\tDEADdead\t0.840"])

    def test_k_2(self):
        header, lines = run_split("--hex-query deadbeef -k 2 --threshold 0.9", 1, SIMPLE_FPS)
        self.assertIn("simple.fps", header.pop("#targets"))
        self.assertEquals(header,
                          {"#num_bits": "32",
                           "#software": SOFTWARE,
                           "#type": "Tanimoto k=2 threshold=0.9"})
        self.assertEquals(lines,
                          ["2\tQuery1\tdeadbeef\t1.000\tDeaf Beef\t0.960"])

    def test_k_1(self):
        header, lines = run_split("--hex-query deadbeef -k 1 -t 0.0", 1, SIMPLE_FPS)
        self.assertIn("simple.fps", header.pop("#targets"))
        self.assertEquals(header,
                          {"#num_bits": "32",
                           "#software": SOFTWARE,
                           "#type": "Tanimoto k=1 threshold=0.0"})
        self.assertEquals(lines,
                          ["1\tQuery1\tdeadbeef\t1.000"])

    def test_knearest_1(self):
        header, lines = run_split("--hex-query deadbeef --k-nearest 1", 1, SIMPLE_FPS)
        self.assertIn("simple.fps", header.pop("#targets"))
        self.assertEquals(header,
                          {"#num_bits": "32",
                           "#software": SOFTWARE,
                           "#type": "Tanimoto k=1 threshold=0.0"})
        self.assertEquals(lines,
                          ["1\tQuery1\tdeadbeef\t1.000"])

    def test_k_10(self):
        # Asked for 10 but only 7 are available
        header, lines = run_split("--hex-query deadbeef -k 10 --threshold 0.0", 1, SIMPLE_FPS)
        self.assertIn("simple.fps", header.pop("#targets"))
        self.assertEquals(header,
                          {"#num_bits": "32",
                           "#software": SOFTWARE,
                           "#type": "Tanimoto k=10 threshold=0.0"})
        self.assertEquals(lines,
                          ["7\tQuery1\tdeadbeef\t1.000\tDeaf Beef\t0.960\tDEADdead\t0.840\t"
                           "several\t0.240\tbit1\t0.042\ttwo_bits\t0.040\tzeros\t0.000"])

    def test_knearest_all(self):
        header, lines = run_split("--hex-query deadbeef --k-nearest all", 1, SIMPLE_FPS)
        self.assertIn("simple.fps", header.pop("#targets"))
        self.assertEquals(header,
                          {"#num_bits": "32",
                           "#software": SOFTWARE,
                           "#type": "Tanimoto k=all threshold=0.0"})
        self.assertEquals(lines,
                          ['7\tQuery1\tzeros\t0.000\tbit1\t0.042\ttwo_bits\t0.040\tseveral\t0.240\tdeadbeef\t1.000\tDEADdead\t0.840\tDeaf Beef\t0.960'])

    def test_threshold(self):
        header, lines = run_split("--hex-query deadbeef --threshold 0.9", 1, SIMPLE_FPS)
        self.assertIn("simple.fps", header.pop("#targets"))
        self.assertEquals(header,
                          {"#num_bits": "32",
                           "#software": SOFTWARE,
                           "#type": "Tanimoto k=all threshold=0.9"})
        self.assertEquals(lines,
                          ["2\tQuery1\tdeadbeef\t1.000\tDeaf Beef\t0.960"])

    def test_threshold_and_k(self):
        header, lines = run_split("--hex-query deadbeef -t 0.9 -k 1", 1, SIMPLE_FPS)
        self.assertIn("simple.fps", header.pop("#targets"))
        self.assertEquals(header,
                          {"#num_bits": "32",
                           "#software": SOFTWARE,
                           "#type": "Tanimoto k=1 threshold=0.9"})
        self.assertEquals(lines,
                          ["1\tQuery1\tdeadbeef\t1.000"])

    def test_threshold_and_k_all(self):
        header, lines = run_split("--hex-query deadbeef --threshold 0.9 --k-nearest all", 1, SIMPLE_FPS)
        self.assertIn("simple.fps", header.pop("#targets"))
        self.assertEquals(header,
                          {"#num_bits": "32",
                           "#software": SOFTWARE,
                           "#type": "Tanimoto k=all threshold=0.9"})
        self.assertEquals(lines,
                          ["2\tQuery1\tdeadbeef\t1.000\tDeaf Beef\t0.960"])

    
    def test_stdin(self):
        header, lines = run_split_stdin("deadbeef\tspam\n", "", 1, SIMPLE_FPS)
        self.assertIn("simple.fps", header.pop("#targets"))
        self.assertEquals(header,
                          {"#num_bits": "32",
                           "#software": SOFTWARE,
                           "#type": "Tanimoto k=3 threshold=0.7"})
        self.assertEquals(lines,
                          ["3\tspam\tdeadbeef\t1.000\tDeaf Beef\t0.960\tDEADdead\t0.840"])

    def test_stdin2(self):
        header, lines = run_split_stdin("deadbeef\tspam\nDEADBEEF\teggs\n",
                                        "-k 3 --threshold 0.6", 2, SIMPLE_FPS)
        self.assertIn("simple.fps", header.pop("#targets"))
        self.assertEquals(header,
                          {"#num_bits": "32",
                           "#software": SOFTWARE,
                           "#type": "Tanimoto k=3 threshold=0.6"})
        self.assertEquals(lines,
                          ["3\tspam\tdeadbeef\t1.000\tDeaf Beef\t0.960\tDEADdead\t0.840",
                           "3\teggs\tdeadbeef\t1.000\tDeaf Beef\t0.960\tDEADdead\t0.840"])

    def test_stdin3(self):
        header, lines = run_split_stdin("deadbeef\tspam\n87654321\tcountdown\n",
                                        "-k 3 --threshold 0.9", 2, SIMPLE_FPS)
        self.assertIn("simple.fps", header.pop("#targets"))
        self.assertEquals(header,
                          {"#num_bits": "32",
                           "#software": SOFTWARE,
                           "#type": "Tanimoto k=3 threshold=0.9"})
        self.assertEquals(lines,
                          ["2\tspam\tdeadbeef\t1.000\tDeaf Beef\t0.960",
                           "0\tcountdown"])

    def test_in_format_deprecated(self):
        s = gzip_compress("deadbeef\tspam\n")
        # You should use use "--query-format" instead of "--in".
        header, lines = run_split_stdin(s, "--in fps.gz", 1, SIMPLE_FPS)
        self.assertIn("simple.fps", header.pop("#targets"))
        self.assertEquals(header,
                          {"#num_bits": "32",
                           "#software": SOFTWARE,
                           "#type": "Tanimoto k=3 threshold=0.7"})
        self.assertEquals(lines,
                          ["3\tspam\tdeadbeef\t1.000\tDeaf Beef\t0.960\tDEADdead\t0.840"])

    def test_query_format_stdin(self):
        s = gzip_compress("deadbeef\tspam\n")
        header, lines = run_split_stdin(s, "--query-format fps.gz", 1, SIMPLE_FPS)
        self.assertIn("simple.fps", header.pop("#targets"))
        self.assertEquals(header,
                          {"#num_bits": "32",
                           "#software": SOFTWARE,
                           "#type": "Tanimoto k=3 threshold=0.7"})
        self.assertEquals(lines,
                          ["3\tspam\tdeadbeef\t1.000\tDeaf Beef\t0.960\tDEADdead\t0.840"])

    def test_query_format_fps_gzip_but_actually_fps(self):
        errmsg = run_exit(["--queries", SIMPLE_FPS, "--query-format", "fps.gz"], SIMPLE_FPS)
        self.assertIn("Not a gzipped file", errmsg)

    def test_query_format_fps_but_actually_fps_gz(self):
        errmsg = run_exit(["--queries", SIMPLE_FPS_GZ, "--query-format", "fps"], SIMPLE_FPS)
        self.assertIn("Line must end with a newline character at line 1 of", errmsg)

    def test_query_format_unknown(self):
        errmsg = run_exit(["--queries", SIMPLE_FPS, "--query-format", "f2s"], SIMPLE_FPS)
        self.assertIn("Unknown fingerprint format 'f2s'", errmsg)

    def test_query_format_compression_unknown(self):
        errmsg = run_exit(["--queries", SIMPLE_FPS, "--query-format", "fps.33"], SIMPLE_FPS)
        self.assertIn("Unsupported compression in format 'fps.33'", errmsg)

    def test_target_format_fps(self):
        header, lines = run_split_stdin("deadbeef\tspam\n", "--target-format fps", 1, SIMPLE_FPS)
        self.assertIn("simple.fps", header.pop("#targets"))
        self.assertEquals(header,
                          {"#num_bits": "32",
                           "#software": SOFTWARE,
                           "#type": "Tanimoto k=3 threshold=0.7"})
        self.assertEquals(lines,
                          ["3\tspam\tdeadbeef\t1.000\tDeaf Beef\t0.960\tDEADdead\t0.840"])

    def test_target_format_fps_gz(self):
        header, lines = run_split_stdin("deadbeef\tspam\n", "--target-format fps.gz", 1, SIMPLE_FPS_GZ)
        self.assertIn("simple.fps", header.pop("#targets"))
        self.assertEquals(header,
                          {"#num_bits": "32",
                           "#software": SOFTWARE,
                           "#type": "Tanimoto k=3 threshold=0.7"})
        self.assertEquals(lines,
                          ["3\tspam\tdeadbeef\t1.000\tDeaf Beef\t0.960\tDEADdead\t0.840"])

    def test_target_format_fps_but_actually_fps_gz(self):
        errmsg = run_exit(["--queries", SIMPLE_FPS, "--target-format", "fps"], SIMPLE_FPS_GZ)
        self.assertIn("Cannot open targets file: Line must end with a newline character at line 1 of", errmsg)

    def test_target_format_unknown(self):
        errmsg = run_exit(["--queries", SIMPLE_FPS, "--target-format", "fp3"], SIMPLE_FPS_GZ)
        self.assertIn("Unknown fingerprint format 'fp3'", errmsg)

    def test_target_format_compression_unknown(self):
        errmsg = run_exit(["--queries", SIMPLE_FPS, "--target-format", "fps.33"], SIMPLE_FPS_GZ)
        self.assertIn("Unsupported compression in format 'fps.33'", errmsg)
        

class _AgainstSelf:
    def test_with_threshold(self):
        header, lines = run_split(
            ["--queries", SIMPLE_FPS, "-k", "3", "--threshold", "0.8"] + self.extra_arg,
            7, SIMPLE_FPS)
        self.assertEquals(lines,
                          ["0\tzeros",
                           "1\tbit1\tbit1\t1.000",
                           "1\ttwo_bits\ttwo_bits\t1.000",
                           "1\tseveral\tseveral\t1.000",
                           "3\tdeadbeef\tdeadbeef\t1.000\tDeaf Beef\t0.960\tDEADdead\t0.840",
                           "3\tDEADdead\tDEADdead\t1.000\tdeadbeef\t0.840\tDeaf Beef\t0.808",
                           "3\tDeaf Beef\tDeaf Beef\t1.000\tdeadbeef\t0.960\tDEADdead\t0.808"])

    def test_with_count_and_threshold(self):
        header, lines = count_run_split(
            ["--queries", SIMPLE_FPS, "--threshold", "0.8", "--count"] + self.extra_arg,
            7, SIMPLE_FPS)
        self.assertEquals(lines,
                          ["0\tzeros",
                           "1\tbit1",
                           "1\ttwo_bits",
                           "1\tseveral",
                           "3\tdeadbeef",
                           "3\tDEADdead",
                           "3\tDeaf Beef"])

class TestAgainstSelf(unittest2.TestCase, _AgainstSelf):
    extra_arg = []

class TestAgainstSelfInMemory(unittest2.TestCase, _AgainstSelf):
    extra_arg = ["--memory"]

class TestAgainstSelfFileScan(unittest2.TestCase, _AgainstSelf):
    extra_arg = ["--scan"]


class TestNxN(unittest2.TestCase):
    def test_default(self):
        header, lines = run_split("--NxN", 7, SIMPLE_FPS)
        self.assertIn("simple.fps", header.pop("#targets"))
        self.assertEquals(header,
                          {"#num_bits": "32",
                          "#software": SOFTWARE,
                          "#type": "Tanimoto k=3 threshold=0.7 NxN=full"})
        self.assertEquals(len(lines), 7, lines)
        self.assertEquals(lines, ['0\tzeros',
                                  '0\tbit1',
                                  '0\ttwo_bits',
                                  '0\tseveral',
                                  '2\tdeadbeef\tDeaf Beef\t0.960\tDEADdead\t0.840',
                                  '2\tDEADdead\tdeadbeef\t0.840\tDeaf Beef\t0.808',
                                  '2\tDeaf Beef\tdeadbeef\t0.960\tDEADdead\t0.808'])
    def test_specify_default_values(self):
        header, lines = run_split("--NxN -k 3 --threshold 0.7", 7, SIMPLE_FPS)
        self.assertIn("simple.fps", header.pop("#targets"))
        self.assertEquals(header,
                          {"#num_bits": "32",
                          "#software": SOFTWARE,
                          "#type": "Tanimoto k=3 threshold=0.7 NxN=full"})
        self.assertEquals(len(lines), 7, lines)
        self.assertEquals(lines, ['0\tzeros',
                                  '0\tbit1',
                                  '0\ttwo_bits',
                                  '0\tseveral',
                                  '2\tdeadbeef\tDeaf Beef\t0.960\tDEADdead\t0.840',
                                  '2\tDEADdead\tdeadbeef\t0.840\tDeaf Beef\t0.808',
                                  '2\tDeaf Beef\tdeadbeef\t0.960\tDEADdead\t0.808'])

    def test_k_2(self):
        # This sets the theshold to 0.0
        header, lines = run_split("--NxN -k 2 ", 7, SIMPLE_FPS)
        self.assertIn("simple.fps", header.pop("#targets"))
        self.assertEquals(header,
                          {"#num_bits": "32",
                          "#software": SOFTWARE,
                          "#type": "Tanimoto k=2 threshold=0.0 NxN=full"})
        self.assertEquals(len(lines), 7, lines)
        self.assertEquals(lines, ['0\tzeros',
                                  '2\tbit1\ttwo_bits\t0.500\tseveral\t0.143',
                                  '2\ttwo_bits\tbit1\t0.500\tseveral\t0.286',
                                  '2\tseveral\ttwo_bits\t0.286\tDEADdead\t0.261',
                                  '2\tdeadbeef\tDeaf Beef\t0.960\tDEADdead\t0.840',
                                  '2\tDEADdead\tdeadbeef\t0.840\tDeaf Beef\t0.808',
                                  '2\tDeaf Beef\tdeadbeef\t0.960\tDEADdead\t0.808'])

    def test_threshold(self):
        header, lines = run_split("--NxN --threshold 0.5 ", 7, SIMPLE_FPS)
        self.assertIn("simple.fps", header.pop("#targets"))
        self.assertEquals(header,
                          {"#num_bits": "32",
                          "#software": SOFTWARE,
                          "#type": "Tanimoto k=all threshold=0.5 NxN=full"})
        self.assertEquals(len(lines), 7, lines)
        self.assertEquals(lines, ['0\tzeros',
                                  '1\tbit1\ttwo_bits\t0.500',
                                  '1\ttwo_bits\tbit1\t0.500',
                                  '0\tseveral',
                                  # The order here is implementation dependent...
                                  '2\tdeadbeef\tDeaf Beef\t0.960\tDEADdead\t0.840',
                                  '2\tDEADdead\tdeadbeef\t0.840\tDeaf Beef\t0.808',
                                  #'2\tDeaf Beef\tdeadbeef\t0.960\tDEADdead\t0.808',
                                  '2\tDeaf Beef\tDEADdead\t0.808\tdeadbeef\t0.960',
            ])

    def test_count_with_threshold(self):
        header, lines = count_run_split("--NxN --count --threshold 0.5 ", 7, SIMPLE_FPS)
        self.assertIn("simple.fps", header.pop("#targets"))
        self.assertEquals(header,
                          {"#num_bits": "32",
                          "#software": SOFTWARE,
                          "#type": "Count threshold=0.5 NxN=full"})
        self.assertEquals(len(lines), 7, lines)
        self.assertEquals(lines, ['0\tzeros',
                                  '1\tbit1',
                                  '1\ttwo_bits',
                                  '0\tseveral',
                                  '2\tdeadbeef',
                                  '2\tDEADdead',
                                  '2\tDeaf Beef',
            ])

    def test_count_with_default_threshold(self):
        header, lines = count_run_split("--NxN --count", 7, SIMPLE_FPS)
        self.assertIn("simple.fps", header.pop("#targets"))
        self.assertEquals(header,
                          {"#num_bits": "32",
                          "#software": SOFTWARE,
                          "#type": "Count threshold=0.7 NxN=full"})
        self.assertEquals(len(lines), 7, lines)
        self.assertEquals(lines, ['0\tzeros',
                                  '0\tbit1',
                                  '0\ttwo_bits',
                                  '0\tseveral',
                                  '2\tdeadbeef',
                                  '2\tDEADdead',
                                  '2\tDeaf Beef',
            ])

    def test_threshold_with_low_batch_size(self):
        header, lines = run_split("--NxN --threshold 0.5 --batch-size 1", 7, SIMPLE_FPS)
        self.assertIn("simple.fps", header.pop("#targets"))
        self.assertEquals(header,
                          {"#num_bits": "32",
                          "#software": SOFTWARE,
                          "#type": "Tanimoto k=all threshold=0.5 NxN=full"})
        self.assertEquals(len(lines), 7, lines)
        self.assertEquals(lines, ['0\tzeros',
                                  '1\tbit1\ttwo_bits\t0.500',
                                  '1\ttwo_bits\tbit1\t0.500',
                                  '0\tseveral',
                                  # The order here is implementation dependent...
                                  '2\tdeadbeef\tDeaf Beef\t0.960\tDEADdead\t0.840',
                                  '2\tDEADdead\tdeadbeef\t0.840\tDeaf Beef\t0.808',
                                  #'2\tDeaf Beef\tdeadbeef\t0.960\tDEADdead\t0.808',
                                  '2\tDeaf Beef\tDEADdead\t0.808\tdeadbeef\t0.960',
            ])
    def test_knearest_with_low_batch_size(self):
        # This is the same as test_k_2 but with a batch-size of 1.
        # This tests a bug where I wasn't incrementing the offset
        # to the start of each batch location in the results.
        header, lines = run_split("--NxN -k 2 --batch-size 1", 7, SIMPLE_FPS)
        self.assertIn("simple.fps", header.pop("#targets"))
        self.assertEquals(header,
                          {"#num_bits": "32",
                          "#software": SOFTWARE,
                          "#type": "Tanimoto k=2 threshold=0.0 NxN=full"})
        self.assertEquals(len(lines), 7, lines)
        self.assertEquals(lines, ['0\tzeros',
                                  '2\tbit1\ttwo_bits\t0.500\tseveral\t0.143',
                                  '2\ttwo_bits\tbit1\t0.500\tseveral\t0.286',
                                  '2\tseveral\ttwo_bits\t0.286\tDEADdead\t0.261',
                                  '2\tdeadbeef\tDeaf Beef\t0.960\tDEADdead\t0.840',
                                  '2\tDEADdead\tdeadbeef\t0.840\tDeaf Beef\t0.808',
                                  '2\tDeaf Beef\tdeadbeef\t0.960\tDEADdead\t0.808'])
        
    def test_count_with_low_batch_size(self):
        header, lines = count_run_split("--NxN --count --batch-size 1", 7, SIMPLE_FPS)
        self.assertIn("simple.fps", header.pop("#targets"))
        self.assertEquals(header,
                          {"#num_bits": "32",
                          "#software": SOFTWARE,
                          "#type": "Count threshold=0.7 NxN=full"})
        self.assertEquals(len(lines), 7, lines)
        self.assertEquals(lines, ['0\tzeros',
                                  '0\tbit1',
                                  '0\ttwo_bits',
                                  '0\tseveral',
                                  '2\tdeadbeef',
                                  '2\tDEADdead',
                                  '2\tDeaf Beef',
            ])

        
        
        
class TestCompatibility(unittest2.TestCase):
    def test_incompatible_fingerprint(self):
        errmsg = run_exit(["--hex-query", "dead"], SIMPLE_FPS)
        self.assertIn("error: query fingerprint contains 2 bytes but", errmsg)
        self.assertIn("simple.fps", errmsg)
        self.assertIn("has 4 byte fingerprints", errmsg)

    def test_targets_is_not_an_fps_file(self):
        errmsg = run_exit(["--queries", SIMPLE_FPS])
        self.assertIn("Cannot open targets file:", errmsg)
        self.assertIn("Unable to determine fingerprint format type from", errmsg)
        self.assertIn("pubchem.sdf", errmsg)

    def test_targets_does_not_exist(self):
        errmsg = run_exit(["--queries", SIMPLE_FPS], "/this/file/does_not_exist_t")
        self.assertIn("Cannot open targets file:", errmsg)
        self.assertIn("No such file or directory", errmsg) # Mac specific?
        self.assertIn("does_not_exist_t", errmsg)

    def test_queries_is_not_an_fps_file(self):
        errmsg = run_exit(["--queries", support.PUBCHEM_SDF], SIMPLE_FPS)
        self.assertIn("Cannot open queries file:", errmsg)
        self.assertIn("Unable to determine fingerprint format type from ", errmsg)
        self.assertIn("pubchem.sdf", errmsg)

    def test_queries_does_not_exist(self):
        errmsg = run_exit(["--queries", "/this/file/does_not_exist_q"], SIMPLE_FPS)
        self.assertIn("Cannot open queries file:", errmsg)
        self.assertIn("No such file or directory", errmsg) # Mac specific?
        self.assertIn("does_not_exist_q", errmsg)

class TestCommandlineErrors(unittest2.TestCase):
    def test_mix_count_and_knearest(self):
        errmsg = count_run_exit("--count --hex-query beefcafe --k-nearest 4", SIMPLE_FPS)
        self.assertIn("--count search does not support --k-nearest", errmsg)
        
    def test_negative_k(self):
        errmsg = run_exit("--hex-query beefcafe -k -1", SIMPLE_FPS)
        self.assertIn("--k-nearest must be non-negative or 'all'", errmsg)

    def test_negative_threshold(self):
        errmsg = run_exit("--hex-query beefcafe --threshold -0.1", SIMPLE_FPS)
        self.assertIn("--threshold must be between 0.0 and 1.0, inclusive", errmsg)
        errmsg = run_exit("--hex-query beefcafe --threshold -1.0", SIMPLE_FPS)
        self.assertIn("--threshold must be between 0.0 and 1.0, inclusive", errmsg)

    def test_too_large_threshold(self):
        errmsg = run_exit("--hex-query beefcafe --threshold 1.1", SIMPLE_FPS)
        self.assertIn("--threshold must be between 0.0 and 1.0, inclusive", errmsg)

    def test_non_positive_batch_size(self):
        errmsg = run_exit("--hex-query beefcafe --batch-size 0", SIMPLE_FPS)
        self.assertIn("--batch-size must be positive", errmsg)
        errmsg = run_exit("--hex-query beefcafe --batch-size -1", SIMPLE_FPS)
        self.assertIn("--batch-size must be positive", errmsg)

    def test_NxN_with_scan(self):
        errmsg = run_exit("--NxN --scan", SIMPLE_FPS)
        self.assertIn("Cannot specify --scan with an --NxN search", errmsg)

    def test_NxN_with_hex_query(self):
        errmsg = run_exit("--NxN --hex-query feedfeed", SIMPLE_FPS)
        self.assertIn("Cannot specify --hex-query with an --NxN search", errmsg)
        
    def test_NxN_with_queries(self):
        errmsg = run_exit("--NxN --queries ignored", SIMPLE_FPS)
        self.assertIn("Cannot specify --queries with an --NxN search", errmsg)

    def test_scan_with_memory(self):
        errmsg = run_exit("--scan --memory", SIMPLE_FPS)
        self.assertIn("Cannot specify both --scan and --memory", errmsg)
        
    def test_hex_query_with_queries(self):
        errmsg = run_exit("--hex-query faceb00c --queries not_important", SIMPLE_FPS)
        self.assertIn("Cannot specify both --hex-query and --queries", errmsg)

    def test_hex_query_with_bad_character(self):
        errmsg = run_exit("--hex-query faceb00k", SIMPLE_FPS)
        self.assertIn("--hex-query is not a hex string: Non-hexadecimal digit found", errmsg)
        
    def test_hex_query_with_bad_length(self):
        errmsg = run_exit("--hex-query deadbeef2", SIMPLE_FPS)
        self.assertIn("--hex-query is not a hex string: Odd-length string", errmsg)

    def test_query_id_with_bad_character(self):
        for (bad_id, name) in (("A\tB", "tab"), ("C\nD", "newline"),
                               ("E\rF", "control-return"), ("G\0H", "NUL")):
            errmsg = run_exit(["--hex-query", "abcd1234", "--query-id", bad_id], SIMPLE_FPS)
            self.assertIn("--query-id must not contain the %s character" % name, errmsg)
        
    
if __name__ == "__main__":
    unittest2.main()
