from __future__ import with_statement
import sys
import unittest2
from cStringIO import StringIO as SIO

import support

# Make sure the module correctly implements __all__
before = after = None
before = set(globals())
from chemfp.sdf_reader import *
after = set(globals())
assert len(after - before) == 4, ("wrong import * count", after-before)

# Needed for access to the experimental FileLocation
from chemfp import sdf_reader

TRYPTOPHAN_SDF = support.fullpath("tryptophan.sdf")
PUBCHEM_SDF = support.fullpath("pubchem.sdf")
PUBCHEM_SDF_GZ = support.fullpath("pubchem.sdf.gz")
STRANGE_SDF = support.fullpath("strange.sdf")

expected_identifiers = ["9425004", "9425009", "9425012", "9425015",
                        "9425018", "9425021", "9425030", "9425031",
                        "9425032", "9425033", "9425034", "9425035",
                        "9425036", "9425037", "9425040", "9425041",
                        "9425042", "9425045", "9425046"]
expected_linenos = [1, 191, 401, 592, 817, 1027, 1236, 1457, 1678,
                    1872, 2086, 2288, 2493, 2700, 2894, 3103, 3305,
                    3507, 3722]
assert len(expected_identifiers) == len(expected_linenos)
expected_locs = [dict(title=title, lineno=lineno) for (title, lineno) in
                     zip(expected_identifiers, expected_linenos)]

class TestReadRecords(unittest2.TestCase):
    def test_reads_the_only_record(self):
        n = sum(1 for x in open_sdf(TRYPTOPHAN_SDF))
        self.assertEquals(n, 1)
    def test_reads_all_records(self):
        n = sum(1 for x in open_sdf(PUBCHEM_SDF))
        self.assertEquals(n, 19)
    def test_reads_all_compressed_records(self):
        n = sum(1 for x in open_sdf(PUBCHEM_SDF_GZ))
        self.assertEquals(n, 19)

    def test_reads_from_stdin(self):
        old_stdin = sys.stdin
        sys.stdin = open(PUBCHEM_SDF, "rb")
        try:
            n = sum(1 for x in open_sdf())
        finally:
            sys.stdin = sys.stdin
        self.assertEquals(n, 19)

    def test_reads_from_gzip_stdin(self):
        old_stdin = sys.stdin
        sys.stdin = open(PUBCHEM_SDF_GZ, "rb")
        try:
            n = sum(1 for x in open_sdf(None, "gzip"))
        finally:
            sys.stdin = sys.stdin
        self.assertEquals(n, 19)

    def test_reads_from_fileobj(self):
        f = open(PUBCHEM_SDF, "rU")
        n = sum(1 for x in open_sdf(f))
        self.assertEquals(n, 19)
        
    def test_reads_from_uncompressed_fileobj(self):
        f = open(PUBCHEM_SDF, "rU")
        n = sum(1 for x in open_sdf(f, "none"))
        self.assertEquals(n, 19)

    def test_reads_from_gzip_fileobj(self):
        f = open(PUBCHEM_SDF_GZ, "rb")
        n = sum(1 for x in open_sdf(f, "gzip"))
        self.assertEquals(n, 19)

    def test_handles_loc(self):
        loc = sdf_reader.FileLocation()
        results = []
        for x in open_sdf(PUBCHEM_SDF_GZ, loc=loc):
            if sys.version_info[:3] > (2, 5, 4):
                # Earlier versions of the gzip library didn't
                # keep track of the .name attribute
                self.assertEquals(loc.name, PUBCHEM_SDF_GZ)
            else:
                self.assertEquals(getattr(loc, "name"), None)

            results.append(dict(title=loc.title,
                                lineno=loc.lineno))
        self.assertEquals(results, expected_locs)

    def test_when_using_wrong_compression(self):
        try:
            n = sum(1 for x in open_sdf(PUBCHEM_SDF, "gzip"))
            raise AssertionError("parsed a gzip'ed file?")
        except IOError:
            pass
            
        

tryptophan = open(TRYPTOPHAN_SDF).read()

class ReadReturnsSmallAmounts(object):
    def __init__(self):
        self.f = open(PUBCHEM_SDF, "rU")
    def read(self, n):
        return self.f.read(1)

class ReadReturnsOneRecord(object):
    def __init__(self):
        self.f = open_sdf(PUBCHEM_SDF)
    def read(self, n):
        try:
            rec = self.f.next()
        except StopIteration:
            return ""
        assert rec.endswith("\n$$$$\n")
        rec = rec[:-6] + ">spam\n" + ("X" * (n-len(rec)-6)) + "\n$$$$\n"
        assert len(rec) == n
        return rec

class ReadReturnsTwoRecords(object):
    def __init__(self):
        self.f = open_sdf(PUBCHEM_SDF)
    def read(self, n):
        try:
            rec = self.f.next()
        except StopIteration:
            return ""
        try:
            rec2 = self.f.next()
        except StopIteration:
            return rec
        rec = rec + rec2
        assert rec.endswith("\n$$$$\n")
        rec = rec[:-6] + ">spam\n" + ("X" * (n-len(rec)-6)) + "\n$$$$\n"
        assert len(rec) == n
        return rec

class TestBoundaryConditions(unittest2.TestCase):
    def test_missing_terminal_newline(self):
        f = SIO(tryptophan.rstrip("\n"))
        n = sum(1 for x in open_sdf(f))
        self.assertEquals(n, 1)
    def test_small_amounts(self):
        # the O(n**2) behavior really hits hard here - this takes a full second to work
        n = sum(1 for x in open_sdf(ReadReturnsSmallAmounts()))
        self.assertEquals(n, 19)
    def test_exact_record_boundary_reads(self):
        loc = sdf_reader.FileLocation()
        titles = [loc.title for x in open_sdf(ReadReturnsOneRecord(), loc=loc)]
        self.assertEquals(titles, expected_identifiers)
    def test_two_record_boundary_reads(self):
        loc = sdf_reader.FileLocation()
        titles = [loc.title for x in open_sdf(ReadReturnsTwoRecords(), loc=loc)]
        self.assertEquals(titles, expected_identifiers)

class TestReadErrors(unittest2.TestCase):
    def test_wrong_format(self):
        f = SIO("Spam\n")
        try:
            for x in open_sdf(f):
                raise AssertionError("Bad parse")
        except sdf_reader.SDFParseError, err:
            self.assertEquals("Could not find a valid SD record" in str(err), True)
            self.assertEquals("line 1" in str(err), True, str(err))
            
    def test_record_too_large(self):
        f = SIO( (tryptophan * ((200000 // len(tryptophan)) + 1)).replace("$$$$", "1234"))
        try:
            for x in open_sdf(f):
                raise AssertionError("should not be able to read the first record")
        except sdf_reader.SDFParseError, err:
            self.assertEquals("too large" in str(err), True, str(err))
            self.assertEquals("at line 1" in str(err), True, str(err))

    def test_has_extra_data(self):
        f = SIO(tryptophan + tryptophan + "blah")
        try:
            for i, x in enumerate(open_sdf(f)):
                if i > 1:
                    raise AssertionError("bad record count")
        except sdf_reader.SDFParseError, err:
            self.assertEquals("unexpected content" in str(err), True)
            expected_lineno = (tryptophan.count("\n")*2) + 1
            expected_lineno_msg = "at line %d" % expected_lineno
            self.assertEquals(expected_lineno_msg in str(err), True, str(err))

    def test_bad_format(self):
        f = SIO(tryptophan + tryptophan.replace("V2000", "V4000"))
        try:
            for i, x in enumerate(open_sdf(f)):
                if i > 0:
                    raise AssertionError("bad record count")
        except sdf_reader.SDFParseError, err:
            self.assertEquals("incorrectly formatted record" in str(err), True, str(err))
            self.assertEquals("at line 70" in str(err), True, str(err))

    def test_my_error_handler(self):
        class CaptureErrors(object):
            def __init__(self):
                self.errors = []
            def __call__(self, msg, loc):
                self.errors.append( (msg, loc.info()) )
        my_error_handler = CaptureErrors()
        loc = sdf_reader.FileLocation()
        f = SIO(tryptophan + tryptophan.replace("V2000", "V4000") + tryptophan)
        titles = [loc.lineno for rec in open_sdf(f, loc=loc,
                                                    errors=my_error_handler)]
        self.assertEquals(titles, [1, 137])
        self.assertEquals(my_error_handler.errors, [("incorrectly formatted record",
                                                     {"name": None,
                                                      "lineno": 70,
                                                      "title": "tryptophan.pdb"})])
        

expected_hbond_donors = ["2","2","2","2","2","2","4","4","2", "2",
                         "2","2","2","3","2","3","3","2","2"]
expected_complexity = ["491", "513", "419", "597", "545", "590", "660", "660",
                       "394", "544", "458", "589", "532", "506", "640", "557",
                       "557", "520", "437"]
expected_xlogp = ["2.8", "1.9", "1", "3.3", "1.5", "2.6", None, "-0.9",
                  "2", "2.1", "2.9", "1.7", "-1.5", "0.4", "0.6",
                  "0.4", "0.4", "2", "2.5"]

assert len(expected_complexity) == len(expected_hbond_donors) == len(expected_linenos)
assert len(expected_xlogp) == len(expected_linenos)

class TestIterTwoTags(unittest2.TestCase):
    def test_read_two_existing_tags(self):
        fields = list(iter_two_tags(open_sdf(PUBCHEM_SDF),
                            "PUBCHEM_CACTVS_HBOND_DONOR", "PUBCHEM_CACTVS_COMPLEXITY"))
        self.assertEquals(fields, zip(expected_hbond_donors, expected_complexity))

    def test_read_tag_missing_data_field1(self):
        fields = list(iter_two_tags(open_sdf(PUBCHEM_SDF),
                                    "PUBCHEM_CACTVS_XLOGP", "PUBCHEM_CACTVS_HBOND_DONOR"))
        self.assertEquals(fields, zip(expected_xlogp, expected_hbond_donors))

    def test_read_tag_missing_data_field2(self):
        fields = list(iter_two_tags(open_sdf(PUBCHEM_SDF),
                                    "PUBCHEM_CACTVS_HBOND_DONOR", "PUBCHEM_CACTVS_XLOGP"))
        self.assertEquals(fields, zip(expected_hbond_donors, expected_xlogp))
        
    def test_edge_conditions1(self):
        fields = list(iter_two_tags(open_sdf(STRANGE_SDF), "noblank", "twolines"))
        self.assertEquals(fields, [("This line is not followed by a blank line",
                                    "This contains two lines"), (None, None)])

    def test_edge_conditions2(self):
        fields = list(iter_two_tags(open_sdf(STRANGE_SDF),
                                    "duplicate", "embedded-tags"))
        self.assertEquals(fields, [("This is the first version.",
                    "I<junk> have<junk> tags  <junk>  in the data line <junk>"),
                                   (None, None)])

    def test_edge_conditions3(self):
        fields = list(iter_two_tags(open_sdf(STRANGE_SDF), "junk", "blank lines"))
        self.assertEquals(fields, [
      ("This line contains some of the strange junk that might exist on the tag line", ""),
            (None, None)])

    def test_edge_conditions4(self):
        fields = list(iter_two_tags(open_sdf(STRANGE_SDF), "nada", "fini"))
        self.assertEquals(fields, [(None, None), ("", "")])

    def test_bad_tags(self):
        for tag in ("<", ">", "\n", "\t", "1<2", "2<1", "blah\t"):
            self.assertRaises(TypeError, iter_two_tags([], tag, "fini"))
            self.assertRaises(TypeError, iter_two_tags([], "fini", tag))
    
class TestReadTitleAndTag(unittest2.TestCase):
    def test_read_existing_tag(self):
        fields = list(iter_title_and_tag(open_sdf(PUBCHEM_SDF),
                                         "PUBCHEM_CACTVS_HBOND_DONOR"))
        self.assertEquals(fields, zip(expected_identifiers, expected_hbond_donors))

    def test_missing_tag(self):
        fields = list(iter_title_and_tag(open_sdf(PUBCHEM_SDF), "PUBCHEM_CACTVS_XLOGP"))
                                         
        self.assertEquals(fields, zip(expected_identifiers, expected_xlogp))

    def test_bad_tags(self):
        for tag in ("<", ">", "\n", "\t", "1<2", "2<1", "blah\t"):
            self.assertRaises(TypeError, iter_title_and_tag([], tag))


if __name__ == "__main__":
    unittest2.main()
