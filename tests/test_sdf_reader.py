import sys
import unittest
from cStringIO import StringIO as SIO

from chemfp import sdf_reader

# At some point make this independent of where the tests are started
TRYPTOPHAN_SDF = "tryptophan.sdf"
PUBCHEM_SDF = "pubchem.sdf"
PUBCHEM_SDF_GZ = "pubchem.sdf.gz"

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

class FakeDecompressor(object):
    def open_filename_universal(self, filename):
        assert filename == "spam.blah"
        return open(PUBCHEM_SDF, "rU")

class TestReadRecords(unittest.TestCase):
    def test_reads_the_only_record(self):
        n = sum(1 for x in sdf_reader.open_sdfile(TRYPTOPHAN_SDF))
        self.assertEquals(n, 1)
    def test_reads_all_records(self):
        n = sum(1 for x in sdf_reader.open_sdfile(PUBCHEM_SDF))
        self.assertEquals(n, 19)
    def test_reads_all_compressed_records(self):
        n = sum(1 for x in sdf_reader.open_sdfile(PUBCHEM_SDF_GZ))
        self.assertEquals(n, 19)

    def test_reads_from_stdin(self):
        old_stdin = sys.stdin
        sys.stdin = open(PUBCHEM_SDF, "rb")
        try:
            n = sum(1 for x in sdf_reader.open_sdfile())
        finally:
            sys.stdin = sys.stdin
        self.assertEquals(n, 19)

    def test_reads_from_gzip_stdin(self):
        old_stdin = sys.stdin
        sys.stdin = open(PUBCHEM_SDF_GZ, "rb")
        try:
            n = sum(1 for x in sdf_reader.open_sdfile(None, "gzip"))
        finally:
            sys.stdin = sys.stdin
        self.assertEquals(n, 19)

    def test_reads_from_fileobj(self):
        f = open(PUBCHEM_SDF, "rU")
        n = sum(1 for x in sdf_reader.open_sdfile(f))
        self.assertEquals(n, 19)
        
    def test_reads_from_gzip_fileobj(self):
        f = open(PUBCHEM_SDF_GZ, "rb")
        n = sum(1 for x in sdf_reader.open_sdfile(f, "gzip"))
        self.assertEquals(n, 19)

    def test_handles_alternate_decompressor(self):
        n = sum(1 for x in sdf_reader.open_sdfile("spam.blah", FakeDecompressor()))
        self.assertEquals(n, 19)

    def test_handles_loc(self):
        loc = sdf_reader.FileLocation()
        results = []
        for x in sdf_reader.open_sdfile(PUBCHEM_SDF_GZ, loc=loc):
            self.assertEquals(loc.filename, PUBCHEM_SDF_GZ)
            results.append(dict(title=loc.title,
                                lineno=loc.lineno))
        self.assertEquals(results, expected_locs)

tryptophan = open(TRYPTOPHAN_SDF).read()

class ReadReturnsSmallAmounts(object):
    def __init__(self):
        self.f = open(PUBCHEM_SDF, "rU")
    def read(self, n):
        return self.f.read(1)

class ReadReturnsOneRecord(object):
    def __init__(self):
        self.f = sdf_reader.open_sdfile(PUBCHEM_SDF)
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
        self.f = sdf_reader.open_sdfile(PUBCHEM_SDF)
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

class TestBoundaryConditions(unittest.TestCase):
    def test_missing_terminal_newline(self):
        f = SIO(tryptophan.rstrip("\n"))
        n = sum(1 for x in sdf_reader.open_sdfile(f))
        self.assertEquals(n, 1)
    def test_small_amounts(self):
        # the O(n**2) behavior really hits hard here - this takes a full second to work
        n = sum(1 for x in sdf_reader.open_sdfile(ReadReturnsSmallAmounts()))
        self.assertEquals(n, 19)
    def test_exact_record_boundary_reads(self):
        loc = sdf_reader.FileLocation()
        titles = [loc.title for x in sdf_reader.open_sdfile(ReadReturnsOneRecord(), loc=loc)]
        self.assertEquals(titles, expected_identifiers)
    def test_two_record_boundary_reads(self):
        loc = sdf_reader.FileLocation()
        titles = [loc.title for x in sdf_reader.open_sdfile(ReadReturnsTwoRecords(),
                                                            loc=loc)]
        self.assertEquals(titles, expected_identifiers)

class TestReadErrors(unittest.TestCase):
    def test_record_too_large(self):
        f = SIO( (tryptophan * ((200000 // len(tryptophan)) + 1)).replace("$$$$", "1234"))
        try:
            for x in sdf_reader.open_sdfile(f):
                raise AssertionError("should not be able to read the first record")
        except TypeError, err:
            self.assertEquals("too large" in str(err), True)
            self.assertEquals("at line 1" in str(err), True)

    def test_has_extra_data(self):
        f = SIO(tryptophan + tryptophan + "blah")
        try:
            for i, x in enumerate(sdf_reader.open_sdfile(f)):
                if i > 1:
                    raise AssertionError("bad record count")
        except TypeError, err:
            self.assertEquals("unexpected content" in str(err), True)
            expected_lineno = (tryptophan.count("\n")*2) + 1
            expected_lineno_msg = "at line %d" % expected_lineno
            self.assertEquals(expected_lineno_msg in str(err), True, str(err))

    def test_bad_format(self):
        f = SIO(tryptophan + tryptophan.replace("V2000", "V4000"))
        try:
            for i, x in enumerate(sdf_reader.open_sdfile(f)):
                if i > 0:
                    raise AssertionError("bad record count")
        except TypeError, err:
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
        titles = [loc.lineno for rec in sdf_reader.open_sdfile(
            f, loc=loc, reader_error=my_error_handler)]
        self.assertEquals(titles, [1, 137])
        self.assertEquals(my_error_handler.errors, [("incorrectly formatted record",
                                                     {"filename": None,
                                                      "lineno": 70,
                                                      "title": "tryptophan.pdb"})])
        
            
class TestReadTitleTagAndFPTag(unittest.TestCase):
    pass
class TestReadTitleAndFPTag(unittest.TestCase):
    pass



if __name__ == "__main__":
    unittest.main()
