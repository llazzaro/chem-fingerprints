import sys
import os
from cStringIO import StringIO
import tempfile

# Ignore the close. io.write_fps1_output() auto-closes its output.
class SIO(object):
    def __init__(self):
        self.sio = StringIO()
    def write(self, s):
        return self.sio.write(s)
    def writelines(self, lines):
        return self.sio.writelines(lines)
    def close(self):
        # Ignore this
        pass
    def getvalue(self):
        return self.sio.getvalue()

# Given a filename in the "tests/" directory, return its full path

_dirname = os.path.dirname(__file__)
def fullpath(name):
    path = os.path.join(_dirname, name)
    assert os.path.exists(path), path
    return path

PUBCHEM_SDF = fullpath("pubchem.sdf")
PUBCHEM_SDF_GZ = fullpath("pubchem.sdf.gz")
PUBCHEM_ANOTHER_EXT = fullpath("pubchem.should_be_sdf_but_is_not")
EXAMPLES_INCHI = fullpath("examples.inchi")

MISSING_TITLE = fullpath("missing_title.sdf")

real_stdin = sys.stdin
real_stdout = sys.stdout
real_stderr = sys.stderr

class Runner(object):
    def __init__(self, main):
        self.main = main

    def pre_run(self):
        pass
    def post_run(self):
        pass

    def run(self, cmdline, source=PUBCHEM_SDF):
        if isinstance(cmdline, basestring):
            args = cmdline.split()
        else:
            args = cmdline
            assert isinstance(args, list) or isinstance(args, tuple)
        if source is not None:
            args = args + [source]
        self.pre_run()

        try:
            sys.stdout = stdout = SIO()
            self.main(args)
        finally:
            sys.stdout = real_stdout

        self.post_run()

        result = stdout.getvalue().splitlines()
        if result:
            self.verify_result(result)
        return result

    def verify_result(self, result):
        assert result[0] == "#FPS1", result[0]
        # TODO: .. verify more more line format ...

    def run_stdin(self, cmdline):
        raise NotImplementedError("Implement in the derived class")

    def run_fps(self, cmdline, expect_length=None, source=PUBCHEM_SDF):
        result = self.run(cmdline, source)
        while result[0].startswith("#"):
            del result[0]
        if expect_length is not None:
            assert len(result) == expect_length, (len(result), expect_length)
        return result

    def run_split(self, cmdline, expect_length=None, source=PUBCHEM_SDF):
        "split into dict of headers and list of values"
        result = self.run(cmdline, source)
        headers = {}
        fps = []
        result_iter = iter(result)
        # I know the first line is correct (it was tested in verify_result)
        # Plus, this lets the SimsearchRunner use run_split
        result_iter.next()
        for line in result_iter:
            if line.startswith("#"):
                k, v = line.split("=", 1)
                assert k not in headers, k
                headers[k] = v
                continue
            fps.append(line)
            break
        fps.extend(result_iter)
        if expect_length is not None:
            assert len(fps) == expect_length, (len(fps), expect_length)
        return headers, fps
            

    def run_exit(self, cmdline, source=PUBCHEM_SDF):
        sys.stderr = stderr = SIO()
        try:
            try:
                self.run(cmdline, source)
            except SystemExit:
                pass
            else:
                raise AssertionError("should have exited: %r" % (cmdline,))
        finally:
            sys.stderr = real_stderr
        return stderr.getvalue()

    def run_split_capture(self, cmdline, expect_length=None, source=PUBCHEM_SDF):
        sys.stderr = stderr = SIO()
        try:
            try:
                headers, fps = self.run_split(cmdline, expect_length, source)
            except SystemExit:
                raise AssertionError("unexpected SystemExit")
        finally:
            sys.stderr = real_stderr
        return headers, fps, stderr.getvalue()
        

####

def can_skip(name):
    s = os.environ.get("TOX_CHEMFP_TEST", "")
    return not (s.startswith(name) or (","+name) in s)

#### fingerprint encoding

def set_bit(n):
    assert n <= 16
    bytes = [0, 0, 0]
    bytes[n//8] = 1<<(n%8)
    return "%02x%02x%02x" % tuple(bytes)

class TestIdAndErrors(object):
    #
    # One of the records doesn't have an XLOGP field
    #
    def test_missing_id_tag(self):
        errmsg = self._runner.run_exit("--id-tag PUBCHEM_CACTVS_XLOGP")
        self.assertIn("ERROR: Missing id tag 'PUBCHEM_CACTVS_XLOGP' for record #7 ", errmsg)
        self.assertIn("pubchem.sdf", errmsg)

    # Should be the same as the previous code.
    def test_missing_id_strict(self):
        errmsg = self._runner.run_exit("--id-tag PUBCHEM_CACTVS_XLOGP --errors strict")
        self.assertIn("ERROR: Missing id tag 'PUBCHEM_CACTVS_XLOGP' for record #7 ", errmsg)
        self.assertIn("pubchem.sdf", errmsg)
    

    def test_missing_id_tag_report(self):
        headers, fps, errmsg = self._runner.run_split_capture("--id-tag PUBCHEM_CACTVS_XLOGP --errors report", 18)
        self.assertIn("ERROR: Missing title for record #1", errmsg)
        self.assertIn("missing_title.sdf", errmsg)
        self.assertEquals(fps[-1], "")

    def test_missing_id_tag_ignore(self):
        headers, fps, errmsg = self._runner.run_split_capture("--id-tag PUBCHEM_CACTVS_XLOGP --errors ignore", 18)
        self.assertNotIn("ERROR: Missing title for record #1", errmsg)
        self.assertNotIn("missing_title.sdf", errmsg)
        ids = [fp.split("\t")[1] for fp in fps]
        self.assertEquals(ids, ['2.8', '1.9', '1', '3.3', '1.5', '2.6', '-0.9', '2', '2.1', 
                                '2.9', '1.7', '-1.5', '0.4', '0.6', '0.4', '0.4', '2', '2.5'])


    #
    # Various ways of having a strange title
    #

    def test_missing_title(self):
        errmsg = self._runner.run_exit("", MISSING_TITLE)
        self.assertIn("ERROR: Missing title for record #1", errmsg)

    def test_missing_title_strict(self):
        errmsg = self._runner.run_exit("--errors strict", MISSING_TITLE)
        self.assertIn("ERROR: Missing title for record #1", errmsg)

    def test_missing_title_report(self):
        headers, fps, errmsg = self._runner.run_split_capture("--errors report", 1, MISSING_TITLE)
        self.assertIn("ERROR: Missing title for record #1", errmsg)
        self.assertNotIn("ERROR: Missing title for record #2", errmsg)
        self.assertIn("ERROR: Missing title for record #3", errmsg)
        self.assertEquals(len(fps), 1)
        self.assertEquals(fps[0].split("\t")[1], "Good")

    def test_missing_title_ignore(self):
        headers, fps, errmsg = self._runner.run_split_capture("--errors ignore", 1, MISSING_TITLE)
        self.assertNotIn("ERROR: Missing title for record #1", errmsg)
        self.assertNotIn("ERROR: Missing title for record #2", errmsg)
        self.assertNotIn("ERROR: Missing title for record #3", errmsg)
        self.assertEquals(len(fps), 1)
        self.assertEquals(fps[0].split("\t")[1], "Good")

    #
    # Various ways of handling a missing id in a tag
    #

    def test_missing_id_tag(self):
        errmsg = self._runner.run_exit("--id-tag Blank", MISSING_TITLE)
        self.assertIn("ERROR: Empty id tag 'Blank' for record #1", errmsg)

    def test_missing_id_tag_strict(self):
        errmsg = self._runner.run_exit("--id-tag Blank --errors strict", MISSING_TITLE)
        self.assertIn("ERROR: Empty id tag 'Blank' for record #1", errmsg)
        self.assertIn("missing_title.sdf", errmsg)

    def test_missing_id_tag_report(self):
        headers, fps, errmsg = self._runner.run_split_capture("--id-tag Blank --errors report", 1, MISSING_TITLE)
        self.assertIn("ERROR: Empty id tag 'Blank' for record #1", errmsg)
        self.assertIn("ERROR: Empty id tag 'Blank' for record #2", errmsg)
        self.assertNotIn("ERROR: Empty id tag 'Blank' for record #3", errmsg)
        self.assertEquals(fps[0].split("\t")[1], "This is not Blank")

    def test_missing_id_tag_ignore(self):
        headers, fps, errmsg = self._runner.run_split_capture("--id-tag Blank --errors ignore", 1, MISSING_TITLE)
        self.assertNotIn("ERROR: Empty id tag 'Blank' for record #1", errmsg)
        self.assertNotIn("ERROR: Empty id tag 'Blank' for record #2", errmsg)
        self.assertNotIn("ERROR: Empty id tag 'Blank' for record #3", errmsg)
        self.assertEquals(fps[0].split("\t")[1], "This is not Blank")

    #
    # Various ways of handling a tab characters in an id tag
    #

    def test_tab_id_tag(self):
        errmsg = self._runner.run_exit("--id-tag Tab", MISSING_TITLE)
        self.assertIn("ERROR: Empty id tag 'Tab' for record #2", errmsg)

    def test_tab_id_tag_strict(self):
        errmsg = self._runner.run_exit("--id-tag Tab --errors strict", MISSING_TITLE)
        self.assertIn("ERROR: Empty id tag 'Tab' for record #2", errmsg)
        self.assertIn("missing_title.sdf", errmsg)

    def test_tab_id_tag_report(self):
        headers, fps, errmsg = self._runner.run_split_capture("--id-tag Tab --errors report", 2, MISSING_TITLE)
        self.assertIn("ERROR: Empty id tag 'Tab' for record #2", errmsg)
        self.assertEquals(fps[0].split("\t")[1], "Leading tab")
        self.assertEquals(fps[1].split("\t")[1], "This does not")

    def test_tab_id_tag_ignore(self):
        headers, fps, errmsg = self._runner.run_split_capture("--id-tag Tab --errors ignore", 2, MISSING_TITLE)
        self.assertNotIn("ERROR: Empty id tag 'Tab'", errmsg)
        self.assertEquals(fps[0].split("\t")[1], "Leading tab")
        self.assertEquals(fps[1].split("\t")[1], "This does not")


    def test_contains_tab_id_tag(self):
        headers, fps = self._runner.run_split("--id-tag ContainsTab", 3, MISSING_TITLE)
        ids = [fp.split("\t")[1] for fp in fps]
        self.assertEquals(ids, ["ThreeTabs", "tabseparated", "twotabs"])

    def test_contains_tab_id_tag_strict(self):
        headers, fps = self._runner.run_split("--id-tag ContainsTab --errors strict", 3, MISSING_TITLE)
        ids = [fp.split("\t")[1] for fp in fps]
        self.assertEquals(ids, ["ThreeTabs", "tabseparated", "twotabs"])

    def test_contains_tab_id_tag_report(self):
        headers, fps, errmsg = self._runner.run_split_capture("--id-tag ContainsTab --errors report", 3, MISSING_TITLE)
        self.assertNotIn("ContainsTab", errmsg)
        self.assertNotIn("ERROR", errmsg)
        ids = [fp.split("\t")[1] for fp in fps]
        self.assertEquals(ids, ["ThreeTabs", "tabseparated", "twotabs"])

    def test_contains_tab_id_tag_ignore(self):
        headers, fps, errmsg = self._runner.run_split_capture("--id-tag ContainsTab --errors ignore", 3, MISSING_TITLE)
        self.assertNotIn("ERROR: Empty id tag 'ContainsTab'", errmsg)
        ids = [fp.split("\t")[1] for fp in fps]
        self.assertEquals(ids, ["ThreeTabs", "tabseparated", "twotabs"])

    #
    # Handling bad files
    #

    def test_handles_missing_filename(self):
        errmsg = self._runner.run_exit("this_file_does_not_exist.sdf", PUBCHEM_SDF)
        self.assertIn("Structure file '", errmsg)
        self.assertIn("this_file_does_not_exist.sdf", errmsg)
        self.assertIn("' does not exist", errmsg)
        self.assertNotIn("pubchem", errmsg)

    def test_handles_missing_filename_at_end(self):
        errmsg = self._runner.run_exit([PUBCHEM_SDF, "this_file_does_not_exist.sdf"])
        self.assertIn("Structure file '", errmsg)
        self.assertIn("this_file_does_not_exist.sdf", errmsg)
        self.assertIn("' does not exist", errmsg)
        self.assertNotIn("pubchem", errmsg)

    def test_unreadable_file(self):
        tf = tempfile.NamedTemporaryFile(suffix="unreadable.sdf")
        try:
            os.chmod(tf.name, 0222)
            errmsg = self._runner.run_exit([PUBCHEM_SDF, tf.name])
            self.assertIn("Problem reading structure fingerprints", errmsg)
            self.assertIn("unreadable.sdf", errmsg)
            self.assertNotIn("pubchem", errmsg)
        finally:
            tf.close()
