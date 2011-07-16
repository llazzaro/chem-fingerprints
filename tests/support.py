import sys
import os
from cStringIO import StringIO

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
        assert result[0] == "#FPS1"
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
                assert k not in headers
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
