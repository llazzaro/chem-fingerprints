import sys
from cStringIO import StringIO as SIO

PUBCHEM_SDF = "pubchem.sdf"
PUBCHEM_SDF_GZ = "pubchem.sdf.gz"
PUBCHEM_ANOTHER_EXT = "pubchem.should_be_sdf_but_is_not"

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
        args = cmdline.split()
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

    def run_exit(self, cmdline, source=PUBCHEM_SDF):
        sys.stderr = stderr = SIO()
        try:
            try:
                self.run(cmdline, source)
            except SystemExit:
                pass
            else:
                raise AssertionError("should have existed: %r" % (cmdline,))
        finally:
            sys.stderr = real_stderr
        return stderr.getvalue()
    

#### fingerprint encoding

def set_bit(n):
    assert n <= 16
    bytes = [0, 0, 0]
    bytes[n//8] = 1<<(n%8)
    return "%02x%02x%02x" % tuple(bytes)
