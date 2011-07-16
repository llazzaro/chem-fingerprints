# This tests that the tox environment is set up the way it should be
# set up.

import sys
import os
import unittest2

versions = os.environ.get("TOX_CHEMFP_TEST", None)
if not versions:
    versions = []
else:
    versions = versions.split(",")

def check_py25():
    version = sys.version_info[:2]
    assert version == (2,5), version

def check_py26():
    version = sys.version_info[:2]
    assert version == (2,6), version

def check_py27():
    version = sys.version_info[:2]
    assert version == (2,7), version

def check_x32():
    assert sys.maxint == 2147483647

def check_x64():
    assert sys.maxint == 9223372036854775807
    
def check_oe161():
    from openeye.oechem import OEChemGetRelease
    version = OEChemGetRelease()
    assert version == "1.6.1", version

def check_oe174():
    from openeye.oechem import OEChemGetRelease
    version = OEChemGetRelease()
    assert version == "1.7.4", version

def check_ob223():
    pass


def _check(required):
    req = required.split()
    for name in versions:
        if name in req:
            return
    raise AssertionError("Missing one of %r: %r" % (required, versions))

class TestToxVersion(unittest2.TestCase):
    def test_enough_specifications(self):
        _check("x32 x64")
        _check("py25 py26 py27 py32")
        _check("oe174 ob223 ob230")
        
    def test_version(self):
        for name in versions:
            func = globals()["check_" + name]
            func()
