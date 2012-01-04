# This tests that the tox environment is set up the way it should be
# set up.

import sys
import os
import unittest2

import support

envstr = os.environ.get("TOX_CHEMFP_TEST", None)
if not envstr:
    versions = []
else:
    versions = envstr.split(",")

def check_py25():
    version = sys.version_info[:2]
    assert version == (2,5), version
    assert not support.can_skip("py25")

def check_py26():
    version = sys.version_info[:2]
    assert version == (2,6), version
    assert not support.can_skip("py26")

def check_py27():
    version = sys.version_info[:2]
    assert version == (2,7), version
    assert not support.can_skip("py27")

def check_x32():
    assert sys.maxint == 2147483647
    assert not support.can_skip("x32")

def check_x64():
    assert sys.maxint == 9223372036854775807
    assert not support.can_skip("x64")
    
def check_oe174():
    from openeye.oechem import OEChemGetRelease
    version = OEChemGetRelease()
    assert version == "1.7.4", version
    assert not support.can_skip("oe")
    assert not support.can_skip("oe174")

def check_ob223():
    import openbabel
    from chemfp import openbabel
    assert openbabel._ob_version == "2.2.3", openbabel._ob_version
    assert not support.can_skip("ob")
    assert not support.can_skip("ob223")

def check_ob230():
    import openbabel
    from chemfp import openbabel
    assert openbabel._ob_version == "2.3.0", openbabel._ob_version
    assert not support.can_skip("ob")
    assert not support.can_skip("ob230")

def check_ob23svn1():
    import openbabel
    from chemfp import openbabel
    assert openbabel._ob_version == "2.3.90", openbabel._ob_version
    assert not support.can_skip("ob")
    assert not support.can_skip("ob23svn1")

def check_rd201106():
    from rdkit.rdBase import rdkitVersion
    assert rdkitVersion[:7] == "2011.06", rdkitVersion

def check_rd201103():
    import rdkit.rdBase
    from rdkit.rdBase import rdkitVersion
    assert rdkitVersion[:7] == "2011.03", rdkitVersion

def check_rd201012():
    import rdkit.rdBase
    from rdkit.rdBase import rdkitVersion
    assert rdkitVersion[:7] == "2010.12", rdkitVersion

def check_rd201112_svn():
    import rdkit.rdBase
    from rdkit.rdBase import rdkitVersion
    assert rdkitVersion[:7] == "2011.12", rdkitVersion

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
        _check("oe174 ob223 ob230 ob23svn1 rd201106 rd201103 rd201012 rd201112_svn")
        
    def test_version(self):
        for name in versions:
            func = globals()["check_" + name]
            func()

TestToxVersion = unittest2.skipUnless(envstr, "Not building under the tox environment")(
    TestToxVersion)
