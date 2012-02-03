from __future__ import with_statement

import unittest2
from chemfp import types

try:
    from chemfp import openeye
except ImportError:
    skip_openeye = True
else:
    skip_openeye = False

try:
    from chemfp import rdkit
except ImportError:
    skip_rdkit = True
else:
    skip_rdkit = False


# These tests are to improve overall test coverage

class ParseType(unittest2.TestCase):
    def test_no_type_string(self):
        with self.assertRaisesRegexp(ValueError, r"Empty fingerprint type \(''\)"):
            types.parse_type("")
        with self.assertRaisesRegexp(ValueError, r"Empty fingerprint type \(' '\)"):
            types.parse_type(" ")

    def test_unknown_fingerprint_family(self):
        with self.assertRaisesRegexp(ValueError, "Unknown fingerprint family 'Blah-Blah/1'"):
            types.parse_type("Blah-Blah/1")

class TestOEFingerprintTypes(unittest2.TestCase):
    def test_missing_equals(self):
        with self.assertRaisesRegexp(ValueError, "Term 'atype' in type 'OpenEye-Path atype' must have one and only one '='"):
            types.parse_type("OpenEye-Path atype")

    def test_extra_equals(self):
        with self.assertRaisesRegexp(ValueError, "Term 'atype=DefaultAtom=DefaultAtom' in type 'OpenEye-Path atype=DefaultAtom=DefaultAtom' must have one and only one '='"):
            types.parse_type("OpenEye-Path atype=DefaultAtom=DefaultAtom")

    def test_duplicate_names(self):
        with self.assertRaisesRegexp(ValueError, "Duplicate name 'atype' in type 'OpenEye-Path atype=DefaultAtom atype=InRing'"):
            types.parse_type("OpenEye-Path atype=DefaultAtom atype=InRing")

    def test_unknown_name(self):
        with self.assertRaisesRegexp(ValueError, "Unknown name 'atomtype' in type 'OpenEye-Path atomtype=DefaultAtom'"):
            types.parse_type("OpenEye-Path atomtype=DefaultAtom")

    def test_negative_minbonds(self):
        with self.assertRaisesRegexp(ValueError, "Unable to parse minbonds value '-1' in type 'OpenEye-Path minbonds=-1'"):
            types.parse_type("OpenEye-Path minbonds=-1")

    def test_make_fingerprinter_from_type_with_empty_string(self):
        family = types.get_fingerprint_family("OpenEye-Path")
        with self.assertRaisesRegexp(ValueError, r"Empty fingerprint type \(''\)"):
            family.make_fingerprinter_from_type("")

    def test_compare_types_for_equality(self):
        t1 = types.parse_type("OpenEye-Path")
        self.assertEquals(t1, t1)
        t2 = types.parse_type("OpenEye-Path")
        self.assertEquals(t1, t2)

    def test_compare_types_for_inequality(self):
        t1 = types.parse_type("OpenEye-Path")
        t2 = types.parse_type("OpenEye-Path atype=InRing")
        self.assertNotEquals(t1, t2)
        
TestOEFingerprintTypes = unittest2.skipIf(skip_openeye, "OEChem not installed")(TestOEFingerprintTypes)

class TestRDKitFingerprintTypes(unittest2.TestCase):
    def test_passing_0_to_positive_int(self):
        with self.assertRaisesRegexp(ValueError, "Unable to parse fpSize value '0' in type 'RDKit-Fingerprint fpSize=0'"):
            types.parse_type("RDKit-Fingerprint fpSize=0")
            
    def test_passing_negative_to_nonnegative_int(self):
        with self.assertRaisesRegexp(ValueError, r"Unable to parse fpSize value '\+2' in type 'RDKit-Fingerprint fpSize=\+2'"):
            types.parse_type("RDKit-Fingerprint fpSize=+2")

    def test_passing_something_other_than_0_or_1(self):
        with self.assertRaisesRegexp(ValueError, r"Unable to parse useHs value '2' in type 'RDKit-Fingerprint useHs=2'"):
            types.parse_type("RDKit-Fingerprint useHs=2")

TestRDKitFingerprintTypes = unittest2.skipIf(skip_rdkit, "OEChem not installed")(TestRDKitFingerprintTypes)

if __name__ == "__main__":
    unittest2.main()
