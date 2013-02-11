from __future__ import with_statement
import sys
import unittest2
import support

try:
    from openeye.oechem import OEGraphMol, OEParseSmiles, OEChemIsLicensed
    if not OEChemIsLicensed():
        raise ImportError
    skip_oechem = False
except ImportError:
    skip_oechem = support.can_skip("oe")

if not skip_oechem:
    from chemfp import openeye_patterns
    from chemfp.openeye import OEGRAPHSIM_API_VERSION

import test_patterns
from chemfp import types

MACCS_SMI = support.fullpath("maccs.smi")

def parse_smiles(smiles):
    mol = OEGraphMol()
    OEParseSmiles(mol, smiles)
    return mol

def _count(it):
    return sum(1 for item in it)

class ReferenceMixin(object):
    def test_reference_data_set(self):
        largest = min(self.reference_limit, max(v for (k,v) in self.reference_cases))

        matcher = self.reference_class(largest)
            
        for (smiles, expected) in self.reference_cases:
            mol = parse_smiles(smiles)

            self.assertEquals(matcher.SingleMatch(mol), bool(expected), smiles)

            expected = min(expected, largest)
            self.assertGreaterEqual(_count(matcher.Match(mol)), expected, smiles)

    def test_match_limit(self):
        largest = min(5, self.reference_limit)
        for max_count in range(1, largest+1):
            matcher = self.reference_class(max_count)
            for (smiles, expected) in self.reference_cases:
                mol = parse_smiles(smiles)
                expected = min(expected, max_count)
                self.assertGreaterEqual(_count(matcher.Match(mol)), expected, smiles)

class TestHydrogenMatcher(ReferenceMixin, unittest2.TestCase):
    if not skip_oechem:
        reference_class = openeye_patterns.HydrogenMatcher
        reference_cases = test_patterns.hydrogen_test_cases
        reference_limit = 100

TestHydrogenMatcher = unittest2.skipIf(skip_oechem, "OEChem not installed")(
    TestHydrogenMatcher)

class TestAromaticRingMatcher(ReferenceMixin, unittest2.TestCase):
    if not skip_oechem:
        reference_class = openeye_patterns.AromaticRings
        reference_cases = test_patterns.aromatic_ring_cases
        reference_limit = 2

TestAromaticRingMatcher = unittest2.skipIf(skip_oechem, "OEChem not installed")(
    TestAromaticRingMatcher)

# XXX These are too low-level. The tests should really be done through the
# file interface and shared with the other implementations.
 
class TestAromaticRings(unittest2.TestCase):
    def _count(self, matcher, smiles):
        return sum(1 for x in matcher.Match(parse_smiles(smiles)))
    def test_at_least_one_aromatic_ring(self):
        matcher = openeye_patterns.AromaticRings(1)
        self.assertFalse(matcher.SingleMatch(parse_smiles("C")))
        self.assertTrue(matcher.SingleMatch(parse_smiles("c1ccccc1")))
        self.assertTrue(matcher.SingleMatch(parse_smiles("c1ccccc1c1ccccc1")))

        self.assertEqual(self._count(matcher, "C"), 0)
        self.assertEqual(self._count(matcher, "c1ccccc1"), 1)
        self.assertIn(self._count(matcher, "c12cccccc1ccccc2"), (1, 2)) # Matcher can return >count; XXX Why?
        self.assertEqual(self._count(matcher, "c1ccccc1c1ccccc1"), 1)
        self.assertEqual(self._count(matcher, "c1ccccc1C1CCCCC1"), 1)

    def test_at_least_two_aromatic_rings(self):
        matcher = openeye_patterns.AromaticRings(2)

        self.assertEqual(self._count(matcher, "C"), 0)
        self.assertEqual(self._count(matcher, "c1ccccc1"), 1)
        self.assertEqual(self._count(matcher, "c12cccccc1ccccc2"), 2)
        self.assertEqual(self._count(matcher, "c1ccccc1c1ccccc1"), 2)
        self.assertEqual(self._count(matcher, "c1ccccc1C1CCCCC1"), 1)
        self.assertGreaterEqual(self._count(matcher, "c1ccc-2cccc-2cc1"), 1)
        self.assertGreaterEqual(self._count(matcher, "c1ccc-2cccc-2c3c1cccn3"), 2)
        self.assertGreaterEqual(self._count(matcher, "c1ccccc1.c1ccccc1"), 2)

    def test_aromatic_rings_failure(self):
        with self.assertRaisesRegexp(NotImplementedError, "No support for >=3 aromatic rings"):
            openeye_patterns.AromaticRings(3)

class TestHeteroAromaticRings(unittest2.TestCase):
    def _count(self, matcher, smiles):
        return sum(1 for x in matcher.Match(parse_smiles(smiles)))
    def test_at_least_one_heteroaromatic_ring(self):
        matcher = openeye_patterns.HeteroAromaticRings(1)
        self.assertFalse(matcher.SingleMatch(parse_smiles("C")))
        self.assertTrue(matcher.SingleMatch(parse_smiles("c1cccnc1")))
        self.assertTrue(matcher.SingleMatch(parse_smiles("c1cccnc1c1ccccc1")))
        self.assertTrue(matcher.SingleMatch(parse_smiles("c1ccccc1c1ncccc1")))

        self.assertEqual(self._count(matcher, "C"), 0)
        self.assertEqual(self._count(matcher, "c1ccccc1"), 0)
        self.assertEqual(self._count(matcher, "c1ccncc1"), 1)
        self.assertEqual(self._count(matcher, "c12cccccc1ccccc2"), 0)
        self.assertEqual(self._count(matcher, "c12ccnc1ccccc2"), 1)
        self.assertEqual(self._count(matcher, "c12cncc1ccccc2"), 1)
        self.assertEqual(self._count(matcher, "c12ccccc1cncc2"), 1)
        self.assertIn(self._count(matcher, "c12cccnc1cncc2"), (1, 2))  # a matcher can return >count. Why? XXX
        self.assertEqual(self._count(matcher, "c1ccccc1c1ccccc1"), 0)
        self.assertEqual(self._count(matcher, "c1ccccc1C1CCCNC1"), 0)
        self.assertEqual(self._count(matcher, "c1cnccc1C1CCCNC1"), 1)
        self.assertEqual(self._count(matcher, "c1cnccc1C1CCCNC1"), 1)

    def test_at_least_two_heteroaromatic_rings(self):
        matcher = openeye_patterns.HeteroAromaticRings(2)
        
        self.assertEqual(self._count(matcher, "C"), 0)
        self.assertEqual(self._count(matcher, "c1ccccc1"), 0)
        self.assertEqual(self._count(matcher, "c1ccncc1"), 1)
        self.assertEqual(self._count(matcher, "c12cccccc1ccccc2"), 0)
        self.assertEqual(self._count(matcher, "c12ccnc1ccccc2"), 1)
        self.assertEqual(self._count(matcher, "c12cncc1ccccc2"), 1)
        self.assertEqual(self._count(matcher, "c12ccccc1cncc2"), 1)
        self.assertEqual(self._count(matcher, "c12cccnc1cncc2"), 2)
        self.assertEqual(self._count(matcher, "c1ccccc1c1ccccc1"), 0)
        self.assertEqual(self._count(matcher, "c1ccccc1C1CCCNC1"), 0)
        self.assertEqual(self._count(matcher, "c1cnccc1C1CCCNC1"), 1)
        self.assertEqual(self._count(matcher, "c1cnccc1C1CCCNC1"), 1)

        self.assertEqual(self._count(matcher, "c1ccccn1.c1ccoc1"), 2)
        self.assertEqual(self._count(matcher, "c1ccc-2cccc-2cc1"), 0)
        self.assertEqual(self._count(matcher, "c1ccn-2ccc-2cc1"), 1)
        self.assertEqual(self._count(matcher, "c1ccccc2c1cncnc2"), 1)

    def test_heteroaromatic_rings_failure(self):
        with self.assertRaisesRegexp(NotImplementedError, "No support for >=3 hetero-aromatic rings"):
            openeye_patterns.HeteroAromaticRings(3)


class TestNumFragments(unittest2.TestCase):
    def _count(self, matcher, smiles):
        return sum(1 for x in matcher.Match(parse_smiles(smiles)))
    def test_single_match(self):
        matcher = openeye_patterns.NumFragments(10)
        self.assertEqual(matcher.SingleMatch(OEGraphMol()), 0)
        self.assertEqual(matcher.SingleMatch(parse_smiles("C")), 1)
        self.assertEqual(matcher.SingleMatch(parse_smiles("C.C")), 1)

    def test_multiple_matches(self):
        matcher = openeye_patterns.NumFragments(10)
        self.assertEqual(self._count(matcher, ""), 0)
        self.assertEqual(self._count(matcher, "C"), 1)
        self.assertEqual(self._count(matcher, "CC"), 1)
        self.assertEqual(self._count(matcher, "C.C"), 2)
        self.assertEqual(self._count(matcher, "C.C.O.[U]"), 4)

reason = None
if skip_oechem:
    reason = "OEChem not installed"
elif OEGRAPHSIM_API_VERSION == "1":
    reason = "OEGraphSim is not new enough"

TestAromaticRings = unittest2.skipIf(reason, reason)(TestAromaticRings)
TestHeteroAromaticRings = unittest2.skipIf(reason, reason)(TestHeteroAromaticRings)
TestNumFragments = unittest2.skipIf(reason, reason)(TestNumFragments)

if __name__ == "__main__":
    unittest2.main()
