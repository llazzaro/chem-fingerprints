import sys
import unittest2
import support

try:
    from openeye.oechem import OEGraphMol, OEParseSmiles
    skip_oechem = False
except ImportError:
    skip_oechem = support.can_skip("oe")
else:
    from chemfp import openeye_patterns

import test_patterns


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

if __name__ == "__main__":
    unittest2.main()
