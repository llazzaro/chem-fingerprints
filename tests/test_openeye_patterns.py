import sys
import unittest2

from openeye.oechem import OEGraphMol, OEParseSmiles


import test_patterns

from chemfp import openeye_patterns

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
    reference_class = openeye_patterns.HydrogenMatcher
    reference_cases = test_patterns.hydrogen_test_cases
    reference_limit = 100


class TestAromaticRingMatcher(ReferenceMixin, unittest2.TestCase):
    reference_class = openeye_patterns.AromaticRings
    reference_cases = test_patterns.aromatic_ring_cases
    reference_limit = 2

if __name__ == "__main__":
    unittest2.main()
