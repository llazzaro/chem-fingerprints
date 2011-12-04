from __future__ import absolute_import

from openbabel import OBSmartsPattern, OBBitVec

from . import openbabel
from . import pattern_fingerprinter
from . import types
from . import __version__ as chemfp_version

SOFTWARE = openbabel.SOFTWARE + ("chemfp/%s" % (chemfp_version,))

# HasMatch was added to OpenBabel 2.3 .
_HAS_HASMATCH = hasattr(OBSmartsPattern(), "HasMatch")

class Matcher(object):
    def __init__(self, ob_matcher):
        self.ob_matcher = ob_matcher

    if _HAS_HASMATCH:
        def HasMatch(self, mol):
            return self.ob_matcher.HasMatch(mol)
    else:
        def HasMatch(self, mol):
            return self.ob_matcher.Match(mol, True)
        
    def NumUniqueMatches(self, mol):
        self.ob_matcher.Match(mol)
        return sum(1 for x in self.ob_matcher.GetUMapList())

class HydrogenMatcher(object):
    def __init__(self, max_count):
        self.max_count = max_count
    def HasMatch(self, mol):
        for i in range(1, mol.NumAtoms()+1):
            atom = mol.GetAtom(i)
            if atom.GetAtomicNum() == 1:
                return True
            if atom.ImplicitHydrogenCount():
                return True
        return False
        
    def NumUniqueMatches(self, mol):
        count = 0
        for i in range(1, mol.NumAtoms()+1):
            atom = mol.GetAtom(i)
            if atom.GetAtomicNum() == 1:
                count += 1
            count += atom.ImplicitHydrogenCount()
            if count >= self.max_count:
                return count
        return count

class AromaticRings(object):
    def __init__(self, max_count):
        self.max_count = max_count
        # The single ring case is easy; if there's an aromatic atom in a ring
        # then there's an aromatic ring.
        self._single_matcher = OBSmartsPattern()
        assert self._single_matcher.Init("[aR]")
        
    if _HAS_HASMATCH:
        def HasMatch(self, mol):
            return self._single_matcher.HasMatch(mol)
    else:
        def HasMatch(self, mol):
            return self._single_matcher.Match(mol, True)
        

    def NumUniqueMatches(self, mol):
        num_rings = 0
        for ring in mol.GetSSSR():
            if ring.IsAromatic():
                num_rings += 1
                if num_rings == self.max_count:
                    break
        return num_rings


###

class HeteroAromaticRings(object):
    def __init__(self, max_count):
        self.max_count = max_count
        # The single ring case is easy; if there's an non-carbon aromatic atom
        # in a ring then there's a hetero-aromatic ring
        self._single_matcher = OBSmartsPattern()
        assert self._single_matcher.Init("[aR;!#6]")

    if _HAS_HASMATCH:
        def HasMatch(self, mol):
            return self._single_matcher.HasMatch(mol)
    else:
        def HasMatch(self, mol):
            return self._single_matcher.Match(mol, True)
        

    def NumUniqueMatches(self, mol):
        num_rings = 0
        for ring in mol.GetSSSR():
            if (ring.IsAromatic() and
                  any(mol.GetAtom(atom_id).GetAtomicNum() != 6 for atom_id in ring._path)):
                num_rings += 1
                if num_rings == self.max_count:
                    break
        return num_rings


class NumFragments(object):
    def __init__(self, max_count):
        if max_count >= 3:
            raise NotImplemented("<fragment> >= 3 not implemented")
    def HasMatch(self, mol):
        # Has at least one fragment, which is the same as having atoms
        return (mol.NumAtoms() > 0)
    def NumUniqueMatches(self, mol):
        # If the largest fragment is not the same size as the entire molecule
        # then there are at least two fragments
        n = mol.NumAtoms()
        if n == 0:
            # No atoms means no fragments
            return 0
        
        v = OBBitVec()
        mol.FindLargestFragment(v)
        if v.CountBits() == n:
            return 1
        
        return 2

_pattern_classes = {
    "<H>": HydrogenMatcher,
    "<aromatic-rings>": AromaticRings,
    "<hetero-aromatic-rings>": HeteroAromaticRings,
    "<fragments>": NumFragments,
    }

def ob_compile_pattern(pattern, max_count):
    if pattern in _pattern_classes:
        return _pattern_classes[pattern](max_count)
    
    if pattern.startswith("<"):
        raise NotImplementedError(pattern)

    matcher = OBSmartsPattern()
    if not matcher.Init(pattern):
        raise pattern_fingerprinter.UnsupportedPatternError(
            pattern, "Uninterpretable SMARTS pattern")
    return Matcher(matcher)

class OBPatternFingerprinter(pattern_fingerprinter.PatternFingerprinter):
    def __init__(self, patterns):
        assert patterns is not None
        super(OBPatternFingerprinter, self).__init__(patterns, ob_compile_pattern)
        
    def fingerprint(self, mol):
        bytes = [0] * self.num_bytes
        for matcher, largest_count, count_info_tuple in self.matcher_definitions:
            #print matcher, largest_count, count_info_tuple
            if largest_count == 1:
                if matcher.HasMatch(mol):
                    count_info = count_info_tuple[0]
                    bytes[count_info.byteno] |= count_info.bitmask
            else:
                # There's no good way to get the number of unique
                # matches without iterating over all of them.
                actual_count = matcher.NumUniqueMatches(mol)
                for count_info in count_info_tuple:
                    if actual_count >= count_info.count:
                        bytes[count_info.byteno] |= count_info.bitmask
                    else:
                        break
        return "".join(map(chr, bytes))


class _CachedFingerprinters(dict):
    def __missing__(self, name):
        patterns = pattern_fingerprinter._load_named_patterns(name)
        fingerprinter = OBPatternFingerprinter(patterns)
        self[name] = fingerprinter
        return fingerprinter
_cached_fingerprinters = _CachedFingerprinters()

# XXX Why are there two "Fingerprinter" classes?
# XX Shouldn't they be merged?

_base = openbabel._base.clone(
    software = SOFTWARE)
    

SubstructOpenBabelFingerprinter_v1 = _base.clone(
    name = "ChemFP-Substruct-OpenBabel/1",
    num_bits = 881,
    make_fingerprinter = lambda: _cached_fingerprinters["substruct"].fingerprint)
    

RDMACCSOpenBabelFingerprinter_v1 = _base.clone(
    name = "RDMACCS-OpenBabel/1",
    num_bits = 166,
    make_fingerprinter = lambda: _cached_fingerprinters["rdmaccs"].fingerprint)
