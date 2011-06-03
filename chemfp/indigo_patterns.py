from __future__ import absolute_import

from indigo import Indigo, IndigoException

from . import indigo

# Yech! It looks like I have to use the same Indigo() instance
# used to make the structure in the first place.
_indigo = indigo._indigo

from . import indigo
from . import pattern_fingerprinter
from . import types

SOFTWARE = indigo.SOFTWARE

class Matcher(object):
    def __init__(self, query):
        self.query = query
    def has_match(self, indigo_matcher, mol):
        for match in indigo_matcher.iterateMatches(self.query):
            return True
        return False

    def num_matches(self, indigo_matcher, mol, max_count):
        return indigo_matcher.countMatches(self.query)

class HydrogenMatcher(object):
    def has_match(self, indigo_matcher, mol):
        for atom in mol.iterateAtoms():
            if atom.countImplicitHydrogens():
                return True
            if atom.atomicNumber() == 1:
                return True
        return False
    def num_matches(self, indigo_matcher, mol, max_count):
        count = 0
        for atom in mol.iterateAtoms():
            count += atom.countImplicitHydrogens()
            if atom.atomicNumber() == 1:
                count += 1
            if count >= max_count:
                return count
        return count

class AromaticRings(object):
    def __init__(self):
        self._single_query = _indigo.loadSmarts("a")
    def has_match(self, indigo_matcher, mol):
        for match in indigo_matcher.match(self._single_query):
            return True
        return False
    def num_matches(self, indigo_matcher, mol, max_count):
        count = 0
        for ring in mol.iterateSSSR():
            if all(bond.bondOrder() == 4 for bond in ring.iterateBonds()):
                count += 1
                if count == max_count:
                    return count
        return count

class HeteroAromaticRings(object):
    def __init__(self):
        self._single_query = _indigo.loadSmarts("[a;!#6]")
    def has_match(self, indigo_matcher, mol):
        for match in indigo_matcher.match(self._single_query):
            return True
        return False
    def num_matches(self, indigo_matcher, mol, max_count):
        count = 0
        for ring in mol.iterateSSSR():
            if all(bond.bondOrder() == 4 for bond in ring.iterateBonds()):
                if any(1 for atom in ring.iterateAtoms() if atom.atomicNumber() != 6):
                    count += 1
                    if count == max_count:
                        return count
        return count

class NumFragments(object):
    def has_match(self, indigo_matcher, mol):
        return mol.countAtoms() > 0
    def num_matches(self, indigo_matcher, mol, max_count):
        # Looks like I need to do this the hard way
        if mol.countAtoms() == 0:
            return 0
        return mol.smiles().count(".") + 1
        

_pattern_classes = {
    "<H>": HydrogenMatcher,
    "<aromatic-rings>": AromaticRings,
    "<hetero-aromatic-rings>": HeteroAromaticRings,
    "<fragments>": NumFragments,
    }

def indigo_compile_pattern(pattern, max_count):
    if pattern in _pattern_classes:
        return _pattern_classes[pattern]()

    if pattern.startswith("<"):
        raise NotImplementedError(pattern)

    try:
        return Matcher(_indigo.loadSmarts(pattern))
    except IndigoException, err:
        raise pattern_fingerprinter.UnsupportedPatternError(
            pattern, "Uninterpretable SMARTS pattern: %s" % (err,))

class IndigoPatternFingerprinter(pattern_fingerprinter.PatternFingerprinter):
    def __init__(self, patterns):
        assert patterns is not None
        super(IndigoPatternFingerprinter, self).__init__(patterns, indigo_compile_pattern)

    def fingerprint(self, mol):
        bytes = [0] * self.num_bytes
        indigo_matcher = _indigo.substructureMatcher(mol)
        for matcher, largest_count, count_info_tuple in self.matcher_definitions:
            if largest_count == 1:
                if matcher.has_match(indigo_matcher, mol):
                    count_info = count_info_tuple[0]
                    bytes[count_info.byteno] |= count_info.bitmask
            else:
                actual_count = matcher.num_matches(indigo_matcher, mol, largest_count)
                if not actual_count:
                    continue
                for count_info in count_info_tuple:
                    if actual_count >= count_info.count:
                        bytes[count_info.byteno] |= count_info.bitmask
                    else:
                        break
        return "".join(map(chr, bytes))
        


class _CachedFingerprinters(dict):
    def __missing__(self, name):
        patterns = pattern_fingerprinter._load_named_patterns(name)
        fingerprinter = IndigoPatternFingerprinter(patterns)
        self[name] = fingerprinter
        return fingerprinter
_cached_fingerprinters = _CachedFingerprinters()

import atexit
def _remove():
    _cached_fingerprinters.clear()
atexit.register(_remove)

# XXX Why are there two "Fingerprinter" classes?
# XX Shouldn't they be merged?

class _PatternFingerprinter(types.Fingerprinter):
    software = SOFTWARE
    def __init__(self, kwargs):
        self._fingerprinter = _cached_fingerprinters[self._pattern_name]
        
        super(_PatternFingerprinter, self).__init__(kwargs)

    def fingerprint(self, mol):
        return self._fingerprinter(mol)

    def describe(self, bitno):
        return self._fingerprinter.describe(bitno)

    def _get_reader(self, source=None, format=None, kwargs={}):
        assert not kwargs
        fingerprinter = _cached_fingerprinters[self._pattern_name].fingerprint
        structure_reader = indigo.read_structures(source, format)

        def read_pattern_fingerprints(fingerprinter, structure_reader):
            for (title, mol) in structure_reader:
                yield fingerprinter(mol), title

        return read_pattern_fingerprints(fingerprinter, structure_reader)


class SubstructIndigoFingerprinter_v1(_PatternFingerprinter):
    name = "ChemFP-Substruct-Indigo/1"
    num_bits = 881
    _pattern_name = "substruct"
    software = SOFTWARE


class RDMACCSIndigoFingerprinter_v1(_PatternFingerprinter):
    name = "RDMACCS-Indigo/1"
    num_bits = 166
    _pattern_name = "rdmaccs"
    software = SOFTWARE
