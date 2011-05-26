from __future__ import absolute_import

from rdkit import Chem

from . import pattern_fingerprinter
from . import rdkit

SOFTWARE = rdkit.SOFTWARE

class HydrogenMatcher(object):
    pass

class InvertedMatcher(object):
    def __init__(self, matcher):
        self.matcher = matcher
    def has_match(self, mol):
        return mol.HasSubstructMatch(self.matcher)
    def num_matches(self, mol):
        return len(mol.GetSubstructMatches(self.matcher))

def rdkit_compile_pattern(pattern, max_count):
    if pattern == "<H>":
        return HydrogenMatcher()

    elif pattern.startswith("<"):
        return NotImplemented

    matcher = Chem.MolFromSmarts(pattern)
    if matcher is None:
        raise pattern_fingerprinter.UnsupportedPatternError(
            pattern, "Can not interpret SMARTS pattern")
    return InvertedMatcher(matcher)

class RDKitPatternFingerprinter(pattern_fingerprinter.PatternFingerprinter):
    def __init__(self, patterns):
        assert patterns is not None
        super(RDKitPatternFingerprinter, self).__init__(patterns, rdkit_compile_pattern)
        
    def fingerprint(self, mol):
        bytes = [0] * self.num_bytes
        for matcher, largest_count, count_info_tuple in self.matcher_definitions:
            if largest_count == 1:
                if matcher.has_match(mol):
                    count_info = count_info_tuple[0]
                    bytes[count_info.byteno] |= count_info.bitmask
            else:
                actual_count = matcher.num_matches(mol)
                if actual_count:
                    for count_info in count_info_tuple:
                        if actual_count >= count_info.count:
                            bytes[count_info.byteno] |= count_info.bitmask
                        else:
                            break
        return "".join(map(chr, bytes))

class _CachedFingerprinters(dict):
    def __missing__(self, name):
        patterns = pattern_fingerprinter._load_named_patterns(name)
        fingerprinter = RDKitPatternFingerprinter(patterns)
        self[name] = fingerprinter
        return fingerprinter
_cached_fingerprinters = _CachedFingerprinters()

def _read_fingerprints(pattern_name, source, format, kwargs):
    assert not kwargs
    # The OEChem interface only handles stdin and filenames
    if not (isinstance(source, basestring) or source is None):  # Why is the check here?
        raise NotImplementedError

    fingerprinter = _cached_fingerprinters[pattern_name].fingerprint
    structure_reader = rdkit.read_structures(source, format)

    def read_pattern_fingerprints():
        for (title, mol) in structure_reader:
            yield fingerprinter(mol), title

    return read_pattern_fingerprints()
    

def read_substruct_fingerprints_v1(source=None, format=None, kwargs={}):
    return _read_fingerprints("substruct", source, format, kwargs)

def read_rdmaccs_fingerprints_v1(source=None, format=None, kwargs={}):
    return _read_fingerprints("rdmaccs", source, format, kwargs)

