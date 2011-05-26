from __future__ import absolute_import

from openbabel import OBSmartsPattern

from . import openbabel
from . import pattern_fingerprinter

SOFTWARE = openbabel.SOFTWARE

class HydrogenMatcher(object):
    def Match(self, mol):
        count = 0
        for i in range(1, mol.NumAtoms()+1):
            atom = mol.GetAtom(i)
            if atom.GetAtomicNum() == 1:
                count += 1
            count += atom.ImplicitHydrogenCount()

        self.count = count
        return self.count  # The OB API returns 0 or non-zero
        
    def NumMatches(self):
        return self.count


def ob_compile_pattern(pattern, max_count):
    if pattern == "<H>":
        return HydrogenMatcher()
    
    if pattern.startswith("<"):
        return NotImplemented

    matcher = OBSmartsPattern()
    if not matcher.Init(pattern):
        raise pattern_fingerprinter.UnsupportedPatternError(
            pattern, "Uninterpretable SMARTS pattern")
    return matcher

class OBPatternFingerprinter(pattern_fingerprinter.PatternFingerprinter):
    def __init__(self, patterns):
        assert patterns is not None
        super(OBPatternFingerprinter, self).__init__(patterns, ob_compile_pattern)
        
    def fingerprint(self, mol):
        bytes = [0] * self.num_bytes
        for matcher, largest_count, count_info_tuple in self.matcher_definitions:
            #print matcher, largest_count, count_info_tuple
            if largest_count == 1:
                if matcher.Match(mol, True):
                    count_info = count_info_tuple[0]
                    bytes[count_info.byteno] |= count_info.bitmask
            else:
                if matcher.Match(mol):
                    actual_count = matcher.NumMatches()
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


def _read_fingerprints(pattern_name, source, format, kwargs):
    assert not kwargs
    fingerprinter = _cached_fingerprinters[pattern_name].fingerprint
    structure_reader = openbabel.read_structures(source, format)

    def read_pattern_fingerprints():
        for (title, mol) in structure_reader:
            yield fingerprinter(mol), title

    return read_pattern_fingerprints()
    

def read_substruct_fingerprints_v1(source=None, format=None, kwargs={}):
    return _read_fingerprints("substruct", source, format, kwargs)

def read_rdmaccs_fingerprints_v1(source=None, format=None, kwargs={}):
    return _read_fingerprints("rdmaccs", source, format, kwargs)
    
