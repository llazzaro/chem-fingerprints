from __future__ import absolute_import

from openeye.oechem import OESubSearch, OEChemGetRelease, OEChemGetVersion

from . import openeye
from . import pattern_fingerprinter
        
class HydrogenMatcher(object):
    pattern = "<H>"
    def __init__(self):
        self.SetMaxMatches(1024)
    def SetMaxMatches(self, max_count):
        self.max_count = max_count
    def Match(self, mol, flg):
        max_count = self.max_count
        count = 0
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 1:
                count += 1
            count += atom.GetImplicitHCount()
            if count > max_count:
                break
        return [0] * count

    def SingleMatch(self, mol):
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 1:
                return 1
            if atom.GetImplicitHCount():
                return 1
        return 0

class TrueMatcher(object):
    def __init__(self):
        self.SetMaxMatches(1024)
    def SetMaxMatches(self, max_count):
        self.max_count = max_count
        self._matches = [0] * max_count
    def Match(self, mol, flg):
        self._matches
    def SingleMatch(self, mol):
        return 1

def oechem_compile_pattern(pattern, max_count):
    if pattern == "<H>":
        return HydrogenMatcher()
    if pattern == "<1>":
        return TrueMatcher(max_count)
    if pattern == "<0>":
        return 0
    elif pattern.startswith("<"):
        return NotImplemented # No other special patterns are supported; set to 0
    else:
        pat = OESubSearch()
        if not pat.Init(pattern):
            raise pattern_fingerprinter.UnsupportedPatternError(
                pattern, "Uninterpretable SMARTS pattern")
        pat.SetMaxMatches(max_count)
        return pat


class OEChemPatternFingerprinter(pattern_fingerprinter.PatternFingerprinter):
    def __init__(self, patterns):
        assert patterns is not None
        super(OEChemPatternFingerprinter, self).__init__(patterns, oechem_compile_pattern)
        
    def fingerprint(self, mol):
        bytes = [0] * self.num_bytes
        for matcher, largest_count, count_info_tuple in self.matcher_definitions:
            #print matcher, largest_count, count_info_tuple
            if largest_count == 1:
                if matcher.SingleMatch(mol):
                    count_info = count_info_tuple[0]
                    bytes[count_info.byteno] |= count_info.bitmask
            else:
                actual_count = sum(1 for ignore in matcher.Match(mol, True)) # unique matches
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
        fingerprinter = OEChemPatternFingerprinter(patterns)
        self[name] = fingerprinter
        return fingerprinter
_cached_fingerprinters = _CachedFingerprinters()
        
#def load_substruct_fingerprinter_v1():
#    return _cached_fingerprinters["substruct"]
#
#def load_rdmaccs_fingerprinter_v1():
#    return _cached_fingerprinters["rdmaccs"]

def _read_fingerprints(pattern_name, source, format, kwargs):
    assert not kwargs
    # The OEChem interface only handles stdin and filenames
    if not (isinstance(source, basestring) or source is None):  # Why is the check here?
        raise NotImplementedError

    fingerprinter = _cached_fingerprinters[pattern_name].fingerprint
    structure_reader = openeye.read_structures(source, format)

    def read_pattern_fingerprints():
        for (title, mol) in structure_reader:
            yield fingerprinter(mol), title

    return read_pattern_fingerprints()
    

def read_substruct_fingerprints_v1(source=None, format=None, kwargs={}):
    return _read_fingerprints("substruct", source, format, kwargs)

def read_rdmaccs_fingerprints_v1(source=None, format=None, kwargs={}):
    return _read_fingerprints("rdmaccs", source, format, kwargs)


# XXX include ChemFP version information? Probably
SOFTWARE = "OEChem/%(release)s (%(version)s)" % dict(
    release = OEChemGetRelease(),
    version = OEChemGetVersion())
