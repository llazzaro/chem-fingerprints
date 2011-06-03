from __future__ import absolute_import

from indigo import Indigo, IndigoException

from . import indigo

# You may not share objects between different Indigo sessions. Rather
# than doing the work to isolate the objects so the Indigo() instance
# here can take molecules from another Indigo() instance, I assume
# that all molecules come from chemfp.indigo._indigo and just use that.
_indigo = indigo._indigo

from . import indigo
from . import pattern_fingerprinter
from . import types

SOFTWARE = indigo.SOFTWARE


# Indigo uses a "matcher" (derived from the molecule) and the "query"
# as part of its API. In use it's something like:
#    matcher = indigo.matcher(molecule)
#    for match in matcher.iterateMatches(query):
#        count += 1
# The matcher can hold (opaque) configuration about the
# molecule which is shared with the queries.

# Some queries, like "number of hydrogens", can't be done in SMARTS. I
# ended up making a new matching API for both kinds of cases.

class Matcher(object):
    def __init__(self, query):
        self.query = query
    def has_match(self, indigo_matcher, mol):
        for match in indigo_matcher.iterateMatches(self.query):
            return True
        return False

    def num_matches(self, indigo_matcher, mol, max_count):
        return indigo_matcher.countMatches(self.query)

# Get the hydrogen count

class HydrogenMatcher(object):
    def has_match(self, indigo_matcher, mol):
        count = mol.countImplicitHydrogens()
        if count:
            return True
        return any(atom.atomicNumber() == 1 for atom in mol.iterateAtoms())
    def num_matches(self, indigo_matcher, mol, max_count):
        # This should all be fast, at the C level
        count = mol.countImplicitHydrogens()
        if count >= max_count:
            return count
        # Perhaps there's a few more hanging around?
        for atom in mol.iterateAtoms():
            if atom.atomicNumber() == 1:
                count += 1
                if count >= max_count:
                    return count
        return count

# Use SSSR to count the number of aromatic rings

class AromaticRings(object):
    def __init__(self):
        self._single_query = _indigo.loadSmarts("[aR]")
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

# Use SSSR to count the number of hetero-aromatic rings

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
            # bond-order == 4 means "aromatic"; all rings bonds must be aromatic
            if all(bond.bondOrder() == 4 for bond in ring.iterateBonds()):
                # Check that there's at least one non-carbon atom in the ring
                if any(1 for atom in ring.iterateAtoms() if atom.atomicNumber() != 6):
                    count += 1
                    if count == max_count:
                        return count
        return count


# Get the number of components/fragments in the molecule

class NumFragments(object):
    def has_match(self, indigo_matcher, mol):
        return mol.countAtoms() > 0
    def num_matches(self, indigo_matcher, mol, max_count):
        return mol.countComponents()

######

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

# There's a bug in Indigo 1.0-beta13 where the __del__ doesn't work
# correctly during module shutdown. I've reported this and it will be
# in the next release. Until then, force a cleanup.
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
