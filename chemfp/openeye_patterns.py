from __future__ import absolute_import

from openeye.oechem import (
    OESubSearch, OEChemGetRelease, OEChemGetVersion, OEGraphMol, OEAndAtom,
    OENotAtom, OEIsAromaticAtom, OEIsCarbon, OEIsAromaticBond, OEHasBondIdx,
    OEFindRingAtomsAndBonds, OEDetermineAromaticRingSystems)
                            
                            


from . import openeye
from . import pattern_fingerprinter
        
class HydrogenMatcher(object):
    def __init__(self, max_count):
        self.max_count = max_count

    def SingleMatch(self, mol):
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 1:
                return 1
            if atom.GetImplicitHCount():
                return 1
        return 0

    def Matches(self, mol, flg=True):
        max_count = self.max_count
        count = 0
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 1:
                count += 1
            count += atom.GetImplicitHCount()
            if count > max_count:
                break
        return [0] * count

# OpenEye famously does not include SSSR functionality in OEChem.
# Search for "Smallest Set of Smallest Rings (SSSR) Considered Harmful"
# After much thought, I agree. But it makes this sort of code harder.

# That's why I only support up to max_count = 2. Then again, I know
# that this code does the right thing, while I'm not sure about the
# SSSR-based implementations.

class AromaticRings(object):
    def __init__(self, max_count):
        if max_count > 2:
            raise NotImplementedError("No support for >=3 aromatic rings")
        self._single_aromatic = OESubSearch("[a]")
        # In OpenEye SMARTS, [a;!R2] will find aromatic atoms in at least two rings
        # The following finds atoms which are members of at least two aromatic rings
        self._multiring_aromatic = OESubSearch("[a;!R2](:a)(:a):a")
        
    def SingleMatch(self, mol):
        # This is easy; if there's one aromatic atom then there's one
        # aromatic ring.
        return self._single_aromatic.SingleMatch(mol)
        
    def Matches(self, mol, flg=True):
        # We're trying to find if there are two aromatic rings.
        if self._multiring_aromatic.SingleMatch(mol):
            # then obviously there are two aromatic rings
            return (1,2)

        # Since no aromatic atom is in two aromatic rings that means
        # the aromatic ring systems are disjoint, so this gives me the
        # number of ring systems
        num_aromatic_systems, parts = OEDetermineAromaticRingSystems(mol)
        if num_aromatic_systems == 0:
            return ()
        if num_aromatic_systems == 1:
            return (1,)
        return (1,2)

_is_hetereo_aromatic = OEAndAtom(OEIsAromaticAtom(), OENotAtom(OEIsCarbon()))
class HeteroAromaticRings(object):
    def __init__(self, max_count):
        if max_count > 2:
            raise NotImplementedError("No support for >=3 hetero-aromatic rings")
        self._single_hetero_aromatic = OESubSearch("[a;!#6]")
        self._single_hetero_aromatic.SetMaxMatches(1)

    def SingleMatch(self, mol):
        return self._single_hetero_aromatic.SingleMatch(mol)

    def Matches(self, mol, flg=True):
        # Find all the hetero-aromatic atoms
        hetero_atoms = [atom for atom in mol.GetAtoms(_is_hetereo_aromatic)]
        if len(hetero_atoms) < 2:
            # I just need an iterable
            return hetero_atoms

        # There are at least two hetero-aromatic atoms.
        # Are there multiple ring systems?
        num_aromatic_systems, parts = OEDetermineAromaticRingSystems(mol)
        assert num_aromatic_systems >= 1
        # Are there hetero-atoms in different systems?
        atom_components = set(parts[atom.GetIdx()] for atom in hetero_atoms)
        if len(atom_components) > 1:
            return (1,2)

        # The answer now is "at least one". But are there two?

        # All of the hetero-aromatic atoms are in the same ring system
        # This is the best answer I could think of, and it only works
        # with the OEChem toolkit: remove one of the bonds, re-find
        # the rings, and see if there's still an aromatic hetero-atom.
        hetero_atom = hetero_atoms[0]

        for bond in hetero_atom.GetBonds(OEIsAromaticBond()):
            newmol = OEGraphMol(mol)
            newmol_bond = newmol.GetBond(OEHasBondIdx(bond.GetIdx()))
            newmol.DeleteBond(newmol_bond)
            OEFindRingAtomsAndBonds(newmol)

            if self._single_hetero_aromatic.SingleMatch(newmol):
                return (1,2)

        return (1,)

class NumFragments(object):
    def __init__(self, max_count):
        pass
    def SingleMatch(self, mol):
        count, parts = OEDetermineComponents(mol)
        return count > 0
    def Matches(self, mol, flg=True):
        count, parts = OEDetermineComponents(mol)
        return parts[1:]
        
        

_pattern_classes = {
    "<H>": HydrogenMatcher,
    "<aromatic-rings>": AromaticRings,
    "<hetero-aromatic-rings>": HeteroAromaticRings,
    "<fragments>": NumFragments,
    }
    

def oechem_compile_pattern(pattern, max_count):
    if pattern in _pattern_classes:
        return _pattern_classes[pattern](max_count)

    elif pattern.startswith("<"):
        raise NotImplementedError(pattern) # No other special patterns are supported
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
