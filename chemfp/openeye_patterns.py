from __future__ import absolute_import

from openeye.oechem import (
    OESubSearch, OEChemGetRelease, OEChemGetVersion, OEGraphMol, OEAndAtom,
    OENotAtom, OEIsAromaticAtom, OEIsCarbon, OEIsAromaticBond, OEAtomIsInRing, OEHasBondIdx,
    OEFindRingAtomsAndBonds, OEDetermineAromaticRingSystems, OEDetermineComponents)


from . import openeye
from . import pattern_fingerprinter
from . import types
from . import __version__ as chemfp_version
        
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

    def Match(self, mol, flg=True):
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
        self.max_count = max_count
        self._single_aromatic = OESubSearch("[aR]")
        # In OpenEye SMARTS, [a;!R2] will find aromatic atoms in at least two rings
        # The following finds atoms which are members of at least two aromatic rings
        self._multiring_aromatic = OESubSearch("[a;!R2](:a)(:a):a")
        
    def SingleMatch(self, mol):
        # This is easy; if there's one aromatic atom then there's one
        # aromatic ring.
        return self._single_aromatic.SingleMatch(mol)
        
    def Match(self, mol, flg=True):
        # We're trying to find if there are two aromatic rings.
        if not self._single_aromatic.SingleMatch(mol):
            return ()
            
        if self._multiring_aromatic.SingleMatch(mol):
            # then obviously there are two aromatic rings
            return (1,2)

        # Since no aromatic atom is in two aromatic rings that means
        # the aromatic ring systems are disjoint, so this gives me the
        # number of ring systems
        num_aromatic_systems, parts = OEDetermineAromaticRingSystems(mol)
        if num_aromatic_systems >= self.max_count:
            return [0]*self.max_count
        
        assert num_aromatic_systems != 0, "there is supposed to be an aromatic ring"
        if num_aromatic_systems == 1:
            return (1,)
        raise AssertionError("Should not get here")
    
_is_hetereo_aromatic = OEAndAtom(OEAndAtom(OEIsAromaticAtom(), OENotAtom(OEIsCarbon())), OEAtomIsInRing())
class HeteroAromaticRings(object):
    def __init__(self, max_count):
        if max_count > 2:
            raise NotImplementedError("No support for >=3 hetero-aromatic rings")
        self.max_count = max_count

    def SingleMatch(self, mol):
        for atom in mol.GetAtoms(_is_hetereo_aromatic):
            return True
        return False

    def Match(self, mol, flg=True):
        # Find all the hetero-aromatic atoms
        hetero_atoms = [atom for atom in mol.GetAtoms(_is_hetereo_aromatic)]
        if len(hetero_atoms) < 2:
            # The caller just needs an iterable
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

            for atom in newmol.GetAtoms(_is_hetereo_aromatic):
                return (1,2)

        return (1,)

class NumFragments(object):
    def __init__(self, max_count):
        pass
    def SingleMatch(self, mol):
        return mol.NumAtoms() > 0
    def Match(self, mol, flg=True):
        count, parts = OEDetermineComponents(mol)
        # parts is a list of component numbers.
        # Turn them into a set to get the unique set of component numbers
        # Sets are iterable, so I don't need to do more for the API
        return set(parts)

# Grrr. The substructure keys want up to 4 aromatic rings. The above
# code only works for up to 2. The API doesn't let me say "I can
# handle up to 2; please set the remainder to 0."
#
# XXX Well, I can change that.

def aromatic_rings(max_count):
    if max_count > 2:
        return pattern_fingerprinter.LimitedMatcher(2, AromaticRings(2))
    return AromaticRings(max_count)

def hetero_aromatic_rings(max_count):
    if max_count > 2:
        return NotImplemented
    return HeteroAromaticRings(max_count)

_pattern_classes = {
    "<H>": HydrogenMatcher,
    "<aromatic-rings>": aromatic_rings,
    "<hetero-aromatic-rings>": hetero_aromatic_rings,
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
            if matcher is NotImplemented:
                continue
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
        

SOFTWARE = "OEChem/%(release)s (%(version)s) chemfp/%(chemfp)s" % dict(
    release = OEChemGetRelease(),
    version = OEChemGetVersion(),
    chemfp = chemfp_version)


# XXX Why are there two "Fingerprinter" classes?
# XX Shouldn't they be merged?

_base = openeye._base.clone(
    software = SOFTWARE)

SubstructOpenEyeFingerprinter_v1 = _base.clone(
    name = "ChemFP-Substruct-OpenEye/1",
    num_bits = 881,
    make_fingerprinter = lambda : _cached_fingerprinters["substruct"].fingerprint)
    
#    def describe(self, bitno):
#        return self._fingerprinter.describe(bitno)

RDMACCSOpenEyeFingerprinter_v1 = _base.clone(
    name = "RDMACCS-OpenEye/1",
    num_bits = 166,
    make_fingerprinter = lambda : _cached_fingerprinters["rdmaccs"].fingerprint)
