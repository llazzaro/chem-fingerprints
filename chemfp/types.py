# Information about fingerprint types

from . import io
from .decoders import import_decoder  # XXX too specific to the decoder module


def check_openbabel_maccs166():
    from chemfp.openbabel import HAS_MACCS, MACCS_VERSION
    assert HAS_MACCS
    if MACCS_VERSION == 1:
        return "OpenBabel-MACCS/1"
    elif MACCS_VERSION == 2:
        return "OpenBabel-MACCS/2"
    raise AssertionError

class FingerprintFamily(object):
    def __init__(self, family_name, path):
        self.family_name = family_name
        self.path = path

    def __call__(self, **kwargs):
        cls = import_decoder(self.path)
        assert cls.name == self.family_name
        return cls(kwargs)

    def from_parameters(self, parameters):
        cls = import_decoder(self.path)
        assert cls.name == self.family_name
        return cls.from_parameters(parameters)

_families = [
    FingerprintFamily("OpenEye-MACCS166/1", "chemfp.openeye.OpenEyeMACCSFingerprinter_v1"),
    FingerprintFamily("OpenEye-Path/1", "chemfp.openeye.OpenEyePathFingerprinter_v1"),
    
    FingerprintFamily("RDKit-MACCS166/1", "chemfp.rdkit.RDKitMACCSFingerprinter_v1"),
    FingerprintFamily("RDKit-Fingerprint/1", "chemfp.rdkit.RDKitFingerprinter_v1"),
    
    FingerprintFamily("OpenBabel-FP2/1", "chemfp.openbabel.OpenBabelFP2Fingerprinter_v1"),
    FingerprintFamily("OpenBabel-FP3/1", "chemfp.openbabel.OpenBabelFP3Fingerprinter_v1"),
    FingerprintFamily("OpenBabel-FP4/1", "chemfp.openbabel.OpenBabelFP4Fingerprinter_v1"),
    FingerprintFamily("OpenBabel-MACCS/1", "chemfp.openbabel.OpenBabelMACCSFingerprinter_v1"),
    FingerprintFamily("OpenBabel-MACCS/2", "chemfp.openbabel.OpenBabelMACCSFingerprinter_v2"),

    # In the future this will likely change to use a parameterized class
    # which can dynamically load fingerprint definitions

    FingerprintFamily("ChemFP-Substruct-OpenEye/1",
                      "chemfp.openeye_patterns.SubstructOpenEyeFingerprinter_v1"),
    FingerprintFamily("RDMACCS-OpenEye/1",
                      "chemfp.openeye_patterns.RDMACCSOpenEyeFingerprinter_v1"),

    FingerprintFamily("ChemFP-Substruct-RDKit/1",
                      "chemfp.rdkit_patterns.SubstructRDKitFingerprinter_v1"),
    FingerprintFamily("RDMACCS-RDKit/1",
                      "chemfp.rdkit_patterns.RDMACCSRDKitFingerprinter_v1"),

    FingerprintFamily("ChemFP-Substruct-OpenBabel/1",
                      "chemfp.openbabel_patterns.SubstructOpenBabelFingerprinter_v1"),
    FingerprintFamily("RDMACCS-OpenBabel/1",
                      "chemfp.openbabel_patterns.RDMACCSOpenBabelFingerprinter_v1"),
]

_alternates = {
    "OpenBabel-MACCS": check_openbabel_maccs166
    }


_family_by_name = {}

def _initialize_families():
    for family in _families:
        # Set both the versioned and non-versioned names
        name = family.family_name
        unversioned_name = name.split("/")[0]
        _family_by_name[name] = _family_by_name[unversioned_name] = family

    # Don't include a (likely non-versioned) name if there's a selector function
    for name in _alternates:
        if name in _family_by_name:
            del _family_by_name[name]

_initialize_families()


def get_fingerprint_family(name):
    try:
        return _family_by_name[name]
    except KeyError:
        if name not in _alternates:
            raise
    alternate = _alternates[name]()
    return _family_by_name[alternate]

class Fingerprinter(object):
    format_string = None
    software = None
    def __init__(self, kwargs):
        self.kwargs = kwargs
        # Some self-test code to make sure preconditions are met
        # This means they must be set before calling super().__init__
        if getattr(self, "name", None) is None:
            raise AssertionError("num_bits not defined (%r)" % (self.__class__,))
        if getattr(self, "num_bits", None) is None:
            raise AssertionError("num_bits not defined (%r)" % (self.name,))
        if getattr(self, "num_bits", None) is None:
            raise AssertionError("num_bits not defined (%r)" % (self.name,))
        if getattr(self, "software", None) is None:
            raise AssertionError("software not defined (%r)" % (self.name,))
        

    @classmethod
    def from_parameters(cls, parameters):
        if parameters:
            raise AssertionError
        return cls({})


    # Subclasses may hook into this
    def _encode_parameters(self):
        # Assume they can be interpreted directly in the format_string
        return self.kwargs

    def get_type(self):
        if self.format_string is None:
            if self.kwargs:
                assert not kwargs, self
            return self.name
        encoded = self.format_string % self._encode_parameters()
        return self.name + " " + encoded

    # Subclasses must hook into this
    def _get_reader(self, source, format, kwargs):
        raise NotImplementedError
    
    def read_structure_fingerprints(self, source, format=None):
        reader = self._get_reader(source, format, self.kwargs)

        source_filename = io.get_filename(source)
        return io.FPIterator(io.Header(num_bits = self.num_bits,
                                       source = source_filename,
                                       software = self.software,
                                       type = self.get_type(),
                                       date = io.utcnow()),
                             reader)
    def describe(self, bitno):
        if 0 <= bitno < self.num_bits:
            return "bit %d (unknown)" % (bitno,)
        raise KeyError(bitno)
        

def parse_type(type):
    terms = type.split()
    if not terms:
        raise TypeError("missing name")

    name = terms[0]
    try:
        family = get_fingerprint_family(name)
    except KeyError:
        raise TypeError(name)

    seen = set()
    parameters = []
    for term in terms[1:]:
        left, right = term.split("=")
        if left in seen:
            raise TypeError("Duplicate name")
        seen.add(left)
        parameters.append((left, right))

    return family.from_parameters(parameters)


def read_structure_fingerprints(type, source=None, format=None):
    structure_fingerprinter = parse_type(type)
    return structure_fingerprinter.read_structure_fingerprints(source, format)
