# Information about fingerprint types

from . import io

def _convert_parameters(parameters, converters):
    kwargs = {}
    for (name, value) in parameters:
        convert = converters[name]
        kwargs[name] = convert(value)
    return kwargs

class _Opener(object):
    def _open(self, software, source, iterator):
        if source is not None:
            if not isinstance(source, basestring):
                # Then it's a Python file object
                source = getattr(source, "name", None)

        return io.FPIterator(io.Header(num_bits = self.num_bits,
                                       source = source,
                                       software = software,
                                       type = self.get_type(),
                                       date = io.utcnow()),
                             iterator)
    
class _NoParameters(_Opener):
    @classmethod
    def from_parameters(cls, parameters):
        assert len(parameters) == 0
        return cls()

    def get_type(self):
        return self.name
    
class OpenEyeMACCS166_v1(_NoParameters):
    name = "OpenEye-MACCS166/1"
    num_bits = 166

    def read_structure_fingerprints(self, source=None, format=None):
        from chemfp.openeye import read_maccs166_fingerprints_v1, SOFTWARE
        reader = read_maccs166_fingerprints_v1(source, format)
        return self._open(SOFTWARE, source, reader)

OpenEyeMACCS166 = OpenEyeMACCS166_v1

class OpenEyePath_v1(_Opener):
    name = "OpenEye-Path/1"
    format_string = ("OpenEye-Path/1 numbits=%(numbits)s minbonds=%(minbonds)s "
                     "maxbonds=%(maxbonds)s atype=%(atype)s btype=%(btype)s")
    
    converters = {"numbits": int,
                  "minbonds": int,
                  "maxbonds": int}

    def __init__(self, kwargs):
        assert len(kwargs) == 5, kwargs
        self.kwargs = kwargs
        self.num_bits = kwargs["numbits"]
    
    @staticmethod
    def from_parameters(parameters):
        converters = OpenEyePath_v1.converters
        if len(converters) == 3:
            from chemfp.openeye import atom_description_to_value, bond_description_to_value
            converters["atype"] = atom_description_to_value
            converters["btype"] = bond_description_to_value
        
        return OpenEyePath_v1(_convert_parameters(parameters, converters))
    
        
    def get_type(self):
        from chemfp.openeye import atom_value_to_description, bond_value_to_description
        kw = self.kwargs
        return self.format_string % dict(numbits = kw["numbits"],
                                         minbonds = kw["minbonds"],
                                         maxbonds = kw["maxbonds"],
                                         atype = atom_value_to_description(kw["atype"]),
                                         btype = bond_value_to_description(kw["btype"]))
                                             

    def read_structure_fingerprints(self, source=None, format=None):
        from chemfp.openeye import read_path_fingerprints_v1, SOFTWARE
        reader = read_path_fingerprints_v1(source, format, self.kwargs)
        return self._open(SOFTWARE, source, reader)

OpenEyePath = OpenEyePath_v1

class RDKitMACCS166_v1(_NoParameters):
    name = "RDKit-MACCS166/1"
    num_bits = 167

    def read_structure_fingerprints(self, source=None, format=None):
        from chemfp.rdkit import read_maccs166_fingerprints_v1, SOFTWARE
        reader = read_maccs166_fingerprints_v1(source, format)
        return self._open(SOFTWARE, source, reader)

RDKitMACCS166 = RDKitMACCS166_v1

class RDKitFingerprint_v1(_Opener):
    name = "RDKit-Fingerprint/1"
    format_string = ("RDKit-Fingerprint/1 minPath=%(minPath)s maxPath=%(maxPath)s fpSize=%(fpSize)s "
                     "nBitsPerHash=%(nBitsPerHash)s useHs=%(useHs)s")

    converters = {"minPath": int,
                  "maxPath": int,
                  "fpSize": int,
                  "nBitsPerHash": int,
                  "useHs": int}
    def __init__(self, kwargs):
        assert len(kwargs) == 5
        self.num_bits = kwargs["fpSize"]
        self.kwargs = kwargs

    @staticmethod
    def from_parameters(parameters):
        return RDKitFingerprint_v1(_convert_parameters(parameters, RDKitFingerprint_v1.converters))

    def get_type(self):
        return self.format_string % self.kwargs

    def read_structure_fingerprints(self, source=None, format=None):
        from chemfp.rdkit import read_rdkit_fingerprints_v1, SOFTWARE
        reader = read_rdkit_fingerprints_v1(source, format, self.kwargs)
        return self._open(SOFTWARE, source, reader)

RDKitFingerprint = RDKitFingerprint_v1

class OpenBabelFP2_v1(_NoParameters):
    name = "OpenBabel-FP2/1"
    num_bits = 1021

    def read_structure_fingerprints(self, source=None, format=None):
        from chemfp.openbabel import read_fp2_fingerprints_v1, SOFTWARE
        reader = read_fp2_fingerprints_v1(source, format)
        return self._open(SOFTWARE, source, reader)

OpenBabelFP2 = OpenBabelFP2_v1

class OpenBabelFP3_v1(_NoParameters):
    name = "OpenBabel-FP3/1"
    num_bits = 55

    def read_structure_fingerprints(self, source=None, format=None):
        from chemfp.openbabel import read_fp3_fingerprints_v1, SOFTWARE
        reader = read_fp3_fingerprints_v1(source, format)
        return self._open(SOFTWARE, source, reader)

OpenBabelFP3 = OpenBabelFP3_v1

class OpenBabelFP4_v1(_NoParameters):
    name = "OpenBabel-FP4/1"
    num_bits = 307

    def read_structure_fingerprints(self, source=None, format=None):
        from chemfp.openbabel import read_fp4_fingerprints_v1, SOFTWARE
        reader = read_fp4_fingerprints_v1(source, format)
        return self._open(SOFTWARE, source, reader)

OpenBabelFP4 = OpenBabelFP4_v1

class OpenBabelMACCS166_v1(_NoParameters):
    name = "OpenBabel-MACCS/1"
    num_bits = 166

    def read_structure_fingerprints(self, source=None, format=None):
        from chemfp.openbabel import read_maccs166_fingerprints_v1, SOFTWARE
        reader = read_maccs166_fingerprints_v1(source, format)
        return self._open(SOFTWARE, source, reader)


class OpenBabelMACCS166_v2(_NoParameters):
    name = "OpenBabel-MACCS/2"
    num_bits = 166

    def read_structure_fingerprints(self, source=None, format=None):
        from chemfp.openbabel import read_maccs166_fingerprints_v2, SOFTWARE
        reader = read_maccs166_fingerprints_v2(source, format)
        return self._open(SOFTWARE, source, reader)

def OpenBabelMACCS166():
    from chemfp.openbabel import HAS_MACCS, MACCS_VERSION
    assert HAS_MACCS
    if MACCS_VERSION == 1:
        return OpenBabelMACCS166_v1()
    elif MACCS_VERSION == 2:
        return OpenBabelMACCS166_v2()
    raise AssertionError

### Substructure fingerprints

def check_hydrogens(name):
    if name not in ("none", "implicit", "all"):
        raise TypeError
    return name

_substructure_converters = {
    "hydrogens": check_hydrogens
    }

class _Substruct_v1(_Opener):
    num_bits = 881
    def __init__(self, kwargs):
        assert len(kwargs) == 1, kwargs
        self.kwargs = kwargs

    def get_type(self):
        return self.format_string % self.kwargs

    @classmethod
    def from_parameters(cls, parameters):
        return cls(_convert_parameters(parameters, _substructure_converters))


class ChemFPSubstructOE_v1(_Substruct_v1):
    name = "ChemFPSubstruct-OE/1"
    num_bits = 881
    format_string = "ChemFP-OESubstruct/1 hydrogens=%(hydrogens)s"
    
    def read_structure_fingerprints(self, source=None, format=None):
        from chemfp.openeye_patterns import read_substruct_fingerprints_v1, SOFTWARE
        reader = read_substruct_fingerprints_v1(source, format)
        return self._open(SOFTWARE, source, reader)

ChemFPSubstructOE = ChemFPSubstructOE_v1

# This class name is just about incomprehensible
class ChemFPRDMACCSOE_v1(_Substruct_v1):
    name = "ChemFP-RDMaccs-OE/1"
    num_bits = 166
    format_string = "ChemFP-OESubstruct/1"
    
    def read_structure_fingerprints(self, source=None, format=None):
        from chemfp.openeye_patterns import read_rdmaccs_fingerprints_v1, SOFTWARE
        reader = read_rdmaccs_fingerprints_v1(source, format)
        return self._open(SOFTWARE, source, reader)

ChemFPRDMACCSOE = ChemFPRDMACCSOE_v1


_fingerprint_classes = [
    OpenEyeMACCS166_v1,
    OpenEyePath_v1,

    RDKitMACCS166_v1,
    RDKitFingerprint_v1,

    OpenBabelFP2_v1,
    OpenBabelFP3_v1,
    OpenBabelFP4_v1,
    OpenBabelMACCS166_v1,
    OpenBabelMACCS166_v2,

    ChemFPSubstructOE_v1,
    ChemFPRDMACCSOE_v1,
    ]
            
def parse_type(type):
    terms = type.split()
    if not terms:
        raise TypeError("missing name")
    
    name = terms[0]

    for cls in _fingerprint_classes:
        if cls.name == name:
            break
    else:
        raise TypeError

    seen = set()
    parameters = []
    for term in terms[1:]:
        left, right = term.split("=")
        if left in seen:
            raise TypeError("Duplicate name")
        seen.add(left)
        parameters.append((left, right))

    return cls.from_parameters(parameters)


def read_structure_fingerprints(type, source=None, format=None):
    structure_fingerprinter = parse_type(type)
    return structure_fingerprinter.read_structure_fingerprints(source, format)
