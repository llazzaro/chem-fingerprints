# Information about fingerprint types

class Header(object):
    pass

from datetime import datetime
def to_header(**kwargs):
    header = Header()
    for name in ("num_bits", "software", "type", "source"):
        value = kwargs.pop(name)
        if value is not None:
            setattr(header, name, value)
    date = kwargs.pop("date")
    if date is None:
        date = datetime.utcnow().isoformat().split(".", 1)[0]
    setattr(header, "date", date)
    return header

class StructureFPReader(object):
    def __init__(self, header, fingerprint_reader):
        self.header = header
        self.fingerprint_reader = fingerprint_reader

    def iter_fingerprints(self):
        return self.fingerprint_reader

    def __iter__(self):
        return self.fingerprint_reader

def _convert_parameters(parameters, converters):
    kwargs = {}
    for (name, value) in parameters:
        to_name, convert = converters[name]
        kwargs[to_name] = convert(value)
    return kwargs

class _Opener(object):
    def _open(self, software, source, reader):
        if source is not None:
            if not isinstance(source, basestring):
                # Then it's a Python file object
                source = getattr(source, "name", None)
        
        return StructureFPReader(to_header(num_bits = self.num_bits,
                                           source = source,
                                           software = software,
                                           type = self.get_type(),
                                           date = None),
                                 reader)
class _NoParameters(_Opener):
    @classmethod
    def from_parameters(cls, parameters):
        assert len(parameters) == 0
        return cls

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
    format_string = "OpenEye-Path/1 numbits={num_bits} minbonds={min_bonds} maxbonds={max_bonds} atype={atype} btype={btype}"
    
    converters = {"numbits": ("num_bits", int),
                  "minbonds": ("min_bonds", int),
                  "maxbonds": ("max_bonds", int)}

    def __init__(self, kwargs):
        assert len(kwargs) == 5, kwargs
        self.kwargs = kwargs
        self.num_bits = kwargs["num_bits"]
    
    @staticmethod
    def from_parameters(parameters):
        converters = OpenEyePath_v1.converters
        if len(converters) == 3:
            from chemfp.openeye import atom_description_to_value, bond_description_to_value
            converters["atype"] = ("atype", atom_description_to_value)
            converters["btype"] = ("btype", bond_description_to_value)
        
        return OpenEyePath_v1(_convert_parameters(parameters, converters))
    
        
    def get_type(self):
        from chemfp.openeye import atom_value_to_description, bond_value_to_description
        kw = self.kwargs
        return self.format_string.format(num_bits = kw["num_bits"],
                                         min_bonds = kw["min_bonds"],
                                         max_bonds = kw["max_bonds"],
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
    format_string = ("RDKit-Fingerprint/1 minPath={min_path} maxPath={max_path} fpSize={num_bits} "
                     "nBitsPerHash={bits_per_hash} useHs={use_Hs}")

    converters = {"minPath": ("min_path", int),
                  "maxPath": ("max_path", int),
                  "fpSize": ("num_bits", int),
                  "nBitsPerHash": ("bits_per_hash", int),
                  "useHs": ("use_Hs", int)}
    def __init__(self, kwargs):
        assert len(kwargs) == 5
        self.num_bits = kwargs["num_bits"]
        self.kwargs = kwargs

    @staticmethod
    def from_parameters(parameters):
        return RDKitFingerprint_v1(_convert_parameters(parameters, RDKitFingerprint_v1.converters))

    def get_type(self):
        return self.format_string.format(**self.kwargs)

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

    print "Call with", parameters
    return cls.from_parameters(parameters)


def read_structure_fingerprints(typeinfo, source=None, format=None):
    structure_fingerprinter = parse_type(typeinfo)
    return structure_fingerprinter.read_structure_fingerprints(source, format)
