# Information about fingerprint types

from . import io
from .decoders import import_decoder  # XXX too specific to the decoder module


def check_openbabel_maccs166():
    from chemfp.openbabel import HAS_MACCS, MACCS_VERSION
    assert HAS_MACCS
    if MACCS_VERSION == 1:
        return "OpenBabel-MACCS/1"
    elif MACCS_VERSION == 2:
        return "OpenBabel-MACCS/1"
    raise AssertionError


_fingerprints = [
    dict(name="OpenEye-MACCS166/1", 
         num_bits=166,
         reader="chemfp.openeye.read_maccs166_fingerprints_v1",
         software="chemfp.openeye.SOFTWARE",
         ),

    dict(name="OpenEye-Path/1",
         reader="chemfp.openeye.read_path_fingerprints_v1",
         software="chemfp.openeye.SOFTWARE",
         num_bits=lambda kwargs: kwargs["numbits"],
         decode_parameters="chemfp.openeye.decode_maccs166_parameters",
         encode_parameters="chemfp.openeye.encode_maccs166_parameters",
         format_string=("numbits=%(numbits)s minbonds=%(minbonds)s "
                        "maxbonds=%(maxbonds)s atype=%(atype)s btype=%(btype)s")),

    dict(name="RDKit-MACCS166/1",
         num_bits=166, # XXX 166? 167?
         reader="chemfp.rdkit.read_maccs166_fingerprints_v1",
         software="chemfp.rdkit.SOFTWARE"),

    dict(name="RDKit-Fingerprint/1",
         reader="chemfp.rdkit.read_rdkit_fingerprints_v1",
         software="chemfp.rdkit.SOFTWARE",
         num_bits=lambda kwargs: kwargs["fpSize"],
         decode_parameters="chemfp.rdkit.decode_fingerprint_parameters",
         # Don't need encoders; the format_string is just fine
         format_string = (
             "RDKit-Fingerprint/1 minPath=%(minPath)s maxPath=%(maxPath)s fpSize=%(fpSize)s "
             "nBitsPerHash=%(nBitsPerHash)s useHs=%(useHs)s") ),


    dict(name="OpenBabel-FP2/1",
         num_bits = 1021,
         reader="chemfp.openbabel.read_fp2_fingerprints_v1",
         software="chemfp.openbabel.SOFTWARE"),

    dict(name="OpenBabel-FP3/1",
         num_bits = 55,
         reader="chemfp.openbabel.read_fp3_fingerprints_v1",
         software="chemfp.openbabel.SOFTWARE"),

    dict(name="OpenBabel-FP4/1",
         num_bits = 307,
         reader="chemfp.openbabel.read_fp4_fingerprints_v1",
         software="chemfp.openbabel.SOFTWARE"),

    dict(name="OpenBabel-MACCS/1",
         num_bits = 166,
         reader="chemfp.openbabel.read_maccs166_fingerprints_v1",
         software="chemfp.openbabel.SOFTWARE"),

    dict(name="OpenBabel-MACCS/2",
         num_bits = 166,
         reader="chemfp.openbabel.read_maccs166_fingerprints_v2",
         software="chemfp.openbabel.SOFTWARE"),


    ## pattern-based fingerprints

    dict(name="ChemFP-Substruct-OpenEye/1",
         num_bits=881,
         reader="chemfp.openeye_patterns.read_substruct_fingerprints_v1",
         software="chemfp.openeye_patterns.SOFTWARE"),

    dict(name="ChemFP-RDMACCS-OpenEye/1",
         num_bits=166,
         reader="chemfp.openeye_patterns.read_rdmaccs_fingerprints_v1",
         software="chemfp.openeye_patterns.SOFTWARE"),


    dict(name="ChemFP-Substruct-RDKit/1",
         num_bits=881,
         reader="chemfp.rdkit_patterns.read_substruct_fingerprints_v1",
         software="chemfp.openeye_patterns.SOFTWARE"),

    dict(name="ChemFP-RDMACCS-RDKit/1",
         num_bits=166,
         reader="chemfp.rdkit_patterns.read_rdmaccs_fingerprints_v1",
         software="chemfp.openeye_patterns.SOFTWARE"),


    dict(name="ChemFP-Substruct-OpenBabel/1",
         num_bits=881,
         reader="chemfp.openbabel_patterns.read_substruct_fingerprints_v1",
         software="chemfp.openbabel_patterns.SOFTWARE"),

    dict(name="ChemFP-RDMACCS-OpenBabel/1",
         num_bits=166,
         reader="chemfp.openbabel_patterns.read_rdmaccs_fingerprints_v1",
         software="chemfp.openbabel_patterns.SOFTWARE"),
    

]

_alternates = {
    "OpenBabel-MACCS": check_openbabel_maccs166
    }



_load_fields = ("reader", "decode_parameters", "encode_parameters", "software")
class ConfigLoader(object):
    def __init__(self, config):
        self.config = config
        self.loaded = {}
        for k,v in config.items():
            if k not in _load_fields:
                self.loaded[k] = v

    def __getitem__(self, name):
        try:
            return self.loaded[name]
        except KeyError:
            if name not in _load_fields:
                raise

        path = self.config[name]
        obj = import_decoder(path)
        self.loaded[name] = obj
        return obj

    def get(self, name, default=None):
        try:
            return self[name]
        except KeyError:
            return default

class FingerprintFamily(object):
    def __init__(self, config):
        self.config = ConfigLoader(config)

    def _decode_parameters(self, parameters):
        decode_parameters = self.config.get("decode_parameters", None)
        if decode_parameters is None:
            assert not parameters
            return {}
        if isinstance(decode_parameters, basestring):
            decode_parameters = import_decoder(decode_parameters)
            # Replace for future use
            self.config["decode_parameters"] = decode_parameters
        return decode_parameters(parameters)

    def __call__(self, **kwargs):
        return self.from_kwargs(kwargs)

    def from_parameters(self, parameters):
        kwargs = self._decode_parameters(parameters)
        return self.from_kwargs(kwargs)

    def from_kwargs(self, kwargs=None):
        if kwargs is None:
            kwargs = {}
        num_bits = self.config["num_bits"]
        if not isinstance(num_bits, int):
            num_bits = num_bits(kwargs)
        return FingerprintType(num_bits, self.config, kwargs)

    
_family_by_name = {}

def _initialize_families():
    for config in _fingerprints:
        # Set both the versioned and non-versioned names
        name = config["name"]
        unversioned_name = name.split("/")[0]
        _family_by_name[name] = _family_by_name[unversioned_name] = FingerprintFamily(config)

    # Don't include a (likely non-versioned) name if there's a selector function
    for name in alternates:
        if name in _family_by_name:
            del _family_by_name[name]

_initialize_families()


def get_fingerprint_family(name):
    try:
        return _family_by_name[name]
    except KeyError:
        if name not in alternates:
            raise
    alternate = _alternates[name]()
    return _family_by_name[alternate]


class FingerprintType(object):
    def __init__(self, num_bits, config, kwargs):
        self.num_bits = num_bits
        self.config = config
        self.kwargs = kwargs

    def _encode_parameters(self, kwargs):
        print "Want", self.config
        encode_parameters = self.config.get("encode_parameters", None)
        if encode_parameters is None:
            assert not kwargs, kwargs
            return {}
        return encode_parameters(kwargs)
    
    def get_type(self):
        format_string = self.config.get("format_string", None)
        if format_string is None:
            assert not self.kwargs
            return self.config["name"]
        encoded = format_string % self._encode_parameters(self.kwargs)
        return self.config["name"] + " " + encoded

    def read_structure_fingerprints(self, source, format=None):
        reader = self.config["reader"](source, format, self.kwargs)

        source_filename = io.get_filename(source)
        return io.FPIterator(io.Header(num_bits = self.num_bits,
                                       source = source_filename,
                                       type = self.get_type(),
                                       date = io.utcnow()),
                             reader)

            
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
