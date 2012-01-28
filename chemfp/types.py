from __future__ import absolute_import
# Information about fingerprint types
from . import argparse
from . import FingerprintIterator, Metadata

from . import io
from .decoders import import_decoder  # XXX too specific to the decoder module


def check_openbabel_maccs166():
    from .openbabel import HAS_MACCS, MACCS_VERSION
    assert HAS_MACCS
    if MACCS_VERSION == 1:
        return "OpenBabel-MACCS/1"
    elif MACCS_VERSION == 2:
        return "OpenBabel-MACCS/2"
    raise AssertionError

class FamilyProxy(object):
    def __init__(self, family_name, path):
        self.family_name = family_name
        self.path = path

    def __call__(self, **kwargs):
        cls = import_decoder(self.path)
        assert cls.name == self.family_name
        return cls(kwargs)

    def make_fingerprinter(self, kwargs):
        family = import_decoder(self.path)
        assert family.name == self.family_name
        return family.make_fingerprinter(kwargs)

### The chemfp fingerprint type API isn't powerful enough

# I have to list all of the possible fingerprint types, even if the
# platform doesn't support that specific type. There's also no support
# for toolkit vendors which do carefully tracks the fingerprint
# version (like OEChem); I should be using their version information
# rather than doing it myself.

_family_config_paths = (
    ("OpenEye-MACCS166/1", "chemfp.openeye.OpenEyeMACCSFingerprintFamily_v1"),
    ("OpenEye-Path/1", "chemfp.openeye.OpenEyePathFingerprintFamily_v1"),

    ("OpenEye-MACCS166/2", "chemfp.openeye.OpenEyeMACCSFingerprintFamily_v2"),
    ("OpenEye-Path/2", "chemfp.openeye.OpenEyePathFingerprintFamily_v2"),
    ("OpenEye-Circular/2", "chemfp.openeye.OpenEyeCircularFingerprintFamily_v2"),
    ("OpenEye-Tree/2", "chemfp.openeye.OpenEyeTreeFingerprintFamily_v2"),
    
    ("RDKit-MACCS166/1", "chemfp.rdkit.RDKitMACCSFingerprintFamily_v1"),
    ("RDKit-Fingerprint/1", "chemfp.rdkit.RDKitFingerprintFamily_v1"),
    ("RDKit-Morgan/1", "chemfp.rdkit.RDKitMorganFingerprintFamily_v1"),
    ("RDKit-Torsion/1", "chemfp.rdkit.RDKitTorsionFingerprintFamily_v1"),
    ("RDKit-Pair/1", "chemfp.rdkit.RDKitPairFingerprintFamily_v1"),
    
    ("OpenBabel-FP2/1", "chemfp.openbabel.OpenBabelFP2FingerprintFamily_v1"),
    ("OpenBabel-FP3/1", "chemfp.openbabel.OpenBabelFP3FingerprintFamily_v1"),
    ("OpenBabel-FP4/1", "chemfp.openbabel.OpenBabelFP4FingerprintFamily_v1"),
    ("OpenBabel-MACCS/1", "chemfp.openbabel.OpenBabelMACCSFingerprintFamily_v1"),
    ("OpenBabel-MACCS/2", "chemfp.openbabel.OpenBabelMACCSFingerprintFamily_v2"),

    ("Indigo-Similarity/1", "chemfp.indigo.IndigoSimilarityFingerprinter_v1"),
    ("Indigo-Substructure/1",
                      "chemfp.indigo.IndigoSubstructureFingerprinter_v1"),
    ("Indigo-ResonanceSubstructure/1",
                      "chemfp.indigo.IndigoResonanceSubstructureFingerprinter_v1"),
    ("Indigo-TautomerSubstructure/1",
                      "chemfp.indigo.IndigoTautomerSubstructureFingerprinter_v1"),
    ("Indigo-Full/1", "chemfp.indigo.IndigoFullFingerprinter_v1"),

    # In the future this will likely change to use a parameterized class
    # which can dynamically load fingerprint definitions

    ("ChemFP-Substruct-OpenEye/1",
                      "chemfp.openeye_patterns.SubstructOpenEyeFingerprinter_v1"),
    ("RDMACCS-OpenEye/1",
                      "chemfp.openeye_patterns.RDMACCSOpenEyeFingerprinter_v1"),

    ("ChemFP-Substruct-RDKit/1",
                      "chemfp.rdkit_patterns.SubstructRDKitFingerprinter_v1"),
    ("RDMACCS-RDKit/1",
                      "chemfp.rdkit_patterns.RDMACCSRDKitFingerprinter_v1"),

    ("ChemFP-Substruct-OpenBabel/1",
                      "chemfp.openbabel_patterns.SubstructOpenBabelFingerprinter_v1"),
    ("RDMACCS-OpenBabel/1",
                      "chemfp.openbabel_patterns.RDMACCSOpenBabelFingerprinter_v1"),

    ("ChemFP-Substruct-Indigo/1",
                      "chemfp.indigo_patterns.SubstructIndigoFingerprinter_v1"),
    ("RDMACCS-Indigo/1",
                      "chemfp.indigo_patterns.RDMACCSIndigoFingerprinter_v1"),
)

_alternates = {
    "OpenBabel-MACCS": check_openbabel_maccs166
    }


def _initialize_families(config_paths):
    d = {}
    for name, path in config_paths:
        # Set both the versioned and non-versioned names
        # The paths must be in order from oldest to newest
        
        unversioned_name = name.split("/")[0]
        d[unversioned_name] = d[name] = path
    return d

# Convert into a dictionary, and include the unversioned named
_family_config_paths = _initialize_families(_family_config_paths)

_loaded_families = {}

# Return the fingerprint family given a name like
# "OpenBabel-FP2" or "RDKit-Morgan/1"

def get_fingerprint_family(name):
    try:
        return _loaded_families[name]
    except KeyError:
        pass
    
    # Let's see if we can load it.

    # Is there a better name for this?
    try:
        func = _alternates[name]
    except KeyError:
        new_name = name
    else:
        new_name = func()

    try:
        path = _family_config_paths[new_name]
    except KeyError:
        raise ValueError("Unknown fingerprint family %r" % (name,))
    
    config = import_decoder(path)
    config.validate()
    
    family = FingerprintFamily(config)
    _loaded_families[name] = _loaded_families[new_name] = family
    return family

class FingerprintFamily(object):
    def __init__(self, config):
        self.config = config
        #name = config.name
        #format_string = config.format_string

    def __call__(self, **kwargs):
        return Fingerprinter(self.config, kwargs)

    def make_fingerprinter_from_type(self, type):
        terms = type.split()
        if not terms:
            raise ValueError("missing name in fingerprint type %r" % (type,))

        required_args = self.config.get_args()

        kwargs = {}
        for term in terms[1:]:
            try:
                left, right = term.split("=")
            except ValueError:
                raise ValueError("Term %r in type %r must have one and only one '='" %
                                 (term, type))
            if left in kwargs:
                raise ValueError("Duplicate name %r in type %r" % (left, type))

            try:
                decoder = required_args[left].decoder
            except KeyError:
                raise ValueError("Unknown name %r in type %r" % (left, type))
            
            try:
                value = decoder(right)
            except ValueError, err:
                raise ValueError("Unable to parse %s value %r in type %r" % (
                    left, right, type))
            
            kwargs[left] = value

        # Fill in any missing default
        for name, arg in required_args.items():
            if name not in kwargs:
                kwargs[name] = arg.default

        # Let the configuration verify the kwargs
        verify_args = self.config.verify_args
        if verify_args is not None:
            verify_args(kwargs)

        return Fingerprinter(self.config, kwargs)


    def __repr__(self):
        1/0

    def __eq__(self):
        2/0

    def __ne__(self):
        3/0

    
    

class Fingerprinter(object):
    def __init__(self, config, fingerprinter_kwargs):
        self.config = config
        if isinstance(config.num_bits, int):
            self.num_bits = config.num_bits
        elif config.num_bits is None:
            raise AssertionError(config.name)
        else:
            self.num_bits = config.num_bits(fingerprinter_kwargs)
        
        self.fingerprinter_kwargs = fingerprinter_kwargs

    def __eq__(self, other):
        return self.get_type() == other.get_type()

    def __ne__(self, other):
        return self.get_type() != other.get_type()

    def get_type(self):
        if self.config.format_string is None:
            assert not self.fingerprinter_kwargs, "kwargs but no format string!"
            return self.config.name
        encoded = self.config.format_string % self._encode_parameters()
        return self.config.name + " " + encoded

    def _encode_parameters(self):
        d = {}
        for k, v in self.fingerprinter_kwargs.items():
            encoder = self.config.args[k].encoder
            if encoder is None:
                d[k] = v
            else:
                d[k] = encoder(v)
        return d

    def read_structure_fingerprints(self, source, format=None, id_tag=None,
                                    errors="strict", metadata=None):
        source_filename = io.get_filename(source)
        if source_filename is None:
            sources = []
        else:
            sources = [source_filename]
            
        if metadata is None:
            # XXX I don't like how the user who wants to pass in aromaticity
            # information needs to create the full Metadata
            metadata = Metadata(num_bits=self.num_bits, type=self.get_type(),
                                software=self.config.software,
                                sources=sources)
            
        structure_reader = self.config.read_structures(metadata, source, format, id_tag, errors)
        fingerprinter = self.config.make_fingerprinter(**self.fingerprinter_kwargs)

        def fingerprint_reader(structure_reader, fingerprinter):
            for (id, mol) in structure_reader:
                yield id, fingerprinter(mol)
        reader = fingerprint_reader(structure_reader, fingerprinter)
        
        return FingerprintIterator(Metadata(num_bits = self.num_bits,
                                            sources = sources,
                                            software = self.config.software,
                                            type = self.get_type(),
                                            date = io.utcnow(),
                                            aromaticity = metadata.aromaticity),
                                   reader)
    
    def describe(self, bitno):
        if not (0 <= bitno < self.num_bits):
            raise KeyError("Bit number %d is out of range" % (bitno,))

        bit_description = self.config.bit_description
        if bit_description is None:
            return "(unknown)"
        return bit_description[bitno]

class Dummy(object):
    def __str__(self):
        return ""
    def __int__(self):
        return 0

class GetArgs(object):
    def __init__(self):
        self.args = []
    def __getitem__(self, name):
        # O(n**2) but n is small.
        if name in self.args:
            raise TypeError
        self.args.append(name)
        return Dummy()

def _get_arg_names(s):
    if not s:
        return []
    args = GetArgs()
    s % args
    return args.args
    
class FingerprintArgument(object):
    def __init__(self, name, decoder, encoder, kwargs):
        self.name = name
        self.decoder = decoder
        self.encoder = encoder
        self.default = kwargs["default"]
        self.kwargs = kwargs

# Can't use "first or second" because some of the 'first' arguments
# can be false value like {}.
def OR(first, second):
    if first is None:
        return second
    return first
    
class FingerprintFamilyConfig(object):
    def __init__(self,
                 name = None,
                 format_string = None,
                 software = None,
                 num_bits = None,
                 read_structures = None,
                 make_fingerprinter = None,
                 verify_args = None,
                 args = None,
                 ):
        self.name = name
        self.format_string = format_string
        self.software = software
        self.num_bits = num_bits
        self.read_structures = read_structures
        self.make_fingerprinter = make_fingerprinter
        if args is None:
            args = {}
        self.verify_args = verify_args
        self.args = args.copy() # This can contain extra args!
        self._exact_args = None

    def validate(self):
        pass

    def get_args(self):
        args = self._exact_args
        if args is None:
            result = {}
            for name in _get_arg_names(self.format_string):
                result[name] = self.args[name]
            args = self._exact_args = result
        return args

    def clone(self, name=None, format_string=None, software=None,
              num_bits=None, read_structures=None, make_fingerprinter=None,
              verify_args=None, args=None):
        return FingerprintFamilyConfig(
            name = OR(name, self.name),
            format_string = OR(format_string, self.format_string),
            software = OR(software, self.software),
            num_bits = OR(num_bits, self.num_bits),
            read_structures = OR(read_structures, self.read_structures),
            make_fingerprinter = OR(make_fingerprinter, self.make_fingerprinter),
            verify_args = OR(verify_args, self.verify_args),
            args = OR(args, self.args))

    def add_argument(self, name, decoder=None, encoder=None, default=None,
                     action=None, metavar=None, help=None):
        #if name in self.args:
        #    raise AssertionError("Argument %r already added" % (name,))
        if default is not None:
            if help is not None:
                help = "%s (default=%s)" % (help, default)
        def parse_argument(s):
            try:
                return decoder(s)
            except ValueError, err:
                raise argparse.ArgumentError(None, "%s %s" % (name, err))
        arg = FingerprintArgument(name, decoder, encoder,
                                  kwargs = dict(type=parse_argument,
                                                default=default,
                                                action=action,
                                                metavar=metavar,
                                                help=help))
        self.args[name] = arg

    def remove_argument(self, name):
        del self.args[name]

    def add_argument_to_argparse(self, name, parser):
        info = self.args[name]
        parser.add_argument("--" + info.name, **info.kwargs)
        

# Helper functions

def positive_int(s):
    # Don't do int(s) because that allows "+3" and " 3 ", which I don't want
    if not s.isdigit():
        raise ValueError("must be 1 or greater")
    i = int(s)
    if i == 0:
        raise ValueError("must be 1 or greater")
    return i

def nonnegative_int(s):
    if not s.isdigit():
        raise ValueError("must be 0 or greater")
    return int(s)

def zero_or_one(s):
    if s == "0":
        return 0
    if s == "1":
        return 1
    raise ValueError("must be 0 or a 1")
        

def parse_type(type):
    terms = type.split()
    if not terms:
        raise ValueError("missing name in fingerprint type %r" % (type,))

    name = terms[0]
    family = get_fingerprint_family(name)
    return family.make_fingerprinter_from_type(type)
