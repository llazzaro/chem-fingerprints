from __future__ import absolute_import

# WARNING! This is a first attempt at Indigo support. There are known problems
import warnings
warnings.warn("Indigo support is incomplete, experimental, and not to be trusted!")

from . import io
from . import types

from indigo import Indigo
_indigo = Indigo()

SOFTWARE = "Indigo/" + _indigo.version()

_expected_sizes = {
    "sim": 3736,
    "sub": 3736,
    "sub-res": 3736,
    "sub-tau": 3736,
    "full": 3736,
    }

def _initialize_sizes():
    mol = _indigo.loadMolecule("c1ccccc1O")
    for name, expected_num_bits in _expected_sizes.items():
        num_bits = len(mol.fingerprint(name).toString())*4
        assert num_bits == expected_num_bits, (name, expected_num_bits, num_bits)
_initialize_sizes()

class UnknownFormat(KeyError):
    def __str__(self):
        return "Unknown format %r" % (self.args[0],)

_formats = {
    "sdf": _indigo.iterateSDFile,
    "mol": _indigo.iterateSDFile,
    "sd": _indigo.iterateSDFile,
    "mdl": _indigo.iterateSDFile,

    "rdf": _indigo.iterateRDFile,


    "smi": _indigo.iterateSmilesFile,
    "can": _indigo.iterateSmilesFile,
    "smiles": _indigo.iterateSmilesFile,
    "ism": _indigo.iterateSmilesFile,
    
    "cml": _indigo.iterateCMLFile,
    }

def read_structures(filename=None, format=None):
    """read_structure(filename, format) -> (title, molecule) iterator

    Iterate over structures from filename, returning the structure
    title and Indigo molecule for each record. If filename is None
    this this reads from /dev/stdin instead of the named file.
    """
    if filename is None:
        if os.path.exists("/dev/stdin"):
            filename = "/dev/stdin"
        else:
            raise TypeError("filename of None (stdin) not supported on this platform")
    elif not isinstance(filename, basestring):
        raise TypeError("'filename' must be None or a string")

    format_name, compression = io.normalize_format(filename, format,
                                                   default=("smi", ""))
    if compression not in ("", ".gz"):
        raise TypeError("Unsupported compression type for %r" % (filename,))

    # Indigo automatically detects gzip compression
    reader_factory = _formats.get(format_name, None)
    if reader_factory is None:
        raise UnknownFormat(format)

    reader = reader_factory(filename)
    def indigo_structure_reader(reader):
        for mol in reader:
            yield mol.name(), mol

    return indigo_structure_reader(reader)

class _IndigoFingerprinter(types.Fingerprinter):
    def __init__(self, kwargs):
        self.num_bits = _expected_sizes[self.indigo_name]
        super(_IndigoFingerprinter, self).__init__(kwargs)

    def fingerprint(self, mol):
        return mol.fingerprint(self.indigo_name)

    def _get_reader(self, source, format, kwargs):
        indigo_name = self.indigo_name
        assert not kwargs, "need to support kwargs here"

        reader = read_structures(source, format)
        
        def read_indigo_fingerprints(indigo_name, reader):
            # XXX does read_structures need to take the Indigo
            # instance with the correct fingerprint sizes?
            for (title, mol) in reader:
                # toBuffer returns an array.array but I want a byte string
                yield mol.fingerprint(indigo_name).toBuffer().tostring(), title
        return read_indigo_fingerprints(indigo_name, reader)

# TODO: These can all be parameterized. Each one needs its own Indigo()
#   with the appropriate setOption()
# TODO: Get the correct subset of the fingerprint

class IndigoSimilarityFingerprinter_v1(_IndigoFingerprinter):
    name = "Indigo-Similarity/1"
    software = SOFTWARE
    indigo_name = "sim"
                              
class IndigoSubstructureFingerprinter_v1(_IndigoFingerprinter):
    name = "Indigo-Substructure/1"
    software = SOFTWARE
    indigo_name = "sub"
                              
class IndigoResonanceSubstructureFingerprinter_v1(_IndigoFingerprinter):
    name = "Indigo-ResonanceSubstructure/1"
    software = SOFTWARE
    indigo_name = "sub-res"
                              
class IndigoTautomerSubstructureFingerprinter_v1(_IndigoFingerprinter):
    name = "Indigo-TautomerSubstructure/1"
    software = SOFTWARE
    indigo_name = "sub-tau"
                              
class IndigoFullFingerprinter_v1(_IndigoFingerprinter):
    name = "Indigo-Full/1"
    software = SOFTWARE
    indigo_name = "full"
                              
