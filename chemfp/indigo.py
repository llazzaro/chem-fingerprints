from __future__ import absolute_import

# WARNING! This is a first attempt at Indigo support. There are known problems
import warnings
warnings.warn("Indigo support is incomplete, experimental, and not to be trusted!")

from . import io
from . import types

from indigo import Indigo
_indigo = Indigo()

SOFTWARE = "Indigo/" + _indigo.version()

# There is no way to get the sizes from the SDK API.  While it is
# possible to change the fingerprint sizes, I decided to not support
# that for now since it appears that things will change.

# All fingerprints are currently stored in a single mega-fingerprint
# object, broken down into parts:
#   3   bytes for "extra" bits
#  25*8 bytes for "ordinary" part
#   8*8 bytes for "similarity" part
#  10*8 bytes for "tautomer" part
#  15*8 bytes for "resonance" part
_ORD = 3
_SIM = 3+(25)*8
_TAU = 3+(25+8)*8
_RES = 3+(25+8+10)*8
_END = 3+(25+8+10+15)*8

_expected_sizes = {
    "sim": (_TAU-_SIM)*8,
    "sub": ((_SIM-_ORD)+(_END-_RES))*8,
    "sub-res": (_END-_RES)*8,
    "sub-tau": (_RES-_TAU)*8,
    "full": (_END-_ORD)*8,
    }

def _verify():
    mol = _indigo.loadMolecule("c1ccccc1OSNOP")
    fp = mol.fingerprint("sim")
    
    for name, expected_num_bits in _expected_sizes.items():
        num_bits = len(mol.fingerprint(name).toString())*4
        assert num_bits == expected_num_bits, (name, expected_num_bits, num_bits)


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
            mol.aromatize()
            yield mol.name(), mol

    return indigo_structure_reader(reader)

class _IndigoFingerprinter(types.Fingerprinter):
    def __init__(self, kwargs):
        self.num_bits = _expected_sizes[self.indigo_name]
        super(_IndigoFingerprinter, self).__init__(kwargs)

    def fingerprint(self, mol):
        return mol.fingerprint(self.indigo_name)

    def _get_reader(self, source, format, kwargs):
        assert not kwargs, "need to support kwargs here"
        reader = read_structures(source, format)
        #for fp, title in self._read_indigo_fingerprints(self.indigo_name, reader):
        #    assert len(fp)*8 == self.num_bits, (len(fp)*8, self.num_bits)
        #    yield fp, title
        return  self._read_indigo_fingerprints(self.indigo_name, reader)

# TODO: These can all be parameterized. Each one needs its own Indigo()
#   with the appropriate setOption()
# TODO: Get the correct subset of the fingerprint

class IndigoSimilarityFingerprinter_v1(_IndigoFingerprinter):
    name = "Indigo-Similarity/1"
    software = SOFTWARE
    indigo_name = "sim"

    def _read_indigo_fingerprints(self, indigo_name, reader):
        START = _SIM
        END = _TAU
        for (title, mol) in reader:
            fp = mol.fingerprint(indigo_name).toBuffer().tostring()
            #assert len(set(fp[3:START])) == 1, fp.encode("hex")
            #assert len(set(fp[END:])) == 1, fp.encode("hex")
            yield fp[START:END], title
                              
class IndigoSubstructureFingerprinter_v1(_IndigoFingerprinter):
    name = "Indigo-Substructure/1"
    software = SOFTWARE
    indigo_name = "sub"

    def _read_indigo_fingerprints(self, indigo_name, reader):
        START1 = _ORD
        END1 = _SIM
        START2 = _RES
        END2 = _END
        for (title, mol) in reader:
            fp = mol.fingerprint(indigo_name).toBuffer().tostring()
            #assert len(set(fp[END1:START2])) == 1, fp.encode("hex")
            yield fp[START1:END1]+fp[START2:END2], title

                              
class IndigoResonanceSubstructureFingerprinter_v1(_IndigoFingerprinter):
    name = "Indigo-ResonanceSubstructure/1"
    software = SOFTWARE
    indigo_name = "sub-res"

    def _read_indigo_fingerprints(self, indigo_name, reader):
        START = _RES
        END = _END
        for (title, mol) in reader:
            fp = mol.fingerprint(indigo_name).toBuffer().tostring()
            #assert len(set(fp[3:START])) == 1, fp.encode("hex")
            yield fp[START:END], title

                              
class IndigoTautomerSubstructureFingerprinter_v1(_IndigoFingerprinter):
    name = "Indigo-TautomerSubstructure/1"
    software = SOFTWARE
    indigo_name = "sub-tau"

    def _read_indigo_fingerprints(self, indigo_name, reader):
        START = _TAU
        END = _RES
        for (title, mol) in reader:
            fp = mol.fingerprint(indigo_name).toBuffer().tostring()
            #assert len(set(fp[3:START])) == 1, fp.encode("hex")
            #assert len(set(fp[END:])) == 1, fp.encode("hex")
            yield fp[START:END], title
                              
class IndigoFullFingerprinter_v1(_IndigoFingerprinter):
    name = "Indigo-Full/1"
    software = SOFTWARE
    indigo_name = "full"
                              
    def _read_indigo_fingerprints(self, indigo_name, reader):
        START = _ORD
        END = _END
        for (title, mol) in reader:
            fp = mol.fingerprint(indigo_name).toBuffer().tostring()
            yield fp[3:], title
