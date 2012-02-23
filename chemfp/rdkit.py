"Create RDKit fingerprints"

# Copyright (c) 2010 Andrew Dalke Scientific, AB (Gothenburg, Sweden)
# See the contents of "__init__.py" for full license details.

from __future__ import absolute_import

import os
import sys
import gzip

import rdkit
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import rdkit.rdBase
from rdkit.Chem.MACCSkeys import GenMACCSKeys

from . import sdf_reader
from .decoders import from_binary_lsb as _from_binary_lsb
from . import io
from . import types

# These are the things I consider to be public
__all__ = ["read_structures", "iter_smiles_molecules", "iter_sdf_molecules"]


# If the attribute doesn't exist then this is an unsupported pre-2010 RDKit distribution
SOFTWARE = "RDKit/" + getattr(rdkit.rdBase, "rdkitVersion", "unknown")

# Used to check for version-dependent fingerprints
_VERSION_PROBE_MOL = Chem.MolFromSmiles(r"CC1=CC(=NN1CC(=O)NNC(=O)\C=C\C2=C(C=CC=C2Cl)F)C")

#########

# Helper function to convert a fingerprint to a sequence of bytes.

from rdkit import DataStructs
if getattr(DataStructs, "BitVectToBinaryText", None):
    _fp_to_bytes = DataStructs.BitVectToBinaryText
else:
    # Support for pre-2012 releases of RDKit
    def _fp_to_bytes(fp):
        return _from_binary_lsb(fp.ToBitString())[1]

#########
_allowed_formats = ["sdf", "smi"]
_format_extensions = {
    ".sdf": "sdf",
    ".mol": "sdf",
    ".sd": "sdf",
    ".mdl": "sdf",

    ".smi": "smi",
    ".can": "smi",
    ".smiles": "smi",
    ".ism": "smi",
}


class SmilesFileLocation(object):
    def __init__(self, name=None):
        self.name = name
        self.lineno = 1
    def where(self):
        s = "at line %(lineno)s"
        if self.name is not None:
            s += " of %(name)s"
        return s % self.__dict__
    

# While RDKit has a SMILES file parser, it doesn't handle reading from
# stdin or from compressed files. I wanted to support those as well, so
# ended up not using Chem.SmilesMolSupplier.

def iter_smiles_molecules(fileobj, name=None, errors="strict"):
    """Iterate over the SMILES file records, returning (title, RDKit.Chem.Mol) pairs

    'fileobj' is an input file or any line iterable
    'name' is the name used to report errors (if not specified, use
       fileobj.name if present)
    'errors' is one of "strict" (default), "log", or "ignore" (other values are experimental)

    Each line of the input must at least one whitespace separated
    fields.  The first field is the SMILES. If there is a second field
    then it is used as the title, otherwise the title is the current
    record number, starting with "1".

    """
    if name is None:
        name = getattr(fileobj, "name", None)
    error_handler = sdf_reader.get_parse_error_handler(errors)

    loc = SmilesFileLocation(name)
    for lineno, line in enumerate(fileobj):
        words = line.split()
        if len(words) <= 1:
            loc.lineno = lineno+1
            if not words:
                error_handler("Unexpected blank line", loc)
            else:
                error_handler("Missing SMILES name (second column)", loc)
                continue

        mol = Chem.MolFromSmiles(words[0])
        if mol is None:
            loc.lineno = lineno+1
            error_handler("Cannot parse the SMILES %r" % (words[0],), loc)
            continue
        
        yield words[1], mol


def iter_sdf_molecules(fileobj, name=None, id_tag=None, errors="strict"):
    """Iterate over the SD file records, returning (id, Chem.Mol) pairs

    fileobj - the input file object
    name - the name to use to report errors. If None, use fileobj.name .
    
    """
    # If there's no explicit filename, see if fileobj has one
    if name is None:
        name = getattr(fileobj, "name", None)
    loc = sdf_reader.FileLocation(name)
    error = sdf_reader.get_parse_error_handler(errors)
    if id_tag is None:
        for i, text in enumerate(sdf_reader.iter_sdf_records(fileobj, errors, loc)):
            mol = Chem.MolFromMolBlock(text)
            if mol is None:
                # This was not a molecule?
                error("Could not parse molecule block", loc)
                continue
            title = mol.GetProp("_Name")
            id = io.remove_special_characters_from_id(title)
            if not id:
                error("Missing title for record #%d" % (i+1), loc)
                continue
            yield id, mol
    else:
        # According to
        #   http://www.mail-archive.com/rdkit-discuss@lists.sourceforge.net/msg01436.html
        # I can make a new SDMolSupplier, then SetData(), get the first record, and
        # get its property names. That's ... crazy.
        sdf_iter = sdf_reader.iter_sdf_records(fileobj, errors, loc)
        for i, (id, text) in enumerate(sdf_reader.iter_tag_and_record(sdf_iter, id_tag)):
            mol = Chem.MolFromMolBlock(text)
            if mol is None:
                # This was not a molecule?
                error("Could not parse molecule block", loc)
                continue
            if id is None:
                error("Missing id tag %r for record #%d" % (id_tag, i+1), loc)
                continue
            id = io.remove_special_characters_from_id(id)
            if not id:
                error("Empty id tag %r for record #%d" % (id_tag, i+1), loc)
                continue
            yield id, mol
            
# this class helps the case when someone is entering structure
# by-hand. (Most likely to occur with SMILES input). They would like
# to see the result as soon as a record is entered. But normal
# interation reader grabs a buffer of input to process, and not a
# line. It's faster that way. The following adapter supports the
# iterator protocol but turns it into simple readlines(). This will be
# slower but since do it only if stdin is a tty, there shouldn't be a
# problem.
## class _IterUsingReadline(object):
##     "Internal class for iterating a line at a time from tty input"
##     def __init__(self, fileobj):
##         self.fileobj = fileobj
##     def __iter__(self):
##         return iter(self.fileobj.readline, "")

## def _open(filename, compressed):
##     "Internal function to open the given filename, which might be compressed"
##     if filename is None:
##         if compressed:
##             return gzip.GzipFile(fileobj=sys.stdin, mode="r")
##         else:
##             # Python's iter reads a block.
##             # When someone types interactively, read only a line.
##             if sys.stdin.isatty():
##                 return _IterUsingReadline(sys.stdin)
##             else:
##                 return sys.stdin

##     if compressed:
##         return gzip.GzipFile(filename, "r")
##     return open(filename, "rU")

def is_valid_format(format):
    if format is None:
        return True
    try:
        format_name, compression = io.normalize_format(None, format, ("smi", None))
    except ValueError:
        return False
    format_name = _format_extensions.get(format_name, format_name)
    return format_name in ("sdf", "smi")
    

def read_structures(source, format=None, id_tag=None, errors="strict"):
    """Iterate the records in the input source as (title, RDKit.Chem.Mol) pairs

    'source' is a filename, a file object, or None for stdin
    'format' is either "sdf" or "smi" with optional ".gz" or ".bz2" extensions.
        If None then the format is inferred from the source extension
    'errors' is one of "strict" (default), "log", or "ignore" (other values are experimental)
    """
    format_name, compression = io.normalize_format(source, format, default=("smi", None))
    format_name = _format_extensions.get(format_name, format_name)
    if format_name == "sdf":
        # I have an old PubChem file Compound_09425001_09450000.sdf .
        #   num. lines = 5,041,475   num. bytes = 159,404,037
        # 
        # Parse times for iter_sdf_records (parsing records in Python)
        #   37.6s (best of 37.6, 38.3, 37.8)
        # Parse times for the RDKit implementation (parsing records in C++)
        #   40.2s (best of 41.7, 41.33, 40.2)
        # 
        # The native RDKit reader is slower than the Python one and does
        # not have (that I can tell) support for compressed files, so
        # I'll go with the Python one. For those interested, here's the
        # RDKit version.
        # 
        #if (not compressed) and (source is not None):
        #    supplier = Chem.SDMolSupplier(source)
        #    def native_sdf_reader():
        #        for mol in supplier:
        #            if mol is None:
        #                print >>sys.stderr, "Missing? after", title
        #            else:
        #                title = mol.GetProp("_Name")
        #                yield title, mol
        #    return native_sdf_reader()

        fileobj = io.open_compressed_input_universal(source, compression)
        # fileobj should always have the .name attribute set.
        return iter_sdf_molecules(fileobj, None, id_tag, errors)

    elif format_name == "smi":
        # I timed the native reader at 31.6 seconds (best of 31.6, 31.7, 31.7)
        # and the Python reader at 30.8 seconds (best of 30.8, 30.9, and 31.0)
        # Yes, the Python reader is faster and using it gives me better consistency
        #
        #if (not compressed) and (source is not None):
        #    supplier = Chem.SmilesMolSupplier(source, delimiter=" \t", titleLine=False)
        #    def native_smiles_reader():
        #        for mol in supplier:
        #            yield mol.GetProp("_Name"), mol
        #    return native_smiles_reader()
        fileobj = io.open_compressed_input_universal(source, compression)
        return iter_smiles_molecules(fileobj, None, errors)

    else:
        if format is None:
            raise ValueError("Unknown structure filename extension: %r" % (source,))
        else:
            raise ValueError("Unknown structure format %r" % (format_name,))
        

########### The topological fingerprinter

# Some constants shared by the fingerprinter and the command-line code.

NUM_BITS = 2048
MIN_PATH = 1
MAX_PATH = 7
BITS_PER_HASH = 4
USE_HS = 1
assert USE_HS == 1, "Don't make this 0 unless you know what you are doing"

# Not supporting the tgtDensity and minSize options.
# This program generates fixed-length fingerprints.

def make_rdk_fingerprinter(minPath=MIN_PATH, maxPath=MAX_PATH, fpSize=NUM_BITS,
                           nBitsPerHash=BITS_PER_HASH, useHs=USE_HS):
    if not (fpSize > 0):
        raise ValueError("fpSize must be positive")
    if not (minPath > 0):
        raise ValueError("minPath must be positive")
    if not (maxPath >= minPath):
        raise ValueError("maxPath must not be smaller than minPath")
    if not (nBitsPerHash > 0):
        raise ValueError("nBitsPerHash must be positive")

    def rdk_fingerprinter(mol):
        fp = Chem.RDKFingerprint(
            mol, minPath=minPath, maxPath=maxPath, fpSize=fpSize,
            nBitsPerHash=nBitsPerHash, useHs=useHs)
        return _fp_to_bytes(fp)
    return rdk_fingerprinter

########### The MACCS fingerprinter


def maccs166_fingerprinter(mol):
    fp = GenMACCSKeys(mol)
    # In RDKit the first bit is always bit 1 .. bit 0 is empty (?!?!)
    bitstring_with_167_bits = fp.ToBitString()
    # I want the bits to start at 0, so I do a manual left shift
    return _from_binary_lsb(bitstring_with_167_bits[1:])[1]

def make_maccs166_fingerprinter():
    return maccs166_fingerprinter


########### The Morgan fingerprinter

# Some constants shared by the fingerprinter and the command-line code.

RADIUS = 2
USE_FEATURES = 0
USE_CHIRALITY = 0
USE_BOND_TYPES = 1

def make_morgan_fingerprinter(fpSize=NUM_BITS,
                              radius=RADIUS,
                              useFeatures=USE_FEATURES,
                              useChirality=USE_CHIRALITY,
                              useBondTypes=USE_BOND_TYPES):
    if not (fpSize > 0):
        raise ValueError("fpSize must be positive")
    if not (radius >= 0):
        raise ValueError("radius must be positive or zero")

    def morgan_fingerprinter(mol):
        fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(
            mol, radius, nBits=fpSize, useChirality=useChirality,
            useBondTypes=useBondTypes,useFeatures=useFeatures)
        return _fp_to_bytes(fp)
    return morgan_fingerprinter


########### Torsion fingerprinter

TARGET_SIZE = 4

def make_torsion_fingerprinter(fpSize=NUM_BITS,
                               targetSize=TARGET_SIZE):
    if not (fpSize > 0):
        raise ValueError("fpSize must be positive")
    if not (targetSize >= 0):
        raise ValueError("targetSize must be positive or zero")

    def torsion_fingerprinter(mol):
        fp = rdMolDescriptors.GetHashedTopologicalTorsionFingerprintAsBitVect(
            mol, nBits=fpSize, targetSize=targetSize)
        return _fp_to_bytes(fp)
    return torsion_fingerprinter

TORSION_VERSION = {
    "\xc2\x10@\x83\x010\x18\xa4,\x00\x80B\xc0\x00\x08\x00": "1",
    "\x13\x11\x103\x00\x007\x00\x00p\x01\x111\x0107": "2",
    }[make_torsion_fingerprinter(128)(_VERSION_PROBE_MOL)]

########### Atom Pair fingerprinter

MIN_LENGTH = 1
MAX_LENGTH = 30

def make_atom_pair_fingerprinter(fpSize=NUM_BITS,
                                 minLength=MIN_LENGTH,
                                 maxLength=MAX_LENGTH):
    if not (fpSize > 0):
        raise ValueError("fpSize must be positive")
    if not (minLength >= 0):
        raise ValueError("minLength must be positive or zero")
    if not (maxLength >= minLength):
        raise ValueError("maxLength must not be less than minLength")

    def pair_fingerprinter(mol):
        fp = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(
            mol, nBits=fpSize, minLength=minLength, maxLength=maxLength)
        return _fp_to_bytes(fp)
    return pair_fingerprinter

try:
    ATOM_PAIR_VERSION = {
        "\xfdB\xfe\xbd\xfa\xdd\xff\xf5\xff\x05\xdf?\xe3\xc3\xff\xfb": "1",
        "w\xf7\xff\xf7\xff\x17\x01\x7f\x7f\xff\xff\x7f\xff?\xff\xff": "2",
        }[make_atom_pair_fingerprinter(128)(_VERSION_PROBE_MOL)]
except Exception, err:
    # RDKit 2011.06 contained a bug
    if "Boost.Python.ArgumentError" in str(type(err)):
        ATOM_PAIR_VERSION = None
    else:
        raise


####################

from .types import FingerprintFamilyConfig, positive_int, nonnegative_int, zero_or_one

def _read_structures(metadata, source, format, id_tag, errors):
    if metadata.aromaticity is not None:
        raise ValueError("RDKit does not support alternate aromaticity models "
                         "(want aromaticity=%r)" % metadata.aromaticity)
    return read_structures(source, format, id_tag, errors)

# Check for metadata.aromaticity
_base = FingerprintFamilyConfig(
    software = SOFTWARE,
    read_structures = _read_structures,
    )

_base.add_argument("fpSize", decoder=positive_int, metavar="INT", default=NUM_BITS,
                   help = "number of bits in the fingerprint (applies to RDK, Morgan, topological torsion, and atom pair fingerprints")

_base.add_argument("minPath", decoder=positive_int, metavar="INT", default=MIN_PATH,
                   help = "minimum number of bonds to include in the subgraph")

_base.add_argument("maxPath", decoder=positive_int, metavar="INT", default=MAX_PATH,
                   help = "maximum number of bonds to include in the subgraph")

_base.add_argument("nBitsPerHash", decoder=positive_int, metavar="INT",
                   default=BITS_PER_HASH, help = "number of bits to set per path")

_base.add_argument("useHs", decoder=zero_or_one, metavar="0|1", default=USE_HS,
                   help = "include information about the number of hydrogens on each atom")              
# Morgan
_base.add_argument("radius", decoder=nonnegative_int, metavar="INT", default=RADIUS,
                   help = "radius for the Morgan algorithm")

_base.add_argument("useFeatures", decoder=zero_or_one, metavar="0|1",
                   default=USE_FEATURES, help = "use chemical-feature invariants")

_base.add_argument("useChirality", decoder=zero_or_one, metavar="0|1",
                   default=USE_CHIRALITY, help = "include chirality information")

_base.add_argument("useBondTypes", decoder=zero_or_one, metavar="0|1",
                   default=USE_BOND_TYPES, help = "include bond type information")


# torsion
_base.add_argument("targetSize", decoder=positive_int, metavar="INT",
                   default=TARGET_SIZE, help = "number of bits in the fingerprint")

# pair
_base.add_argument("minLength", decoder=nonnegative_int, metavar="INT",
                   default=MIN_LENGTH, help = "minimum bond count for a pair")

_base.add_argument("maxLength", decoder=nonnegative_int, metavar="INT",
                   default=MAX_LENGTH, help = "maximum bond count for a pair")

#########

RDKitMACCSFingerprintFamily_v1 = _base.clone(
    name = "RDKit-MACCS166/1",
    num_bits = 166,
    make_fingerprinter = make_maccs166_fingerprinter,
    )


# The number of bits depends on the parameters
def _get_num_bits(d):
    return d["fpSize"]

RDKitFingerprintFamily_v1 = _base.clone(
    name = "RDKit-Fingerprint/1",
    format_string = ("minPath=%(minPath)s maxPath=%(maxPath)s fpSize=%(fpSize)s "
                     "nBitsPerHash=%(nBitsPerHash)s useHs=%(useHs)s"),
    num_bits = _get_num_bits,
    make_fingerprinter = make_rdk_fingerprinter,
    )

###

RDKitMorganFingerprintFamily_v1 = _base.clone(
    name = "RDKit-Morgan/1",
    format_string = (
             "radius=%(radius)d fpSize=%(fpSize)s useFeatures=%(useFeatures)d "
             "useChirality=%(useChirality)d useBondTypes=%(useBondTypes)d"),
    num_bits = _get_num_bits,
    make_fingerprinter = make_morgan_fingerprinter,
    )

###

def _check_torsion_version(version):
    def make_fingerprinter(*args, **kwargs):
        if TORSION_VERSION != version:
            raise TypeError("This version of RDKit does not support the RDKit-Torsion/%s fingerprint" % (version,))
        return make_torsion_fingerprinter(*args, **kwargs)
    return make_fingerprinter

RDKitTorsionFingerprintFamily_v1 = _base.clone(
    name = "RDKit-Torsion/1",
    format_string = "fpSize=%(fpSize)s targetSize=%(targetSize)d",
    num_bits = _get_num_bits,
    make_fingerprinter = _check_torsion_version("1"),
    )

RDKitTorsionFingerprintFamily_v2 = _base.clone(
    name = "RDKit-Torsion/1",
    format_string = "fpSize=%(fpSize)s targetSize=%(targetSize)d",
    num_bits = _get_num_bits,
    make_fingerprinter = _check_torsion_version("2"),
    )

###

def _check_atom_pair_version(version):
    def make_fingerprinter(*args, **kwargs):
        if ATOM_PAIR_VERSION != version:
            raise TypeError("This version of RDKit does not support the RDKit-AtomPair/%s fingerprint" % (version,))
        return make_atom_pair_fingerprinter(*args, **kwargs)
    return make_fingerprinter

RDKitAtomPairFingerprintFamily_v1 = _base.clone(
    name = "RDKit-AtomPair/1",
    format_string = "fpSize=%(fpSize)s minLength=%(minLength)d maxLength=%(maxLength)d",
    num_bits = _get_num_bits,
    make_fingerprinter = _check_atom_pair_version("1"),
    )

RDKitAtomPairFingerprintFamily_v2 = _base.clone(
    name = "RDKit-AtomPair/2",
    format_string = "fpSize=%(fpSize)s minLength=%(minLength)d maxLength=%(maxLength)d",
    num_bits = _get_num_bits,
    make_fingerprinter = _check_atom_pair_version("2"),
    )
