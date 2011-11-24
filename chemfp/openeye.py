"""Create OpenEye fingerprints


"""

# Copyright (c) 2010 Andrew Dalke Scientific, AB (Gothenburg, Sweden)
# Licensed under "the MIT license"
# See the contents of COPYING or "__init__.py" for full license details.

from __future__ import absolute_import

import sys
import os
import errno
import ctypes
import warnings
import errno

from openeye.oechem import *
from openeye.oegraphsim import *

from . import ParseError
from . import types
from . import io
from . import error_handlers
from . import argparse

__all__ = ["read_structures", "get_path_fingerprinter", "get_maccs_fingerprinter"]


class UnknownFormat(KeyError):
    def __str__(self):
        return "Unknown format %r" % (self.args[0],)

############# Used when generate the FPS header

SOFTWARE = "OEGraphSim/%(release)s (%(version)s)" % dict(
    release = OEGraphSimGetRelease(),
    version = OEGraphSimGetVersion())



##### Handle the atom and bond type flags for path fingerprints

# The atom and bond type flags can be specified on the command-line
# 
#   --atype=DefaultAtom --btype=BondOrder,InRing
#   --atype AtomicNumber,InRing  --btype DefaultBond,InRing
#
# The type fields may be separated by either a "," or a "|".
# The relevant OpenEye function (OEGetFPAtomType() and
# OEGetFPBondType()) use a "|" but that requires escaping for
# the shell, so I support "," as well.

# There's another conversion of the integer type values into a string
# representation, used when generating the canonical form of the
# generation parameters for the FPS output. That case uses "|"
# (and not ",") and omits the DefaultAtom and DefaultBond name.
# The result is easier to parse with the OpenEye API functions.

_atype_flags = [(OEGetFPAtomType(atype), atype) for atype in
                (OEFPAtomType_Aromaticity,
                OEFPAtomType_AtomicNumber,
                OEFPAtomType_Chiral,
                OEFPAtomType_EqAromatic,
                OEFPAtomType_EqHalogen,
                OEFPAtomType_FormalCharge,
                OEFPAtomType_HvyDegree,
                OEFPAtomType_Hybridization,
                OEFPAtomType_InRing)]
_btype_flags = [(OEGetFPBondType(btype), btype) for btype in
                (OEFPBondType_BondOrder,
                OEFPBondType_Chiral,
                OEFPBondType_InRing)]

_atypes = dict(_atype_flags)
_btypes = dict(_btype_flags)

# I support the DefaultAtom and DefaultBond values but I'm worried
# that OpenEye will change the composition of those in the future. I
# talked with Krisztina Boda and she says that that won't change, but
# I don't necessarily trust future maintainers. To minimize problems,
# I'll hard-code them the current Default* values.

_atype_flags.insert(0, ("DefaultAtom", 191))
_atypes["DefaultAtom"] = 191 # OEFPAtomType_DefaultAtom

_btype_flags.insert(0, ("DefaultBond", 3))
_btypes["DefaultBond"] = 3 # OEFPBondType_DefaultBond


## Go from a "," or "|" separated text field to an integer value
# Removes extra whitespace, but none should be present.

def _get_type_value(a_or_b, table, description):
    value = 0
    # Allow both "|" and "," as separators
    # (Consistent with OEChem)
    description = description.replace("|", ",")
    for word in description.split(","):
        word = word.strip()
        try:
            value |= table[word]
        except KeyError:
            if not word:
                raise ValueError("Missing %s flag" % (a_or_b,))
            raise ValueError("Unknown %s type %r" % (a_or_b, word))
    return value

def atom_description_to_value(description):
    """atom_description_to_value(description) -> integer

    Convert an atom description like FormalCharge,EqHalogen
    or FormalCharge|EqHalogen into its atom type value.

    This is similar to OEGetFPAtomType except both "|" and "," are
    allowed seperators and "AtomDefault" is an allowed term.
    """
    return _get_type_value("atom", _atypes, description)

def bond_description_to_value(description):
    """bond_description_to_value(description) -> integer

    Convert an bond description like BondOrder,Chiral
    or BondOrder|Chiral into its bond type value.

    This is similar to OEGetFPBondType except both "|" and "," are
    allowed seperators and "BondDefault" is an allowed term.
    """
    return _get_type_value("bond", _btypes, description)

## Go from an integer value into a canonical description

# I could use OEGetFPAtomType() and OEGetFPBondType() but I wanted
# something which has a fixed sort order even for future releases,
# which isn't part of those functions.

def _get_type_description(a_or_b, flags, value):
    words = []
    for (word, flag) in flags:
        if flag & value == flag:
            # After over 12 years of full-time use of Python,
            # I finally have a chance to use the xor operator.
            value = value ^ flag
            words.append(word)
    if value != 0:
        raise AssertionError("Unsupported %s value" % (a_or_b,))
    return "|".join(words)


def atom_value_to_description(value):
    """atom_value_to_description(value) -> string

    Convert from an atom type string into its text description,
    separated by "|"s. The result are compatible with
    OEGetFPAtomType and are in canonical order.
    """
    return _get_type_description("atom", _atype_flags, value)

def bond_value_to_description(value):
    """bond_value_to_description(value) -> string

    Convert from a bond type string into its text description,
    separated by "|"s. The result are compatible with
    OEGetFPBontType and are in canonical order.
    """
    return _get_type_description("bond", _btype_flags, value)


def decode_path_parameters(parameters):
    fingerprinter_kwargs = _maccs_defaults.copy()
    for name, value in parameters:
        if name not in _maccs_decoders:
            raise ValueError("Unknown OpenEye-Path parameter %r" % (name,))
        decoder = _maccs_decoders[name]
        fingerprinter_kwargs[name] = decoder(value)
    return fingerprinter_kwargs

_maccs_encoders = {"numbits": str,
                   "minbonds": str,
                   "maxbonds": str,
                   "atype": atom_value_to_description,
                   "btype": bond_value_to_description}

def encode_path_parameters(fingerprinter_kwargs):
    assert len(fingerprinter_kwargs) == len(_maccs_encoders)
    parameters = {}
    for name, encoder in _maccs_encoders.items():
        value = fingerprinter_kwargs[name]
        parameters[name] = encoder(value)
    return parameters


##### Create a function which generate fingerprints

# I use functions which return functions because it was a nice way to
# hide the differences between the two fingerprinters. I also found
# that I can save a bit of time by not creating a new fingerprint each
# time. The measured speedup is about 2% for MACCS166 and 6% for path
# fingerprints.

# Just like the OEGraphMol, these fingerprints must not be reused or
# stored. They are mutated EVERY TIME. They are NOT thread-safe.
# If you need to use these in multiple threads, then make multiple
# fingerprinters.

# Believe it or not, reusing the preallocated fingerprint measurably
# helps the performance.

def get_path_fingerprinter(numbits, minbonds, maxbonds, atype, btype):
    # Extra level of error checking since I expect people will think
    # of this as part of the public API.
    if not (16 <= numbits <= 65536):
        raise ValueError("numbits must be between 16 and 65536 (inclusive)")
    if not (0 <= minbonds):
        raise ValueError("minbonds must be 0 or greater")
    if not (minbonds <= maxbonds):
        raise ValueError("maxbonds must not be smaller than minbonds")

    # XXX validate the atype and type values?
    # It's a simple mask against the | of all possible value, then test for 0.
    # However, I'm not sure what to report as the error message.
    
    fp = OEFingerPrint()
    fp.SetSize(numbits)
    data_location = int(fp.GetData())
    num_bytes = (numbits+7)//8

    def path_fingerprinter(mol):
        OEMakePathFP(fp, mol, numbits, minbonds, maxbonds, atype, btype)
        return ctypes.string_at(data_location, num_bytes)

    return path_fingerprinter

def get_maccs_fingerprinter():
    fp = OEFingerPrint()
    # Call SetSize() now to force space allocation, so I only need one GetData()
    fp.SetSize(166)
    data_location = int(fp.GetData())
    num_bytes = (166+7)//8
    
    def maccs_fingerprinter(mol):
        OEMakeMACCS166FP(fp, mol)
        return ctypes.string_at(data_location, num_bytes)

    return maccs_fingerprinter

### A note on fingerprints and ctypes.string_at

# The FPS format and OEFingerPrint.GetData() values used identical bit
# and byte order. Bytes are in little-endian order and bits are in
# big-endian order. That means I can use GetData() to get the
# underlying C storage area, use ctypes to turn that into a Python
# string, which I then hex encode.

# The other option is to use OEFingerPrint.ToHexString(). But that's
# pure little endian, so I would need a transposition to make the bits
# be what I want them to be. OEChem's hex strings also end with a flag
# which says how many extra bits to trim, which I don't need since I
# handle it a different way.

# Here's some info about the bit order, which I tested by setting a
# few bits though the API then seeing what changed in the hex string
# and in the underlying GetData() field.

# The bit pattern
#   01234567 89ABCDEF  pure little endian
#   10011100 01000011
#   
#         93 2C    using ToHexString()  (pure little endian)
#       0x39 c2    using hex(ord(GetData())) (litle endian byte, big endian bit)
#
#   76543210 FEDCBA98
#   00111001 11000010   little endian byte, big endian bit


################ Handle formats

# Make format names to OEChem format types
_formats = {
    "smi": OEFormat_SMI,
    "ism": OEFormat_ISM,
    "can": OEFormat_CAN,

    "sdf": OEFormat_SDF,
    "mol": OEFormat_SDF,

    "skc": OEFormat_SKC,
    "mol2": OEFormat_MOL2,
    "mmod": OEFormat_MMOD,
    "oeb": OEFormat_OEB,
    "bin": OEFormat_BIN,
}

# Some trickiness to verify that the format specification is
# supported, but without doing anything (yet) to set those flags.

# I return a function which will set the file stream parameters
# correctly.

def _do_nothing(ifs):
    pass

# Format is something like ".sdf.gz" or "pdb" or "smi.gz"

def _get_format_setter(format=None):
    if format is None:
        return _do_nothing
    fmt = format.lower()

    is_compressed = 0
    if fmt.endswith(".gz"):
        is_compressed = 1
        fmt = fmt[:-3]  # Should be something like ".sdf" or "sdf" or "smi"

    format_flag = _formats.get(fmt, None)
    if format_flag is None:
        raise ValueError("Unsupported format %r" % (format,))

    def set_format(ifs):
        ifs.SetFormat(format_flag)
        if is_compressed:
            ifs.Setgz(is_compressed)
    return set_format


def _open_stdin(set_format, aromaticity_flavor):
    ifs = oemolistream()
    ifs.open()
    set_format(ifs)
    
    if aromaticity_flavor is not None:
        flavor = ifs.GetFlavor(ifs.GetFormat())
        flavor |= aromaticity_flavor
        ifs.SetFlavor(ifs.GetFormat(), flavor)

    return ifs

def _open_ifs(filename, set_format, aromaticity_flavor):
    ifs = oemolistream()

    if not ifs.open(filename):
        # Let Python try to do better error reporting.
        open(filename).close()
        # If that didn't work, give up and fake it.
        # (Did manual coverage testing for this. The test cases I can
        # think of, like tricky timing, are too tricky.)
        raise IOError(errno.EIO, "OEChem cannot open the file", filename)

    set_format(ifs)

    if aromaticity_flavor is not None:
        flavor = ifs.GetFlavor(ifs.GetFormat())
        flavor |= aromaticity_flavor
        ifs.SetFlavor(ifs.GetFormat(), flavor)
    
    return ifs

# This code is a bit strange. It needs to do eager error checking but
# lazy parsing. That is, it needs to check right away that the file
# can be opened (if it exists) and the format is understood. But it
# can wait until later to actually parse the files.

_aromaticity_sorted = (
#    ("default", None),
    ("openeye", OEIFlavor_Generic_OEAroModelOpenEye),
    ("daylight", OEIFlavor_Generic_OEAroModelDaylight),
    ("tripos", OEIFlavor_Generic_OEAroModelTripos),
    ("mdl", OEIFlavor_Generic_OEAroModelMDL),
    ("mmff", OEIFlavor_Generic_OEAroModelMMFF),
    )
_aromaticity_flavors = dict(_aromaticity_sorted)
_aromaticity_flavor_names = [pair[0] for pair in _aromaticity_sorted]
_aromaticity_flavors[None] = OEIFlavor_Generic_OEAroModelOpenEye
del _aromaticity_sorted, pair

# If unspecified, use "openeye" (this is what OEChem does internally)
_aromaticity_flavors[None] = _aromaticity_flavors["openeye"]   # Allow "None"

def is_valid_format(filename, format):
    format_name, compression = io.normalize_format(filename, format,
                                                   default=("smi", ""))
    if compression not in ("", ".gz"):
        return False
    try:
        _get_format_setter(format_name + compression)
        return True
    except ValueError:
        return False

def is_valid_aromaticity(aromaticity):
    return aromaticity in _aromaticity_flavors

# Part of the code (parameter checking, opening the file) are eager.
# Actually reading the structures is lazy.

def read_structures(filename=None, format=None, id_tag=None, aromaticity=None, errors="strict"):
    try:
        aromaticity_flavor = _aromaticity_flavors[aromaticity]
    except KeyError:
        raise ValueError("Unsupported aromaticity name %r" % (aromaticity,))
    error_handler = error_handlers.get_parse_error_handler(errors)

    # Check that that the format is known
    format_name, compression = io.normalize_format(filename, format,
                                                   default=("smi", ""))

    if compression not in ("", ".gz"):
        raise ValueError("Unsupported compression type for %r" % (filename,))

    set_format = _get_format_setter(format_name + compression)

    # Input is from a file
    if filename is None:
        ifs = _open_stdin(set_format, aromaticity_flavor)
        filename_repr = "<stdin>"
    else:
        ifs = _open_ifs(filename, set_format, aromaticity_flavor)
        filename_repr = repr(filename)

    # Only SD files can take the id_tag
    if ifs.GetFormat() != OEFormat_SDF:
        id_tag = None

    # Lazy structure reader
    return _iter_structures(ifs, id_tag, filename_repr, error_handler)

def _iter_structures(ifs, id_tag, filename_repr, error_handler):
    def where():
        return " for record #%d of %s" % (recno+1, filename_repr)
        
    if id_tag is None:
        for recno, mol in enumerate(ifs.GetOEGraphMols()):
            title = mol.GetTitle()
            id = io.remove_special_characters_from_id(title)
            if not id:
                error_handler("Missing title" + where())
                continue
            yield id, mol
    else:
        for recno, mol in enumerate(ifs.GetOEGraphMols()):
            dirty_id = OEGetSDData(mol, id_tag)
            if not dirty_id:
                if not OEHasSDData(mol, id_tag):
                    error_handler("Missing id tag %r%s" % (id_tag, where()))
                    continue
            id = io.remove_special_characters_from_id(dirty_id)
            if not id:
                msg = "Empty id tag %r" % (id_tag,)
                error_handler(msg + where())
                continue
            yield id, mol


def _read_fingerprints(structure_reader, fingerprinter):
    for (id, mol) in structure_reader:
        yield id, fingerprinter(mol)

from .types import FingerprintFamilyConfig, positive_int, nonnegative_int, zero_or_one

def _read_structures(metadata, source, format, id_tag, errors):
    return read_structures(source, format, id_tag=id_tag,
                           aromaticity=metadata.aromaticity, errors=errors)

def _correct_numbits(s):
    try:
        if not s.isdigit():
            raise ValueError
        i = int(s)
        if not (16 <= i <= 65536):
            raise ValueError
    except ValueError:
        raise argparse.ArgumentError(None, "numbits must be between 16 and 65536 bits")
    return i

_base = FingerprintFamilyConfig(
    software = SOFTWARE,
    read_structures = _read_structures)

OpenEyePathFingerprintFamily_v1 = _base.clone(
    name = "OpenEye-Path/1",
    format_string = ("numbits=%(numbits)s minbonds=%(minbonds)s "
                     "maxbonds=%(maxbonds)s atype=%(atype)s btype=%(btype)s"),
    num_bits = lambda d: d["numbits"],
    make_fingerprinter = get_path_fingerprinter)

_path = OpenEyePathFingerprintFamily_v1
_path.add_argument("numbits", type=_correct_numbits, metavar="INT", default=4096,
                   help="number of bits in the path fingerprint")

_path.add_argument("minbonds", type=nonnegative_int("minbonds"), metavar="INT", default=0,
                   help="minimum number of bonds in the path")

_path.add_argument("maxbonds", type=nonnegative_int("maxbonds"), metavar="INT", default=5,
                   help="maximum number of bonds in the path")

def _atom_description(s):
    try:
        return atom_description_to_value(s)
    except ValueError, x:
        raise argparse.ArgumentError(None, str(x))

_path.add_argument("atype", type=_atom_description, encoder=atom_value_to_description,
                   help="atom type", default="DefaultAtom")

def _bond_description(s):
    try:
        return bond_description_to_value(s)
    except ValueError, x:
        raise argparse.ArgumentError(None, str(x))



_path.add_argument("btype", type=_bond_description, encoder=bond_value_to_description,
                   help="bond type", default="DefaultBond")

    
OpenEyeMACCSFingerprintFamily_v1 = _base.clone(
    name = "OpenEye-MACCS166/1",
    num_bits = 166,
    make_fingerprinter = get_maccs_fingerprinter)
