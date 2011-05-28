"""Create OpenEye fingerprints


"""

# Copyright (c) 2010 Andrew Dalke Scientific, AB (Gothenburg, Sweden)
# Licensed under "the MIT license"
# See the contents of COPYING or "__init__.py" for full license details.

from __future__ import absolute_import

import sys
import os
import errno
import select
import ctypes
import warnings
import errno

from openeye.oechem import *
from openeye.oegraphsim import *

__all__ = ["read_structures", "get_path_fingerprinter", "get_maccs_fingerprinter"]

from . import types

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

# Mapping from atype string to flag (includes AtomDefault)
_atypes = dict(_atype_flags)
_btypes = dict(_btype_flags)

_atypes["DefaultAtom"] = OEFPAtomType_DefaultAtom
_btypes["DefaultBond"] = OEFPBondType_DefaultBond

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
                raise TypeError("Missing %s flag" % (a_or_b,))
            raise TypeError("Unknown %s type %r" % (a_or_b, word))
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

# I could make the final string sorter by using DefaultAtom/
# DefaultBond but OpenEye might change that value in the future. I
# think this is slightly safer.

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


_maccs_decoders = {"numbits": int,
                   "minbonds": int,
                   "maxbonds": int,
                   "atype": atom_description_to_value,
                   "btype": bond_description_to_value}

def decode_path_parameters(parameters):
    assert len(parameters) == len(_maccs_decoders)
    kwargs = {}
    for name, decoder in _maccs_decoders.items():
        value = parameters[name]
        kwargs[name] = decoder(value)
    return kwargs

_maccs_encoders = {"numbits": str,
                   "minbonds": str,
                   "maxbonds": str,
                   "atype": atom_value_to_description,
                   "btype": bond_value_to_description}

def encode_path_parameters(kwargs):
    assert len(kwargs) == len(_maccs_encoders)
    parameters = {}
    for name, encoder in _maccs_encoders.items():
        value = kwargs[name]
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

def get_path_fingerprinter(numbits, minbonds, maxbonds, atype, btype):
    # Extra level of error checking since I expect people will think
    # of this as part of the public API.
    if not (16 <= numbits <= 65536):
        raise TypeError("numbits must be between 16 and 65536 (inclusive)")
    if not (0 <= minbonds):
        raise TypeError("minbonds must be 0 or greater")
    if not (minbonds <= maxbonds):
        raise TypeError("maxbonds must not be smaller than minbonds")

    # XXX validate the atype and type values? Should just
    # be a simple bitwise-and and test for 0.
    
    fp = OEFingerPrint()
    fp.SetSize(numbits)
    data_location = int(fp.GetData())
    num_bytes = (numbits+7)//8

    def path_fingerprinter(mol):
        OEMakePathFP(fp, mol, numbits, minbonds, maxbonds, atype, btype)
        return ctypes.string_at(data_location, num_bytes)

    return path_fingerprinter

# Believe it or not, reusing the preallocated fingerprint does help
# the performance.
def get_maccs_fingerprinter():
    fp = OEFingerPrint()
    # SetSize() now to force space allocation, so I only need one GetData()
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
    if format_flag is None and fmt.startswith("."):
        # Some OE tools allow ".sdf" as the format. In the interests
        # of compatibility, I support that as well as the more
        # acceptable "sdf".
        format_flag = _formats.get(fmt[1:], None)
        if format_flag is not None:
            warnings.warn("format name %(format)r should be written %(better)r" %
                          dict(format=format, better=format[1:]), DeprecationWarning)

    if format_flag is None:
        raise UnknownFormat(format)

    def _apply_format(ifs):
        ifs.SetFormat(format_flag)
        if is_compressed:
            ifs.Setgz(is_compressed)
    return _apply_format

def _open_ifs(filename, _apply_format):
    ifs = oemolistream()
    _apply_format(ifs)
    if not ifs.open(filename):
        # Let Python try to do better error reporting.
        open(filename).close()
        # If that didn't work, give up and fake it.
        # (Did manual coverage testing for this. The test cases I can
        # think of, like tricky timing, are too tricky.)
        raise IOError(errno.EIO, "OEChem cannot open the file", filename)
    return ifs

# This code is a bit strange. It needs to do eager error checking but
# lazy parsing. That is, it needs to check right away that the file
# can be opened (if it exists) and the format is understood. But it
# can wait until later to actually parse the files.

def read_structures(filename=None, format=None):
    # Check that that the format is known
    _apply_format = _get_format_setter(format)

    # Input is from a file
    if filename is not None:
        ifs = _open_ifs(filename, _apply_format)
        return _iter_structures(ifs)
    else:
        # Input is from stdin
        return _stdin_check(_apply_format)

def _iter_structures(ifs):
    for mol in ifs.GetOEGraphMols():
        yield mol.GetTitle(), mol

# The reason for this is to allow ^C to work if there hasn't yet been
# any input. Why? When Python calls out to a C++ function it blocks
# control-C. It's not possible to get a KeyboardInterrupt until the
# C++ function returns. Consider someone using this script who omitted
# a filename by accident. If OEChem has control of stdin then ^C does
# not work. I found that frustrating, so this is a workaround. When
# reading from stdin, don't dispatch to OEChem until there's input.

_USE_SELECT = True
def _stdin_check(_apply_format):
    if _USE_SELECT:
        try:
            select.select([sys.stdin], [], [sys.stdin])
        except KeyboardInterrupt:
            raise SystemExit()
    ifs = oemolistream()
    ifs.open()
    _apply_format(ifs)
    for mol in ifs.GetOEGraphMols():
        yield mol.GetTitle(), mol

############# Methods to get the right structure readers

def read_maccs166_fingerprints_v1(source=None, format=None, kwargs={}):
    assert not kwargs, kwargs
    # The OEChem interface only handles stdin and filenames
    if not (isinstance(source, basestring) or source is None):
        raise NotImplementedError

    fingerprinter = get_maccs_fingerprinter()
    structure_reader = read_structures(source, format)

    def read_oechem_maccs_structure_fingerprints():
        for (title, mol) in structure_reader:
            yield fingerprinter(mol), title
    return read_oechem_maccs_structure_fingerprints()

def read_path_fingerprints_v1(source=None, format=None, kwargs={}):
    # The OEChem interface only handles stdin and filenames
    if not (isinstance(source, basestring) or source is None):
        raise NotImplementedError
    
    fingerprinter = get_path_fingerprinter(**kwargs)
    structure_reader = read_structures(source, format)
    
    def read_oechem_path_structure_fingerprints():
        for (title, mol) in structure_reader:
            yield fingerprinter(mol), title
    return read_oechem_path_structure_fingerprints()

############# Used when generate the FPS header
class OpenEyePathFingerprinter_v1(types.Fingerprinter):
    name = "OpenEye-Path/1"
    format_string = ("numbits=%(numbits)s minbonds=%(minbonds)s "
                     "maxbonds=%(maxbonds)s atype=%(atype)s btype=%(btype)s")
    software = SOFTWARE
    def __init__(self, kwargs):
        self.num_bits = kwargs["numbits"]
        super(OpenEyePathFingerprinter_v1, self).__init__(kwargs)

    @classmethod
    def from_parameters(cls, parameters):
        return cls(decode_path_parameters(parameters))

    def _encode_parameters(self):
        return encode_path_parameters(self.kwargs)

    _get_reader = staticmethod(read_path_fingerprints_v1)

    
class OpenEyeMACCSFingerprinter_v1(types.Fingerprinter):
    name = "OpenEye-MACCS166/1"
    num_bits = 166
    software = SOFTWARE

    _get_reader = staticmethod(read_maccs166_fingerprints_v1)
