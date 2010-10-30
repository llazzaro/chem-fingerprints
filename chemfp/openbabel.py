"Create OpenBabel fingerprints"

# Copyright (c) 2010 Andrew Dalke Scientific, AB (Gothenburg, Sweden)
# See the contents of "__init__.py" for full license details.

from __future__ import absolute_import

import sys
import os
import struct

import openbabel as ob

# This is the only thing which I consider to be public
__all__ = ["read_structures"]

# This is a "standard" size according to the struct module
# documentation, so the following is an excess of caution
if struct.calcsize("<I") != 4:
    raise AssertionError("The chemfp.ob module assumes 32 bit integers")


# OpenBabel doesn't currently have a function to return the version.
# I've brought this up on the list. I'm going to assume there will be
# one called "GetVersion()" in the future. Until then, I can fake it
# by reading the PDB output text.
try:
    from openbabel import GetVersion
except ImportError:
    def GetVersion():
        "GetVersion() -> the version string for the OpenBabel toolkit"
        obconversion = ob.OBConversion()
        obconversion.SetInFormat("smi")
        obconversion.SetOutFormat("pdb")
        obmol = ob.OBMol()

        obconversion.ReadString(obmol, "C")
        for line in obconversion.WriteString(obmol).splitlines():
            if "GENERATED BY OPEN BABEL" in line:
                return line.split()[-1]
        return "<unknown>"
_ob_version = GetVersion()

SOFTWARE = "OpenBabel/" + _ob_version


# OpenBabel fingerprints are stored as vector<unsigned int>.  On all
# the machines I use, ints have 32 bits.

# Bit lengths must be at least sizeof(int)*8 bits long and must be a
# factor of two. I have no idea why this is required.

# OpenBabel supports new fingerprints through a plugin system.  I got
# it working thanks to Noel O'Boyle's excellent work with Cinfony. I
# then found out that the OB API doesn't have any way to get the
# number of bits in the fingerprint. The size is rounded up to the
# next power of two, so FP4 (307 bits) needs 512 bits (16 ints)
# instead of 320 bits (10 ints). That means I can't even get close to
# guessing the bitsize.

# In the end, I hard-coded the supported fingerprints into the system.

class FingerprinterInfo(object):
    def __init__(self, name, description, ob_fingerprinter, num_bits):
        self.name = name
        self.description = description
        self.ob_fingerprinter = ob_fingerprinter
        self.num_bits = num_bits

        self.fps_params = "OpenBabel-" + name + "/1"

        # This is filled in later with the function which computes the
        # fingerprint as a string of bytes.
        self.calc_fp = None
        
_fingerprinter_table = {}

def _init():
    for name, num_bits in ( ("FP2", 1021),
                            ("FP3", 55),
                            ("FP4", 307) ):
        ob_fingerprinter = ob.OBFingerprint.FindFingerprint(name)
        description = ob_fingerprinter.Description()

        _fingerprinter_table[name] = FingerprinterInfo(
            name, description, ob_fingerprinter, num_bits)

    # Verify that OpenBabel was compiled with 32 bit integers
    n = _fingerprinter_table["FP2"].ob_fingerprinter.Getbitsperint()
    if n != 32:
        raise AssertionError(
            "The chemfp.ob module assumes OB fingerprints have 32 bits per integer")
_init()


############

# I could have written a more general function which created these but
# there's only three fingerprints lengths to worry about.

# This needs 128 bytes, for 1024 bits
# vectorUnsignedInt will contain 32 32-bit words = 1024 bits
def calc_FP2(mol, fp=None,
             get_fingerprint=_fingerprinter_table["FP2"].ob_fingerprinter.GetFingerprint,
             _pack_1024 = struct.Struct("<" + "I"*32).pack):
    if fp is None:
        fp = ob.vectorUnsignedInt()
    get_fingerprint(mol, fp)
    return _pack_1024(*fp)

# This needs 7 bytes, for 56 bits.
# vectorUnsignedInt will contain 2 32-bit words = 64 bits
def calc_FP3(mol, fp=None,
             get_fingerprint=_fingerprinter_table["FP3"].ob_fingerprinter.GetFingerprint,
             _pack_64 = struct.Struct("<II").pack):
    if fp is None:
        fp = ob.vectorUnsignedInt()
    get_fingerprint(mol, fp)
    return _pack_64(*fp)[:7]

# This needs 39 bytes, for 312 bits
# vectorUnsignedInt will contain 16 32-bit words = 512 bits
def calc_FP4(mol, fp=None,
             get_fingerprint=_fingerprinter_table["FP4"].ob_fingerprinter.GetFingerprint,
             _pack_512 = struct.Struct("<" + "I"*16).pack):
    if fp is None:
        fp = ob.vectorUnsignedInt()
    get_fingerprint(mol, fp)
    return _pack_512(*values)[:39]

_fingerprinter_table["FP2"].calc_fp = calc_FP2
_fingerprinter_table["FP3"].calc_fp = calc_FP3
_fingerprinter_table["FP4"].calc_fp = calc_FP4

#########


def _guess_format(filename):
    "This is an internal function"
    # Guess the format based on the filename extension
    base, ext = os.path.splitext(filename.lower())
    if ext in ("", "."):
        # No extension? Assume smiles
        return "smi"

    if ext != ".gz":
        # Not compressed? Use whatever was given.
        return ext[1:] # remove leading "."

    # Compressed? Get the next extension
    base, ext = os.path.splitext(base)
    if ext in ("", "."):
        # No extension? Assume smiles
        return "smi"

    # Use the 2nd level extension
    return ext[1:] # remove leading "."

def normalize_format(filename, format):
    """normalize_format(filename, format) -> normalized format name

    This is perhaps easiestly explained by example:
        ("input.dat", "smi.gz") -> "smi"   (format takes precedence)
        ("input.smi", "sdf") -> "sdf"
        ("INPUT.SDF", None) -> "sdf"       (otherwise use the filename)
        ("input.pdb.gz", None) -> "pdb"
        (None, None) -> "smi"              (otherwise it's SMILES)

    In words, get an OpenBabel format name given an input filename and
    format. The rules are:
      - Convert input strings to lower-case.
      - Remove the ".gz" suffix from each, if present.
      - If the format is not None, return it as the normalized format.
      - If the filename is not None, return the suffix.
      - If the filename has an extension, use that as the format.
      - Otherwise, it's in "smi" format.
    """
    if format:
        format = format.lower()
        # Ignore .gz - Babel handles those automatically
        if format.endswith(".gz"):
            format = format[:-3]
        return format

    if filename is not None:
        return _guess_format(filename)

    # Use SMILES by default for stdin
    return "smi"


def _get_ob_error(log):
    msgs = log.GetMessagesOfLevel(ob.obError)
    return "".join(msgs)

def read_structures(filename=None, format=None):
    """read_structures(filename, format) -> (title, OBMol) iterator 
    
    Iterate over structures from filename, returning the structure
    title and OBMol for each reacord. The structure is assumed to be
    in normalized_format(filename, format) format. If filename is None
    then this reads from stdin instead of the named file.
    """
    obconversion = ob.OBConversion()

    format = normalize_format(filename, format)
    if not obconversion.SetInFormat(format):
        # XXX better exception
        raise Exception("Unsupported format: %r" % (format,))
    
    obmol = ob.OBMol()

    if not filename:
        # OpenBabel's Python libary has no way to read from std::cin
        # Fake it through /dev/stdin for those OSes which support it.
        if not os.path.exists("/dev/stdin"):
            raise Exception("Unable to read from stdin on this platform")

        return _stdin_reader(obconversion, obmol)

    # Deal with OpenBabel's logging
    ob.obErrorLog.ClearLog()
    lvl = ob.obErrorLog.GetOutputLevel()
    ob.obErrorLog.SetOutputLevel(-1) # Suppress messages to stderr

    notatend = obconversion.ReadFile(obmol, filename)

    if ob.obErrorLog.GetErrorMessageCount():
        # The OB error messages are not that helpful. Do
        # some probing of my own before going to OB's message.
        try:
            open(filename).close()
        except IOError, err:
            raise SystemExit("Unable to open structure file %r: %s" %
                             (filename, err.strerror))

        # Okay, don't know what's going on. Report OB's error
        errmsg = _get_ob_error(ob.obErrorLog)
        raise SystemExit("Unable to get structures from %s:\n%s" %
                         (filename, errmsg))

    ob.obErrorLog.SetOutputLevel(lvl) # Revert to normal logging

    # We've opened the file. Switch to the iterator.
    return _file_reader(obconversion, obmol, notatend)

def _stdin_reader(obconversion, obmol):
    "This is an internal function"
    # Read structures from stdin.

    # The current release of scripting for OpenBabel doesn't let me
    # use C++'s std::cin as the input stream so I need to fake it
    # using "/dev/stdin". That works on Macs, Linux, and FreeBSD but
    # not Windows.

    # Python says that it's in charge of checking for ^C. When I pass
    # control over to OpenBabel, Python is still in charge of checking
    # for ^C, but it won't do anything until control returns to Python.

    # I found it annoying that if I started the program, which by
    # default expects SMILES from stdin, then I couldn't hit ^C to
    # stop it. My solution was to stay in Python until there's
    # information on stdin, and once that's happened, call OpenBabel.

    import select
    try:
        select.select([sys.stdin], [], [sys.stdin])
    except KeyboardInterrupt:
        raise SystemExit()

    # There's data. Pass parsing control into OpenBabel.
    notatend = obconversion.ReadFile(obmol, "/dev/stdin")
    while notatend:
        yield obmol.GetTitle(), obmol
        obmol.Clear()
        notatend = obconversion.Read(obmol)

def _file_reader(obconversion, obmol, notatend):
    while notatend:
        yield obmol.GetTitle(), obmol
        obmol.Clear()
        notatend = obconversion.Read(obmol)


