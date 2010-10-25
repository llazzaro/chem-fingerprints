"Create RDKit fingerprints"

# Copyright (c) 2010 Andrew Dalke Scientific, AB (Gothenburg, Sweden)
# See the contents of "__init__.py" for full license details.

from __future__ import absolute_import

import os
import sys
import gzip

import rdkit
from rdkit import Chem
from rdkit.Chem.MACCSkeys import GenMACCSKeys


from . import sdf_reader, decoders

# These are the things I consider to be public
__all__ = ["read_structures", "read_smiles_structures", "read_sdf_structures"]

######## Some shenanigans to get a version field

if hasattr(rdkit, "version"):
    # This will be available in the next release of RDKit
    SOFTWARE = "RDKit/" + rdkit.version()
else:
    # This distribution does not have version().
    # Guess based on the evidence.

    # This was added in Release_Q22009_1
    import rdkit.Chem.AtomPairs.Torsions
    if not hasattr(Chem.AtomPairs.Torsions, "GetHashedTopologicalTorsionFingerprint"):
        SOFTWARE = "RDKit/unknown"
    else:
        # Either Release_Q22009_1 or Release_Q32009_1
        # Q3 removed the DaylightFingerprint function.
        if hasattr(Chem, "DaylightFingerprint"):
            # In 2010 the version numbers were put into lexical order.
            # Keep with the same principle.
            SOFTWARE = "RDKit/2009Q2_1"
        else:
            SOFTWARE = "RDKit/2009Q3_1"

##### Convert from RDKit explicit (dense) fingerprints to byte strings



#########

def normalize_format(filename, format):
    """normalize_format(filename, format) -> (normalized_format, is_compressed)

    This is perhaps easiestly explained by example:
        ("input.dat", "smi.gz") -> ("smi", 1)   (format takes precedence)
        ("input.smi", "sdf") -> ("sdf", 0)
        ("INPUT.SDF.GZ", None) -> ("sdf", 1)    (otherwise use the filename)
        ("input.pdb", None) -> ("pdb", 0)
        (None, None) -> ("smi", 0)              (otherwise it's uncompressed SMILES)

    """
    if format is not None:
        compressed = 0
        if format.endswith(".gz"):
            compressed = 1
            format = format[:-3]
        return format, compressed

    if filename is None:
        # Reading from stdin, with no specified format
        # Assume SMILES
        return "smi", 0

    base, ext = os.path.splitext(filename.lower())
    compressed = 0
    if ext == ".gz":
        compressed = 1
        base, ext = os.path.splitext(base)
        
    if ext in (".sdf", ".mol", ".sd", ".mdl"):
        return "sdf", compressed

    elif ext in (".smi", ".can", ".smiles", ".ism"):
        return "smi", compressed

    # When unknown, guess SMILES
    return "smi", compressed

# While RDKit has a SMILES file parser, it doesn't handle reading from
# stdin or from compressed file. I wanted to support those as well, so
# ended up not using Chem.SmilesMolSupplier.

def read_smiles_structures(infile):
    """read_smiles_structures(infile) -> (title, Chem.Mol) iterator

    Iterate through each record in the SMILES file, returning the
    2-ple of (title, rdkit.Chem.Mol). Each record is a single line of
    fields, separated by whitespace. The first column is the SMILES
    used to create the molecule. It must be present. The second column
    is the title. If not present, the record number (starting with
    "1") is used as the title.
    """
    for lineno, line in enumerate(infile):
        words = line.split()
        if not words:
            continue
        mol = Chem.MolFromSmiles(words[0])
        if len(words) == 1:
            yield str(lineno+1), mol
        else:
            yield words[1], mol

def bad_record(message):
    raise TypeError(message)
    #if message[-1:] != "\n":
    #    message = message + "\n"
    #sys.stderr.write(message)

def read_sdf_structures(infile, filename=None, bad_record=bad_record):
    """read_sdf_structure(infile) -> (title, Chem.Mol) iterator

    Iterate through each record in the SD file, returning the 2-ple of
    (title, rdkit.Chem.Mol). The title comes from the "_Name" property
    of the molecule object.
    """
    if filename is None:
        filename = getattr(infile, "name", None)
    if filename is None:
        errtxt = "Could not parse {what}"
    else:
        errtxt = "Could not parse {what} of {filename!r}"
    for record_lines in sdf_reader.read_sdf_records(infile):
        text = "".join(record_lines)
        mol = Chem.MolFromMolBlock(text)
        if mol is None:
            # This was not a molecule?
            bad_record(errtxt.format(what=record_lines.what(), filename=filename))
        else:
            yield mol.GetProp("_Name"), mol

# This class helps the case when someone is entering structure
# by-hand. (Most likely to occur with SMILES input). They would like
# to see the result as soon as a record is entered. But normal
# interation reader grabs a buffer of input to process, and not a
# line. It's faster that way. The following adapter supports the
# iterator protocol but turns it into simple readlines(). This will be
# slower but since do it only if stdin is a tty, there shouldn't be a
# problem.
class _IterUsingReadline(object):
    "Internal class for iterating a line at a time from tty input"
    def __init__(self, fileobj):
        self.fileobj = fileobj
    def __iter__(self):
        return iter(self.fileobj.readline, "")

def _open(filename, compressed):
    "Internal function to open the given filename, which might be compressed"
    if filename is None:
        if compressed:
            return gzip.GzipFile(fileobj=sys.stdin, mode="r")
        else:
            # Python's iter reads a block.
            # When someone types interactively, read only a line.
            if sys.stdin.isatty():
                return _IterUsingReadline(sys.stdin)
            else:
                return sys.stdin

    if compressed:
        return gzip.GzipFile(filename, "r")
    return open(filename, "rU")
    

def read_structures(filename, format):
    """read_structures(filename, format) -> (title, rdkit.Chem.Mol) iterator 
    
    Iterate over structures from filename, returning the structure
    title and RDKit.Chem.Mol for each reacord. The structure is
    assumed to be in normalized_format(filename, format) format and
    compression. If filename is None then this reads from stdin
    instead of the named file.
    """
    format, compressed = normalize_format(filename, format)
    if format == "sdf":
        # I timed the native reader at 25.8 seconds (best of 25.8, 26.0, 25.8)
        # and the Python reader at 26.0 seconds (best of 26.0, 26.0, 26.1)
        # While the native reader was 1% faster, I prefer the consistancy
        # of having only one reader.
        #
        #if (not compressed) and (filename is not None):
        #    supplier = Chem.SDMolSupplier(filename)
        #    def native_sdf_reader():
        #        for mol in supplier:
        #            if mol is None:
        #                print >>sys.stderr, "Missing? after", title
        #            else:
        #                title = mol.GetProp("_Name")
        #                yield title, mol
        #    return native_sdf_reader()

        return read_sdf_structures(_open(filename, compressed))

    elif format == "smi":
        # I timed the native reader at 31.6 seconds (best of 31.6, 31.7, 31.7)
        # and the Python reader at 30.8 seconds (best of 30.8, 30.9, and 31.0)
        # Yes, the Python reader is faster and using it gives me better consistency
        #
        #if (not compressed) and (filename is not None):
        #    supplier = Chem.SmilesMolSupplier(filename, delimiter=" \t", titleLine=False)
        #    def native_smiles_reader():
        #        for mol in supplier:
        #            yield mol.GetProp("_Name"), mol
        #    return native_smiles_reader()
        return read_smiles_structures(_open(filename, compressed))

    else:
        raise TypeError("Unknown format")

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

def make_rdk_fingerprinter(num_bits=NUM_BITS, min_path=MIN_PATH, max_path=MAX_PATH,
                           bits_per_hash=BITS_PER_HASH, use_Hs=USE_HS):
    # I've already done checks for this, but I'm assuming that people
    # will use this as part of the public API.
    if not (num_bits > 0): raise TypeError("num_bits must be positive")
    if not (min_path > 0): raise TypeError("min_path must be positive")
    if not (max_path >= min_path):
        raise TypeError("max_path cannot be smaller than min_path")
    if not (bits_per_hash > 0): raise TypeError("bits_per_hash must be positive")

    def rdk_fingerprinter(mol):
        fp = Chem.RDKFingerprint(
            mol, minPath=min_path, maxPath=max_path, fpSize=num_bits,
            nBitsPerHash=bits_per_hash, useHs=use_Hs)
        return decoders.from_binary_lsb(fp.ToBitString())[1]
    return rdk_fingerprinter

# Use this to make the parameters for the topological fingerprints
RDK_PARAMS = ("RDKit-Fingerprint/1 minPath={min_path} maxPath={max_path} fpSize={num_bits} "
              "nBitsPerHash={bits_per_hash} useHs={use_Hs}")

format_rdk_params = RDK_PARAMS.format

########### The MACCS fingerprinter

def maccs166_fingerprinter(mol):
    fp = GenMACCSKeys(mol)
    return decoder.from_binary_lsb(fp)

