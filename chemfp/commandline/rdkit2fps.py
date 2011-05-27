# Copyright (c) 2010 Andrew Dalke Scientific, AB (Gothenburg, Sweden)
import sys
from chemfp import argparse, io, rdkit, types
from . import cmdsupport

########### Configure the command-line parser

epilog = """\

This program guesses the input structure format based on the filename
extension. If the data comes from stdin, or the extension name us
unknown, then use "--in" to change the default input format. The
supported format extensions are:

  File Type      Valid FORMATs (use gz if compressed)
  ---------      ------------------------------------
   SMILES        smi, ism, can, smi.gz, ism.gz, can.gz
   SDF           sdf, mol, sd, mdl, sdf.gz, mol.gz, sd.gz, mdl.gz
"""

parser = argparse.ArgumentParser(
    description="Generate FPS fingerprints from a structure file using RDKit",
    epilog=epilog,
    formatter_class=argparse.RawDescriptionHelpFormatter,    
    )
rdk_group = parser.add_argument_group("RDKit topological fingerprints")
rdk_group.add_argument("--RDK", action="store_true",
                       help="generate RDK fingerprints (default)")
rdk_group.add_argument(
    "--fpSize", type=int, metavar="INT", default=rdkit.NUM_BITS,
    help="number of bits in the fingerprint (default=%d)" % rdkit.NUM_BITS)
rdk_group.add_argument(
    "--minPath", type=int, metavar="INT", default=rdkit.MIN_PATH,
    help="minimum number of bonds to include in the subgraphs (default=%d)" % rdkit.MIN_PATH)
rdk_group.add_argument(
    "--maxPath", type=int, metavar="INT", default=rdkit.MAX_PATH,
    help="maximum number of bonds to include in the subgraphs (default=%d)" % rdkit.MAX_PATH)
rdk_group.add_argument(
    "--nBitsPerHash", type=int, metavar="INT", default=rdkit.BITS_PER_HASH,
    help="number of bits to set per path (default=%d)" % rdkit.BITS_PER_HASH)
rdk_group.add_argument(
    "--useHs", type=int, default=1,
    help="information about the number of hydrogens on each atom")

maccs_group = parser.add_argument_group("166 bit MACCS substructure keys")
maccs_group.add_argument(
    "--maccs166", action="store_true", help="generate MACCS fingerprints")

substruct_group = parser.add_argument_group("881 bit substructure keys")
substruct_group.add_argument(
    "--substruct", action="store_true", help="generate ChemFP substructure fingerprints")

rdmaccs_group = parser.add_argument_group("ChemFP version of the 166 bit RDKit/MACCS keys")
rdmaccs_group.add_argument(
    "--rdmaccs", action="store_true", help="generate ChemFP RDKit/MACCS")

parser.add_argument(
    "--in", metavar="FORMAT", dest="format",
    help="input structure format (default guesses from filename)")
parser.add_argument(
    "-o", "--output", metavar="FILENAME",
    help="save the fingerprints to FILENAME (default=stdout)")
parser.add_argument(
    "filename", nargs="?", help="input structure file (default is stdin)", default=None)


def main(args=None):
    args = parser.parse_args(args)

    cmdsupport.mutual_exclusion(parser, args, "RDK",
                                ("maccs166", "RDK", "substruct", "rdmaccs"))

    if args.maccs166:
        opener = types.get_fingerprint_family("RDKit-MACCS166")()
    elif args.RDK:
        fpSize = args.fpSize
        minPath = args.minPath
        maxPath = args.maxPath
        nBitsPerHash = args.nBitsPerHash
        if fpSize < 1:
            parser.error("--fpSize must be positive")
        if nBitsPerHash < 1:
            parser.error("--nBitsPerHash must be a positive value")
        if minPath < 1:
            parser.error("--minPath must be a positive value")
        if maxPath < minPath:
            parser.error("--minPath must not be greater than --maxPath")

        useHs = args.useHs
        if useHs not in (0, 1):
            parser.error("--useHs parameter must be 0 or 1")

        opener = types.get_fingerprint_family("RDKit-Fingerprint")(
            minPath=minPath,
            maxPath=maxPath,
            fpSize=fpSize,
            nBitsPerHash=nBitsPerHash,
            useHs=useHs)

    elif args.substruct:
        opener = types.get_fingerprint_family("Substruct-RDKit")()
    elif args.rdmaccs:
        opener = types.get_fingerprint_family("RDMACCS-RDKit")()
    try:
        reader = opener.read_structure_fingerprints(args.filename, args.format)
    except (TypeError, IOError), err:
        sys.stderr.write("Cannot read structure fingerprints: %s\n" % err)
        raise SystemExit(1)

    io.write_fps1_output(reader, args.output)

if __name__ == "__main__":
    main()

