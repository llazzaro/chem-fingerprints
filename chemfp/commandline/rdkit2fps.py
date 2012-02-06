# Copyright (c) 2010 Andrew Dalke Scientific, AB (Gothenburg, Sweden)
import sys

from .. import ParseError
from .. import argparse, io, rdkit, types
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
    conflict_handler="resolve",
    )

_base = rdkit._base

# --RDK and --morgan both have fpSize but argparse doesn't allow the
# same option in different groups. Especially with different defaults.

_base.add_argument_to_argparse("fpSize", parser)

rdk_group = parser.add_argument_group("RDKit topological fingerprints")
rdk_group.add_argument("--RDK", action="store_true",
                       help="generate RDK fingerprints (default)")
_base.add_argument_to_argparse("minPath", rdk_group)
_base.add_argument_to_argparse("maxPath", rdk_group)
_base.add_argument_to_argparse("nBitsPerHash", rdk_group)
_base.add_argument_to_argparse("useHs", rdk_group)


morgan_group = parser.add_argument_group("RDKit Morgan fingerprints")

morgan_group.add_argument("--morgan", action="store_true",
                          help="generate Morgan fingerprints")

_morgan = rdkit.RDKitMorganFingerprintFamily_v1
_morgan.add_argument_to_argparse("radius", morgan_group)
_morgan.add_argument_to_argparse("useFeatures", morgan_group)
_morgan.add_argument_to_argparse("useChirality", morgan_group)
_morgan.add_argument_to_argparse("useBondTypes", morgan_group)

torsion_group = parser.add_argument_group("RDKit Topological Torsion fingerprints")
torsion_group.add_argument("--torsions", action="store_true",
                           help="generate Topological Torsion fingerprints")
rdkit.RDKitTorsionFingerprintFamily_v1.add_argument_to_argparse(
    "targetSize", torsion_group)

pair_group = parser.add_argument_group("RDKit Atom Pair fingerprints")
pair_group.add_argument("--pairs", action="store_true",
                        help="generate Atom Pair fingerprints")
rdkit.RDKitTorsionFingerprintFamily_v1.add_argument_to_argparse(
    "minLength", pair_group)
rdkit.RDKitTorsionFingerprintFamily_v1.add_argument_to_argparse(
    "maxLength", pair_group)



maccs_group = parser.add_argument_group("166 bit MACCS substructure keys")
maccs_group.add_argument(
    "--maccs166", action="store_true", help="generate MACCS fingerprints")

substruct_group = parser.add_argument_group("881 bit substructure keys")
substruct_group.add_argument(
    "--substruct", action="store_true", help="generate ChemFP substructure fingerprints")

rdmaccs_group = parser.add_argument_group("ChemFP version of the 166 bit RDKit/MACCS keys")
rdmaccs_group.add_argument(
    "--rdmaccs", action="store_true", help="generate 166 bit RDKit/MACCS fingerprints")

parser.add_argument(
    "--id-tag", metavar="NAME",
    help="tag name containing the record id (SD files only)")
    
parser.add_argument(
    "--in", metavar="FORMAT", dest="format",
    help="input structure format (default guesses from filename)")
parser.add_argument(
    "-o", "--output", metavar="FILENAME",
    help="save the fingerprints to FILENAME (default=stdout)")

parser.add_argument(
    "--errors", choices=["strict", "report", "ignore"], default="strict",
    help="how should structure parse errors be handled? (default=strict)")

parser.add_argument(
    "filenames", nargs="*", help="input structure files (default is stdin)")


def main(args=None):
    args = parser.parse_args(args)

    cmdsupport.mutual_exclusion(parser, args, "RDK",
                                ("maccs166", "RDK", "substruct", "rdmaccs", "morgan",
                                 "torsions", "pairs"))

    if args.maccs166:
        opener = types.get_fingerprint_family("RDKit-MACCS166")()
    elif args.RDK:
        fpSize = args.fpSize or rdkit.NUM_BITS
        minPath = args.minPath
        maxPath = args.maxPath
        nBitsPerHash = args.nBitsPerHash
        if maxPath < minPath:
            parser.error("--minPath must not be greater than --maxPath")

        useHs = args.useHs

        opener = types.get_fingerprint_family("RDKit-Fingerprint")(
            minPath=minPath,
            maxPath=maxPath,
            fpSize=fpSize,
            nBitsPerHash=nBitsPerHash,
            useHs=useHs)

    elif args.substruct:
        opener = types.get_fingerprint_family("ChemFP-Substruct-RDKit")()
    elif args.rdmaccs:
        opener = types.get_fingerprint_family("RDMACCS-RDKit")()
    elif args.morgan:
        opener = types.get_fingerprint_family("RDKit-Morgan")(
            radius=args.radius,
            fpSize=args.fpSize,
            useFeatures=args.useFeatures,
            useChirality=args.useChirality,
            useBondTypes=args.useBondTypes)

    elif args.torsions:
        opener = types.get_fingerprint_family("RDKit-Torsion")(
            fpSize=args.fpSize,
            targetSize=args.targetSize)
    elif args.pairs:
        minLength = args.minLength
        maxLength = args.maxLength
        if maxLength < minLength:
            parser.error("--minLength must not be greater than --maxLength")
        opener = types.get_fingerprint_family("RDKit-Pair")(
            fpSize=args.fpSize,
            minLength=minLength,
            maxLength=maxLength)
        
    else:
        raise AssertionError("Unknown fingerprinter")

    if not rdkit.is_valid_format(args.format):
        parser.error("Unsupported format specifier: %r" % (args.format,))

    if not cmdsupport.is_valid_tag(args.id_tag):
        parser.error("Invalid id tag: %r" % (args.id_tag,))

    missing = cmdsupport.check_filenames(args.filenames)
    if missing:
        parser.error("Structure file %r does not exist" % (missing,))

    metadata, reader = cmdsupport.read_multifile_structure_fingerprints(
        opener, args.filenames, format=args.format,
        id_tag=args.id_tag, aromaticity=None, errors=args.errors)

    try:
        io.write_fps1_output(reader, args.output, metadata)
    except ParseError, err:
        sys.stderr.write("ERROR: %s. Exiting." % (err,))
        raise SystemExit(1)

if __name__ == "__main__":
    main()

