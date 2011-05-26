import sys
from chemfp import openbabel as ob
from chemfp import argparse, io, types

from . import cmdsupport


############ Command-line parser definition

epilog = """\

OpenBabel autodetects the input structure format based on the filename
extension. The default format for structures read from stdin is
SMILES. Use"--in FORMAT" to select an alternative, where FORMAT is
one of the extensions at http://openbabel.org/wiki/List_of_extensions .
For a short list of some common formats:


  File Type      Valid FORMATs
  ---------      -------------
   SMILES        smi, can, smiles
   SDF           sdf, mol, sd, mdl
   MOL2          mol2, ml2
   PDB           pdb, ent
   MacroModel    mmod

If OpenBabel is compiled with zlib support then it will automatically
handle gzip'ed input data if the filename ends with ".gz". You may
optionally include that suffix in the format name.

"""

parser = argparse.ArgumentParser(
    description="Generate FPS fingerprints from a structure file using OpenBabel",
    )
group = parser.add_mutually_exclusive_group()
group.add_argument("--FP2", action="store_true",
#                    help=ob._fingerprinter_table["FP2"].description + "(default)"
                   )
group.add_argument("--FP3", action="store_true",
#                    help=ob._fingerprinter_table["FP3"].description
                   )
group.add_argument("--FP4", action="store_true",
#                    help=ob._fingerprinter_table["FP4"].description
                   )

if ob.HAS_MACCS:
    # Added in OpenBabel 2.3
    group.add_argument("--MACCS", action="store_true",
#                       help=ob._fingerprinter_table["MACCS"].description
                       )
else:
    group.add_argument("--MACCS", action="store_true",
                       help="(Not available using your version of OpenBabel)")

group.add_argument(
    "--substruct", action="store_true", help="generate ChemFP substructure fingerprints")

group.add_argument(
    "--rdmaccs", action="store_true", help="generate ChemFP 166 bit RDKit/MACCS fingerprints")


parser.add_argument(
    "--in", metavar="FORMAT", dest="format",
    help="input structure format (default autodetects from the filename extension)")
parser.add_argument(
    "-o", "--output", metavar="FILENAME",
    help="save the fingerprints to FILENAME (default=stdout)")
parser.add_argument(
    "filename", nargs="?", help="input structure file (default is stdin)", default=None)


#########

def main(args=None):
    args = parser.parse_args(args)
    outfile = sys.stdout

    cmdsupport.mutual_exclusion(parser, args, "FP2",
                                ("FP2", "FP3", "FP4", "MACCS", "substruct", "rdmaccs"))

    if args.FP2:
        opener = types.get_fingerprint_family("OpenBabel-FP2")()
    elif args.FP3:
        opener = types.get_fingerprint_family("OpenBabel-FP3")()
    elif args.FP4:
        opener = types.get_fingerprint_family("OpenBabel-FP4")()
    elif args.MACCS:
        if not ob.HAS_MACCS:
            parser.error(
                "--MACCS is not supported in your OpenBabel installation (%s)" % (
                    ob.GetReleaseVersion(),))
        opener = types.get_fingerprint_family("OpenBabel-MACCS")()
    elif args.substruct:
        opener = types.get_fingerprint_family("ChemFP-Substruct-OpenBabel")()
    elif args.rdmaccs:
        opener = types.get_fingerprint_family("ChemFP-Substruct-OpenBabel")()
    else:
        parser.error("should not get here")

    # Ready the input reader/iterator
    try:
        reader = opener.read_structure_fingerprints(args.filename, args.format)
    except IOError, err:
        sys.stderr.write("Cannot read structures: %s" % (err,))
        raise SystemExit(1)
    except TypeError, err:
        msg = str(err)
        if "Unknown structure format" in msg:
            sys.stderr.write(msg)
            raise SystemExit(1)
        raise

    io.write_fps1_output(reader, args.output)
    
if __name__ == "__main__":
    main()
