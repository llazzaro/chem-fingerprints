from __future__ import absolute_import
import sys
import textwrap

from .. import argparse, types, io
from .. import openeye as oe

##### Handle command-line argument parsing

# Build up some help text based on the atype and btype fields
atype_options = "\n  ".join(textwrap.wrap(" ".join(sorted(oe._atypes))))
btype_options = "\n  ".join(textwrap.wrap(" ".join(sorted(oe._btypes))))

# Extra help text after the parameter descriptions
epilog = """\

ATYPE is one or more of the following, separated by commas
  %(atype_options)s
Examples:
  --atype DefaultAtom
  --atype AtomicNumber,HvyDegree

BTYPE is one or more of the following, separated by commas
  %(btype_options)s
Examples:
  --btype DefaultBond,Chiral
  --btype BondOrder

OEChem guesses the input structure format based on the filename
extension and assumes SMILES for structures read from stdin.
Use "--in FORMAT" to select an alternative, where FORMAT is one of:
 
  File Type      Valid FORMATs (use gz if compressed)
  ---------      ------------------------------------
   SMILES        smi, ism, can, smi.gz, ism.gz, can.gz
   SDF           sdf, mol, sdf.gz, mol.gz
   SKC           skc, skc.gz
   CDK           cdk, cdk.gz
   MOL2          mol2, mol2.gz
   PDB           pdb, ent, pdb.gz, ent.gz
   MacroModel    mmod, mmod.gz
   OEBinary v2   oeb, oeb.gz
   old OEBinary  bin
""" % dict(atype_options=atype_options,
           btype_options=btype_options)

parser = argparse.ArgumentParser(
    description="Generate FPS fingerprints from a structure file using OEChem",
    epilog=epilog,
    formatter_class=argparse.RawDescriptionHelpFormatter,    
    )
path_group = parser.add_argument_group("path fingerprints")
path_group.add_argument(
    "--path", action="store_true", help="generate path fingerprints (default)")
path_group.add_argument(
    "--numbits", type=int, metavar="INT",
    help="number of bits in the path fingerprint (default=4096)", default=4096)
path_group.add_argument(
    "--minbonds", type=int, metavar="INT",
    help="minimum number of bonds in path (default=0)", default=0)
path_group.add_argument(
    "--maxbonds", type=int, metavar="INT",
    help="maximum number of bonds in path (default=5)", default=5)
path_group.add_argument(
    "--atype", help="atom type (default=DefaultAtom)", default="DefaultAtom")
path_group.add_argument(
    "--btype", help="bond type (default=DefaultBond)", default="DefaultBond")

maccs_group = parser.add_argument_group("166 bit MACCS substructure keys")
maccs_group.add_argument(
    "--maccs166", action="store_true", help="generate MACCS fingerprints")

substruct_group = parser.add_argument_group("881 bit substructure keys")
substruct_group.add_argument(
    "--substruct", action="store_true", help="generate ChemFP substructure fingerprints")
substruct_group.add_argument(
    "--hydrogens", action="store", help="one of 'none', 'no-implicit', 'all'",
    default="all")

rdmaccs_group = parser.add_argument_group("RDKit version of the 166 bit MACCS keys")
substruct_group.add_argument(
    "--rdmaccs", action="store_true", help="generate ChemFP RDKit/MACCS")
    
parser.add_argument(
    "--in", metavar="FORMAT", dest="format",
    help="input structure format (default guesses from filename)")
parser.add_argument(
    "-o", "--output", metavar="FILENAME",
    help="save the fingerprints to FILENAME (default=stdout)")
parser.add_argument(
    "filename", nargs="?", help="input structure file (default is stdin)", default=None)


#######

def main(args=None):
    args = parser.parse_args(args)
    outfile = sys.stdout

    groups = ("maccs166", "path", "substruct", "rdmaccs")
    for i, g1 in enumerate(groups[:-1]):
        if not getattr(args, g1):
            continue
        for g2 in groups[i+1:]:
            if getattr(args, g1):
                parser.error("Cannot specify both --%s and --%s" % (g1, g2))

    if args.maccs166:
        # Create the MACCS keys fingerprinter
        opener = types.OpenEyeMACCS166()
    elif args.path:
        if not (16 <= args.numbits <= 65536):
            parser.error("--numbits must be between 16 and 65536 bits")

        if not (0 <= args.minbonds):
            parser.error("--minbonds must be 0 or greater")
        if not (args.minbonds <= args.maxbonds):
            parser.error("--maxbonds must not be smaller than --minbonds")

        # Parse the arguments
        try:
            atype = oe.atom_description_to_value(args.atype)
            btype = oe.bond_description_to_value(args.btype)
        except TypeError, err:
            parser.error(str(err))


        opener = types.OpenEyePath({"numbits": args.numbits,
                                    "minbonds": args.minbonds,
                                    "maxbonds": args.maxbonds,
                                    "atype": atype,
                                    "btype": btype})
    elif args.substruct:
        try:
            types.check_hydrogens(args.hydrogens)
        except TypeError:
            parser.error("--hydrogens must be one of 'none', 'no-implicit' or 'all'; not %r" %
                         (args.hydrogens,))
        opener = types.ChemFPOESubstruct({"hydrogens": args.hydrogens})

    elif args.rdmaccs:
        opener = types.ChemFPRDMACCSOE({"hydrogens": "all"})
        
    else:
        parser.error("???")

    # Ready the input reader/iterator
    try:
        reader = opener.read_structure_fingerprints(args.filename, args.format)
    except (IOError, oe.UnknownFormat), err:
        sys.stderr.write("Cannot read structure fingerprints: %s\n" % err)
        raise SystemExit(1)

    io.write_fps1_output(reader, args.output)

if __name__ == "__main__":
    main()
