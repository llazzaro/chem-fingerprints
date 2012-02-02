from __future__ import absolute_import
import sys
import itertools
import textwrap

from .. import ParseError
from .. import argparse, types, io
from .. import openeye as oe
from . import cmdsupport

##### Handle command-line argument parsing

# Build up some help text based on the atype and btype fields
atype_options = "\n  ".join(textwrap.wrap(" ".join(sorted(dict(oe._atype_flags)))))
btype_options = "\n  ".join(textwrap.wrap(" ".join(sorted(dict(oe._btype_flags)))))
if oe.OEGRAPHSIM_API_VERSION == "1":
    from openeye.oegraphsim import (OEGetFPAtomType, OEFPAtomType_DefaultAtom, OEFPAtomType_DefaultAtom,
                                    OEGetFPBondType, OEFPBondType_DefaultBond, OEFPBondType_DefaultBond)
    type_help = """\
ATYPE is one or more of the following, separated by the '|' character.
  %(atype_options)s
The terms 'Default' and 'DefaultAtom' are expanded to OpenEye's
suggested default of %(defaultatom)s.
Examples:
  --atype Default
  --atype AtomicNumber|HvyDegree
(Note that most atom type names change in OEGraphSim 2.0.0.)

BTYPE is one or more of the following, separated by the '|' character
  %(btype_options)s
The terms 'Default' and 'DefaultBond' are expanded to OpenEye's
suggested default of %(defaultbond)s.
Examples:
  --btype Default
  --btype BondOrder
(Note that "BondOrder" changes to "Order" in OEGraphSim 2.0.0.)

For simpler Unix command-line compatibility, a comma may be used
instead of a '|' to separate different fields. Example:
  --atype AtomicNumber,HvyDegree
""" % dict(atype_options=atype_options,
           btype_options=btype_options,
           defaultatom = OEGetFPAtomType(OEFPAtomType_DefaultAtom),
           defaultbond = OEGetFPBondType(OEFPBondType_DefaultBond))
else:
    from openeye.oegraphsim import (
        OEGetFPAtomType, OEFPAtomType_DefaultPathAtom,
        OEFPAtomType_DefaultCircularAtom, OEFPAtomType_DefaultTreeAtom,
        OEGetFPBondType, OEFPBondType_DefaultPathBond,
        OEFPBondType_DefaultCircularBond, OEFPBondType_DefaultTreeBond,
        )
    type_help = """\
ATYPE is one or more of the following, separated by the '|' character
  %(atype_options)s
The following shorthand terms and expansions are also available:
 DefaultPathAtom = %(defaultpathatom)s
 DefaultCircularAtom = %(defaultcircularatom)s
 DefaultTreeAtom = %(defaulttreeatom)s
and 'Default' selects the correct value for the specified fingerprint.
Examples:
  --atype Default
  --atype Arom|AtmNum|FCharge|HCount

BTYPE is one or more of the following, separated by the '|' character
  %(btype_options)s
The following shorthand terms and expansions are also available:
 DefaultPathBond = %(defaultpathbond)s
 DefaultCircularBond = %(defaultcircularbond)s
 DefaultTreeBond = %(defaulttreebond)s
and 'Default' selects the correct value for the specified fingerprint.
Examples:
   --btype Default
   --btype Order|InRing

To simplify command-line use, a comma may be used instead of a '|' to
separate different fields. Example:
  --atype AtmNum,HvyDegree
""" % dict(atype_options=atype_options,
           btype_options=btype_options,
           defaultpathatom=OEGetFPAtomType(OEFPAtomType_DefaultPathAtom),
           defaultcircularatom=OEGetFPAtomType(OEFPAtomType_DefaultCircularAtom),
           defaulttreeatom=OEGetFPAtomType(OEFPAtomType_DefaultTreeAtom),
           defaultpathbond=OEGetFPBondType(OEFPBondType_DefaultPathBond),
           defaultcircularbond=OEGetFPBondType(OEFPBondType_DefaultCircularBond),
           defaulttreebond=OEGetFPBondType(OEFPBondType_DefaultTreeBond),

)


# Extra help text after the parameter descriptions
epilog = type_help + """\

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
"""

parser = argparse.ArgumentParser(
    description="Generate FPS fingerprints from a structure file using OEChem",
    epilog=epilog,
    formatter_class=argparse.RawDescriptionHelpFormatter,    
    )
if oe.OEGRAPHSIM_API_VERSION == "1":
    PathFamily = oe.OpenEyePathFingerprintFamily_v1
    path_group = parser.add_argument_group("path fingerprints")
    path_group.add_argument(
        "--path", action="store_true", help="generate path fingerprints (default)")
    PathFamily.add_argument_to_argparse("numbits", path_group)
    PathFamily.add_argument_to_argparse("minbonds", path_group)
    PathFamily.add_argument_to_argparse("maxbonds", path_group)
else:
    CircularFamily = oe.OpenEyeCircularFingerprintFamily_v2
    path_group = parser.add_argument_group("path, circular, and tree fingerprints")
    path_group.add_argument(
        "--path", action="store_true", help="generate path fingerprints (default)")
    path_group.add_argument(
        "--circular", action="store_true", help="generate circular fingerprints")
    path_group.add_argument(
        "--tree", action="store_true", help="generate tree fingerprints")

    path_group.add_argument(
        "--numbits", action="store", type=int, metavar="INT", default=4096,
        help="number of bits in the fingerprint (default=4096)")
    path_group.add_argument(
        "--minbonds", action="store", type=int, metavar="INT", default=0,
        help="minimum number of bonds in the path or tree fingerprint (default=0)")
    path_group.add_argument(
        "--maxbonds", action="store", type=int, metavar="INT", default=None,
        help="maximum number of bonds in the path or tree fingerprint (path default=5, tree default=4)")
    CircularFamily.add_argument_to_argparse("minradius", path_group)
    CircularFamily.add_argument_to_argparse("maxradius", path_group)

# The expansion of 'Default' differs based on the fingerprint type
path_group.add_argument(
    "--atype", metavar="ATYPE", default="Default",
    help="atom type flags, described below (default=Default)")
path_group.add_argument(
    "--btype", metavar="BTYPE", default="Default",
    help="bond type flags, described below (default=Default)")

maccs_group = parser.add_argument_group("166 bit MACCS substructure keys")
maccs_group.add_argument(
    "--maccs166", action="store_true", help="generate MACCS fingerprints")

substruct_group = parser.add_argument_group("881 bit ChemFP substructure keys")
substruct_group.add_argument(
    "--substruct", action="store_true", help="generate ChemFP substructure fingerprints")

rdmaccs_group = parser.add_argument_group("ChemFP version of the 166 bit RDKit/MACCS keys")
rdmaccs_group.add_argument(
    "--rdmaccs", action="store_true", help="generate 166 bit RDKit/MACCS fingerprints")

parser.add_argument(
    "--aromaticity", metavar="NAME", choices=oe._aromaticity_flavor_names,
    default="openeye",
    help="use the named aromaticity model")

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

def _get_atype_and_btype(args, atom_description_to_value, bond_description_to_value):
    try:
        atype = atom_description_to_value(args.atype)
    except ValueError, err:
        parser.error("--atype must contain '|' separated atom terms: %s" % (err,))
    try:
        btype = bond_description_to_value(args.btype)
    except ValueError, err:
        parser.error("--btype must contain '|' separated atom terms: %s" % (err,))
    return atype, btype

#######

def main(args=None):
    args = parser.parse_args(args)

    supported_fingerprints = ("maccs166", "path", "substruct", "rdmaccs")
    if oe.OEGRAPHSIM_API_VERSION != "1":
        supported_fingerprints += ("circular", "tree")
    cmdsupport.mutual_exclusion(parser, args, "path", supported_fingerprints)

    if args.maccs166:
        # Create the MACCS keys fingerprinter
        opener = types.get_fingerprint_family("OpenEye-MACCS166")()
    elif args.path:
        if not (16 <= args.numbits <= 65536):
            parser.error("--numbits must be between 16 and 65536 bits")

        if not (0 <= args.minbonds):
            parser.error("--minbonds must be 0 or greater")
        if args.maxbonds is None:
            args.maxbonds = 5
        if not (args.minbonds <= args.maxbonds):
            parser.error("--maxbonds must not be smaller than --minbonds")
        atype, btype = _get_atype_and_btype(args, oe.path_atom_description_to_value,
                                            oe.path_bond_description_to_value)
        opener = types.get_fingerprint_family("OpenEye-Path")(
            numbits = args.numbits,
            minbonds = args.minbonds,
            maxbonds = args.maxbonds,
            atype = atype,
            btype = btype)
    elif args.circular:
        if not (16 <= args.numbits <= 65536):
            parser.error("--numbits must be between 16 and 65536 bits")

        if not (0 <= args.minradius):
            parser.error("--minradius must be 0 or greater")
        if not (args.minradius <= args.maxradius):
            parser.error("--maxradius must not be smaller than --minradius")
        atype, btype = _get_atype_and_btype(args, oe.circular_atom_description_to_value,
                                            oe.circular_bond_description_to_value)

        opener = types.get_fingerprint_family("OpenEye-Circular")(
            numbits = args.numbits,
            minradius = args.minradius,
            maxradius = args.maxradius,
            atype = atype,
            btype = btype)
    elif args.tree:
        if not (16 <= args.numbits <= 65536):
            parser.error("--numbits must be between 16 and 65536 bits")

        if not (0 <= args.minbonds):
            parser.error("--minbonds must be 0 or greater")
        if args.maxbonds is None:
            args.maxbonds = 4
        if not (args.minbonds <= args.maxbonds):
            parser.error("--maxbonds must not be smaller than --minbonds")
        atype, btype = _get_atype_and_btype(args, oe.tree_atom_description_to_value,
                                            oe.tree_bond_description_to_value)

        opener = types.get_fingerprint_family("OpenEye-Tree")(
            numbits = args.numbits,
            minbonds = args.minbonds,
            maxbonds = args.maxbonds,
            atype = atype,
            btype = btype)
        
    elif args.substruct:
        opener = types.get_fingerprint_family("ChemFP-Substruct-OpenEye")()
    elif args.rdmaccs:
        opener = types.get_fingerprint_family("RDMACCS-OpenEye")()
        
    else:
        parser.error("ERROR: fingerprint not specified?")

    if args.format is not None:
        if args.filenames:
            filename = args.filenames[0]
        else:
            filename = None
        if not oe.is_valid_format(filename, args.format):
            parser.error("Unsupported format specifier: %r" % (args.format,))

    if not oe.is_valid_aromaticity(args.aromaticity):
        parser.error("Unsupported aromaticity specifier: %r" % (args.aromaticity,))

    if not cmdsupport.is_valid_tag(args.id_tag):
        parser.error("Invalid id tag: %r" % (args.id_tag,))

    missing = cmdsupport.check_filenames(args.filenames)
    if missing:
        parser.error("Structure file %r does not exist" % (missing,))

    # Ready the input reader/iterator
    metadata, reader = cmdsupport.read_multifile_structure_fingerprints(
        opener, args.filenames, args.format, args.id_tag, args.aromaticity, args.errors)
    
    try:
        io.write_fps1_output(reader, args.output, metadata)
    except ParseError, err:
        sys.stderr.write("ERROR: %s. Exiting.\n" % (err,))
        raise SystemExit(1)
    
if __name__ == "__main__":
    main()
