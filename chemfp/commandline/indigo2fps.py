from __future__ import absolute_import

import sys
from chemfp import indigo
from chemfp import argparse, io, types

from . import cmdsupport

epilog = """\

indigo2fps autodetects the input structure format based on the filename
extension. The default format for structures read from stdin is
SMILES. Use"--in FORMAT" to select an alternative, where FORMAT is
one of:

  File Type      Valid FORMATs
  ---------      -------------
   SMILES        smi, can, ism, smiles
   SDF           sdf, mol, sd, mdl
   RDF           rdf
   CML           cml

Indigo will automatically handle gzip'ed input data if the filename
ends with".gz". You may optionally include that suffix in the format name.

"""

parser = argparse.ArgumentParser(
    description="Generate FPS fingerprints from a structure file using Indigo",
    )
group = parser.add_mutually_exclusive_group()
group.add_argument("--sim", action="store_true",
                   help="similarity fingerprints")
group.add_argument("--sub", action="store_true",
                   help="substructure fingerprints")
group.add_argument("--sub-res", action="store_true",
                   help="resonance substructure fingerprints")
group.add_argument("--sub-tau", action="store_true",
                   help="tautomer substructure fingerprints")
group.add_argument("--full", action="store_true",
                   help="full fingerprints")

#group.add_argument(
#    "--substruct", action="store_true", help="generate ChemFP substructure fingerprints")
#
#group.add_argument(
#    "--rdmaccs", action="store_true", help="generate 166 bit RDKit/MACCS fingerprints")


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

    cmdsupport.mutual_exclusion(parser, args, "sim",
                                ("sim", "sub", "sub_res", "sub_tau", "full"
                                 # , "substruct", "rdmaccs"
                                 ))
    if args.sim:
        opener = types.get_fingerprint_family("Indigo-Similarity")()
    elif args.sub:
        opener = types.get_fingerprint_family("Indigo-Substructure")()
    elif args.sub_res:
        opener = types.get_fingerprint_family("Indigo-ResonanceSubstructure")()
    elif args.sub_tau:
        opener = types.get_fingerprint_family("Indigo-TautomerSubstructure")()
    elif args.full:
        opener = types.get_fingerprint_family("Indigo-Full")()
#    elif args.substruct:
#        opener = types.get_fingerprint_family("ChemFP-Substruct-Indigo")()
#    elif args.rdmaccs:
#        opener = types.get_fingerprint_family("RDMACCS-Indigo")()
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
