import sys
from chemfp import openbabel as ob
from chemfp import argparse, shared


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
                    help=ob._fingerprinter_table["FP2"].description + "(default)")
group.add_argument("--FP3", action="store_true",
                    help=ob._fingerprinter_table["FP3"].description)
group.add_argument("--FP4", action="store_true",
                    help=ob._fingerprinter_table["FP4"].description)

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

    if args.FP2:
        fp_name = "FP2"
    elif args.FP3:
        fp_name = "FP3"
    elif args.FP4:
        fp_name = "FP4"
    else:
        # Default
        fp_name = "FP2"

    fp_info = ob._fingerprinter_table[fp_name]
    calc_fp = fp_info.calc_fp

    reader = ob.read_structures(args.filename, args.format)

    shared.generate_fpsv1_output(dict(num_bits=fp_info.num_bits,
                                      software=ob.SOFTWARE,
                                      source=args.filename,
                                      params=fp_info.fps_params),
                                 reader,
                                 calc_fp,
                                 args.output)
    
if __name__ == "__main__":
    main()
