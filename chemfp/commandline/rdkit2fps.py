# Copyright (c) 2010 Andrew Dalke Scientific, AB (Gothenburg, Sweden)

from chemfp import argparse, shared, rdkit

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
    "--num-bits", type=int, metavar="INT", default=rdkit.NUM_BITS,
    help="number of bits in the fingerprint (default=%d)" % rdkit.NUM_BITS)
rdk_group.add_argument(
    "--min-path", type=int, metavar="INT", default=rdkit.MIN_PATH,
    help="minimum number of bonds to include in the subgraphs (default=%d)" % rdkit.MIN_PATH)
rdk_group.add_argument(
    "--max-path", type=int, metavar="INT", default=rdkit.MAX_PATH,
    help="maximum number of bonds to include in the subgraphs (default=%d)" % rdkit.MAX_PATH)
rdk_group.add_argument(
    "--bits-per-hash", type=int, metavar="INT", default=rdkit.BITS_PER_HASH,
    help="number of bits to set per path (default=%d)" % rdkit.BITS_PER_HASH)
rdk_group.add_argument(
    "--ignore-Hs", action="store_true",
    help="do not include information about the number of hydrogens on each atom")

maccs_group = parser.add_argument_group("166 bit MACCS substructure keys")
maccs_group.add_argument(
    "--maccs166", action="store_true", help="generate MACCS fingerprints")

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

    if args.maccs166:
        if args.RDK:
            parser.error("Cannot specify both --maccs166 and --RDK")
        fingerprinter = maccs166_fingerprinter
        num_bits = 166
        params = "RDKit-MACCS166/1"
    else:
        num_bits = args.num_bits
        min_path = args.min_path
        max_path = args.max_path
        bits_per_hash = args.bits_per_hash
        if num_bits < 1:
            parser.error("--num-bits must be positive")
        if bits_per_hash < 1:
            parser.error("--bits-per-hash must be a positive value")
        if min_path < 1:
            parser.error("--min-path must be a positive value")
        if max_path < min_path:
            parser.error("--min-path must not be greater than --max-path")

        use_Hs = not args.ignore_Hs

        fingerprinter = rdkit.make_rdk_fingerprinter(
            num_bits=num_bits, min_path=min_path, max_path=max_path, 
            bits_per_hash = bits_per_hash, use_Hs = use_Hs)
        params = rdkit.format_rdk_params(
            num_bits=num_bits, min_path=min_path, max_path=max_path, 
            bits_per_hash = bits_per_hash, use_Hs = use_Hs)

    reader = rdkit.read_structures(args.filename, args.format)

    shared.generate_fpsv1_output(dict(num_bits=num_bits,
                                      software=rdkit.SOFTWARE,
                                      source=args.filename,
                                      params=params),
                                 reader,
                                 fingerprinter,
                                 args.output)

if __name__ == "__main__":
    main()

