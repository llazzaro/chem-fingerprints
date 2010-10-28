from __future__ import absolute_import

import sys

from .. import argparse, decoders, sdf_reader, shared


def _check_num_bits(num_bits,  # from the user
                    fp_num_bits, # not None if the fp decoder know it exactly
                    num_bytes, # length of decoded fp in bytes
                    parser):
    """Check that the number of fingerprint bits and bytes match the user input

    Difficulties: some fingerprints have only a byte length, and the user
    doesn't have to specify the input.

    Returns the number of bits, or calls parser.error if there are problems
    """
    if fp_num_bits is not None:
        # The fingerprint knows exactly how many bits it contains
        if num_bits is None:
            # The user hasn't specified, so go with the exact number
            return fp_num_bits

        # If the user gave a value, make sure it matches
        if num_bits != fp_num_bits:
            parser.error(
                ("the first fingerprint has {fp_num_bits} bits which "
                 "is not the same as the --num-bits value of {num_bits}").format(
                    num_bits=num_bits, fp_num_bits=fp_num_bits))
            raise AssertionError("should not get here")
        return fp_num_bits

    # If the number of bits isn't specified, assume it's exactly
    # enough to fill up the fingerprint bytes.
    if num_bits is None:
        return num_bytes * 8

    # The user specified the number of bits. The first fingerprint
    # has a number of bytes. This must be enough to hold the bits,
    # but only up to 7 bits larger.
    if (num_bits+7)//8 != num_bytes:
        parser.error(
            ("The byte length of the first fingerprint is {num_bytes} so --num-bits "
             "must be {min} <= num-bits <= {max}, not {num_bits}").format(
                num_bytes=num_bytes, min=num_bytes*8-7, max=num_bytes*8,
                num_bits=num_bits))
        raise AssertError("should not get here")

    # Accept what the user suggested
    return num_bits

parser = argparse.ArgumentParser(
    description="Extract a fingerprint tag from an SD file and generate FPS fingerprints",
    #epilog=epilog,
    #formatter_class=argparse.RawDescriptionHelpFormatter,
    )

parser.add_argument(
    "filename", nargs="?", help="input SD file (default is stdin)", default=None)

parser.add_argument("--title-tag", metavar="TAG", default=None,
            help="get the record title from TAG instead of the first line of the record")
parser.add_argument("--fp-tag", metavar="TAG", 
                    help="get the fingerprint from tag TAG (required)")

parser.add_argument("--num-bits", metavar="INT", type=int,
                    help="use the first INT bits of the input. Use only when the "
                    "last 1-7 bits of the last byte are not part of the fingerprint. "
                    "Unexpected errors will occur if these bits are not all zero.")

parser.add_argument("-o", "--output", metavar="FILENAME",
                    help="save the fingerprints to FILENAME (default=stdout)")
parser.add_argument("--software", metavar="TEXT",
                    help="use TEXT as the software description")
parser.add_argument("--type", metavar="TEXT",
                    help="use TEXT as the fingerprint type description")

# Do I want "--gzip", "--auto", "--none", "--bzip2", and "--decompress METHOD"?
# Do I want to support encoding of the fps output?
parser.add_argument(
    "--decompress", action="store", metavar="METHOD", default="auto",
    help="use METHOD to decompress the input (default='auto', 'none', 'gzip', 'bzip2')")
#parser.add_argument(
#    "--compress", action="store", metavar="METHOD", default="auto",
#    help="use METHOD to compress the output (default='auto', 'none', 'gzip', 'bzip2')")


# This adds --cactvs, --base64 and other decoders to the command-line arguments
decoders._add_decoding_group(parser)

# Support the "--pubchem" option
shortcuts_group = parser.add_argument_group("shortcuts")

class AddSubsKeys(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        namespace.cactvs=True
        namespace.software="PubChem-SubsKeys/1.3",
        namespace.fp_tag="PUBCHEM_CACTVS_SUBSKEYS"

shortcuts_group.add_argument("--pubchem", nargs=0, action=AddSubsKeys,
   help = ("decode CACTVS substructure keys used in PubChem. Same as "
           " --software=PubChem-SubsKeys/1.3 --fp-tag=PUBCHEM_CACTVS_SUBSKEYS --cactvs"))

###############

def main(args=None):
    args = parser.parse_args(args)

    if not args.fp_tag:
        parser.error("argument --fp-tag is required")

    fp_decoder_name, fp_decoder = decoders._extract_decoder(parser, args)

    location = sdf_reader.FileLocation()
    records = sdf_reader.open_sdf(args.filename, args.decompress, loc=location)
    
    if args.title_tag is not None:
        reader = sdf_reader.iter_two_tags(records, args.title_tag, args.fp_tag)
        MISSING_TITLE = "Missing title tag {tag}, ".format(tag=args.title_tag)
        MISSING_TITLE += "line {loc.lineno}. Skipping.\n"
        
    else:
        reader = sdf_reader.iter_title_and_tag(records, args.fp_tag)
        MISSING_TITLE = "Empty record title at line {loc.lineno}. Skipping.\n"

    MISSING_FP = ("Missing fingerprint tag {tag} in record {loc.title!r} line {loc.lineno}. "
                  "Skipping.\n")

    num_bits = args.num_bits

    # I need to get some information from the first record
    first_time = True
    outfile = None       # Don't open it until I'm ready to write the first record
    num_bytes = None     # Will need to get (or at least check) the fingerprint byte length
    expected_fp_num_bits = -1   # 

    def skip(skip_count=[0]):
        if first_time:
            if skip_count[0] > 100:
                raise SystemExit(
                    "ERROR: No fingerprints found in the first 100 records. Exiting.")
            skip_count[0] += 1

    for title, encoded_fp in reader:
        if not title:
            sys.stderr.write(MISSING_TITLE.format(loc=location))
            skip()
            continue
        if not encoded_fp:
            sys.stderr.write(MISSING_FP.format(tag=args.fp_tag, loc=location))
            skip()
            continue
        try:
            fp_num_bits, fp = fp_decoder(encoded_fp)
        except TypeError, err:
            sys.stderr.write(
                ("Could not {decoder_name} decode {tag} value {encoded_fp!r}: {err}\n"
                 "  Skipping record {message}\n").format(
                    decoder_name=fp_decoder_name, tag=args.fp_tag,
                    message=location.message(), err=err, encoded_fp=encoded_fp))
            skip()
            continue
        
        if first_time:
            first_time = False
            num_bytes = len(fp)
            num_bits = _check_num_bits(num_bits, fp_num_bits, num_bytes, parser)
            expected_fp_num_bits = fp_num_bits
            expected_num_bytes = num_bytes

            # Now I know num_bits and num_bytes
            # Time to create output!
            outfile = shared.open_output(args.output)
            shared.write_to_pipe(outfile,
                                 shared.format_fpsv1_header(
                    num_bits=num_bits,
                    software=None, # args.software XXX
                    params=None, # args.params, XXX
                    source=args.filename))

        else:
            if (fp_num_bits != expected_fp_num_bits or
                len(fp) != expected_num_bytes):
                raise SystemExit(
                    ("ERROR: The {message}, tag {tag} has an inconsistent "
                     "fingerprint length".format(
                            message=location.message(), tag=args.fp_tag)))
            
        shared.write_to_pipe(outfile,
                             "%s %s\n" % (fp.encode("hex"), title))

if __name__ == "__main__":
    main()
