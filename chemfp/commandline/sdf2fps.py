from __future__ import absolute_import

import sys

from chemfp import argparse, decoders, sdf_reader, shared


def _check_num_bits(num_bits,  # from the user
                    fp_num_bits, # not None if the fp decoder know it exactly
                    num_bytes, # length of decoded fp in bytes
                    parser):
    if fp_num_bits is not None:
        # The fingerprint knows exactly how many bits it contains
        if num_bits is None:
            return fp_num_bits

        if num_bits != fp_num_bits:
            parser.error(
                ("--num-bits value of {num_bits} does not match "
                 "the first fingerprint size of {fp_num_bits}").format(
                    num_bits=num_bits, fp_num_bits=fp_num_bits))
        raise AssertionError("should not get here")

    # The fingerprint size may be rounded up to the nearest
    # factor of 8 from the expected size.
    if num_bits is None:
        return num_bytes * 8

    if (num_bits+7)//8 != num_bytes:
        parser.error(
            ("The first fingerprint has {num_bytes} bytes so --num-bits "
             "must be between {min} <= x <= {max}, not {num_bits}").format(
                num_bytes=num_bytes, min=num_bytes*8-7, max=num_bytes*8,
                num_bits=num_bits))
        raise AssertError("should not get here")

    # Accept what the user suggested
    return num_bits

parser = argparse.ArgumentParser(
    description="Generate FPS fingerprints from tags in an SD file",
    #epilog=epilog,
    #formatter_class=argparse.RawDescriptionHelpFormatter,
    )

parser.add_argument(
    "filename", nargs="?", help="input SD file (default is stdin)", default=None)

parser.add_argument("--fp-tag", metavar="TAG", 
                    help="tag TAG contains the fingerprint (required)")

parser.add_argument("--num-bits", metavar="INT",
                    help="use the first INT bits of the input")
parser.add_argument("-o", "--output", metavar="FILENAME",
                    help="save the fingerprints to FILENAME (default=stdout)")
parser.add_argument("--software", metavar="TEXT",
                    help="use TEXT as the software description")
parser.add_argument("--title-tag", metavar="TAG", default=None,
                    help="tag TAG contains the title, instead of the record title")
parser.add_argument("--type", metavar="TEXT",
                    help="use TEXT as the fingerprint type description")

parser.add_argument(
    "-d", "--decompress", action="store", metavar="METHOD", default="auto",
    help="use METHOD to decompress the input (default='auto', 'none', 'gzip', 'bzip2')")

# This adds --cactvs, --base64 and other decoders to the command-line arguments
decoders._add_decoding_group(parser)

shortcuts_group = parser.add_argument_group("shortcuts")

class AddSubsKeys(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        namespace.cactvs=True
        namespace.software="PubChem-SubsKeys/1.3",
        namespace.fp_tag="PUBCHEM_CACTVS_SUBSKEYS"

shortcuts_group.add_argument("--pubchem-subskeys", nargs=0, action=AddSubsKeys,
   help = ("decode CACTVS substructure keys. Same as "
           " --software=PubChem-SubsKeys/1.3 --fp-tag=PUBCHEM_CACTVS_SUBSKEYS --cactvs"))


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
    outfile = None
    num_bytes = None
    expected_fp_num_bits = -1

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
