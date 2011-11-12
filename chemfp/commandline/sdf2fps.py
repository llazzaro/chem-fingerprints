from __future__ import absolute_import

import sys
import re
import itertools

from .. import Metadata, FingerprintIterator, ParseError
from .. import argparse
from .. import decoders
from .. import sdf_reader
from .. import io
from .. import error_handlers

from . import cmdsupport

# Backwards compatibility support for Python 2.5
try:
    next
except NameError:
    def next(it):
        return it.next()
    

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
                ("the first fingerprint has %(fp_num_bits)s bits which "
                 "is not the same as the --num-bits value of %(num_bits)s") % dict(
                    num_bits=num_bits, fp_num_bits=fp_num_bits))
            raise AssertionError("should not get here")
        
        return num_bits

    # If the number of bits isn't specified, assume it's exactly
    # enough to fill up the fingerprint bytes.
    if num_bits is None:
        return num_bytes * 8

    # The user specified the number of bits. The first fingerprint
    # has a number of bytes. This must be enough to hold the bits,
    # but only up to 7 bits larger.
    if (num_bits+7)//8 != num_bytes:
        parser.error(
            ("The byte length of the first fingerprint is %(num_bytes)s so --num-bits "
             "must be %(min)s <= num-bits <= %(max)s, not %(num_bits)s") % dict(
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
    "filenames", nargs="*", help="input SD files (default is stdin)", default=None)

parser.add_argument("--id-tag", metavar="TAG", default=None,
            help="get the record id from TAG instead of the first line of the record")
parser.add_argument("--fp-tag", metavar="TAG", 
                    help="get the fingerprint from tag TAG (required)")

parser.add_argument("--num-bits", metavar="INT", type=int,
                    help="use the first INT bits of the input. Use only when the "
                    "last 1-7 bits of the last byte are not part of the fingerprint. "
                    "Unexpected errors will occur if these bits are not all zero.")

parser.add_argument(
    "--errors", choices=["strict", "report", "ignore"], default="strict",
    help="how should structure parse errors be handled? (default=strict)")

parser.add_argument("-o", "--output", metavar="FILENAME",
                    help="save the fingerprints to FILENAME (default=stdout)")
parser.add_argument("--software", metavar="TEXT",
                    help="use TEXT as the software description")
parser.add_argument("--type", metavar="TEXT",
                    help="use TEXT as the fingerprint type description")

# TODO:
# Do I want "--gzip", "--auto", "--none", "--bzip2", and "--decompress METHOD"?
# Do I want to support encoding of the fps output?
# Or, why support all these? Why not just "--in gz", "--in bz2" and be done
#  with it (do I really need to specify the 'auto' and 'none' options?)
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
        # the 1.3 is solely based on the version of the document at
        #  ftp://ftp.ncbi.nlm.nih.gov/pubchem/specifications/pubchem_fingerprints.txt
        namespace.software="CACTVS/unknown"
        namespace.type="CACTVS-E_SCREEN/1.0 extended=2"
        namespace.fp_tag="PUBCHEM_CACTVS_SUBSKEYS"

shortcuts_group.add_argument("--pubchem", nargs=0, action=AddSubsKeys,
   help = ("decode CACTVS substructure keys used in PubChem. Same as "
           "--software=CACTVS/unknown --type 'CACTVS-E_SCREEN/1.0 extended=2' "
           "--fp-tag=PUBCHEM_CACTVS_SUBSKEYS --cactvs"))

###############

_illegal_value_pat = re.compile(r"[\000-\037]")

def main(args=None):
    args = parser.parse_args(args)

    if not args.fp_tag:
        parser.error("argument --fp-tag is required")
    if args.num_bits is not None and args.num_bits <= 0:
        parser.error("--num-bits must be a positive integer")

    fp_decoder_name, fp_decoder = decoders._extract_decoder(parser, args)

    missing = cmdsupport.check_filenames(args.filenames)
    if missing:
        parser.error("Structure file %r does not exist" % (missing,))

    for attr in ("software", "type"):
        description = getattr(args, attr, None)
        if description is None:
            continue
        m = _illegal_value_pat.search(description)
        if m is None:
            continue
        parser.error("--%(attr)s description may not contain the character %(c)r" % dict(
                attr=attr, c = m.group(0)))

    error_handler = error_handlers.get_parse_error_handler(args.errors)

    # What follows is a bit tricky. I set up a chain of iterators:
    #   - iterate through the SDF iterators
    #   -   iterate through the (id, encoded_fp) pairs in each SDF iterator
    #   -     convert to (id, fp, num_bits) 3-element tuples
    #   -       use the first element to figure out the right metadata
    #   -       send to (id, fp) information to the io.write_fps1_output function


    # Iterate through each of the filenames, yielding the corresponding SDF iterator
    location = sdf_reader.FileLocation()
    def get_sdf_iters():
        if not args.filenames:
            yield sdf_reader.open_sdf(None, args.decompress, location=location)
        else:
            for filename in args.filenames:
                location.filename = filename
                location.lineno = 1
                yield sdf_reader.open_sdf(filename, args.decompress, location=location)

    # Set up the error messages for missing id or fingerprints.
    if args.id_tag is None:
        MISSING_ID = "Missing title in the record starting %(where)s"
        MISSING_FP = "Missing fingerprint tag %(tag)r in record starting %(where)s"
    else:
        MISSING_ID = "Missing id tag %(tag)r in the record starting %(where)s"
        MISSING_FP = "Missing fingerprint tag %(tag)r in record %(id)r starting %(where)s"

    # For each SDF iterator, yield the (id, encoded_fp) pairs
    if args.id_tag is None:
        def iter_encoded_fingerprints(sdf_iters):
            counter = itertools.count(1)
            for sdf_iter in sdf_iters:
                for id, fp in sdf_reader.iter_title_and_tag(sdf_iter, args.fp_tag):
                    if id:
                       id = io.remove_special_characters_from_id(id)
                    yield id, fp
    else:
        def iter_encoded_fingerprints(sdf_iters):
            counter = itertools.count(1)
            for sdf_iter in sdf_iters:
                for id, fp in sdf_reader.iter_two_tags(sdf_iter, args.id_tag, args.fp_tag):
                    if id:
                       id = io.remove_special_characters_from_id(id)
                    yield id, fp


    # This is either None or a user-specified integer
    num_bits = args.num_bits

    # At this point I don't have enough information to generate the metadata.
    # I won't get that until I've read the first record.
    outfile = None       # Don't open it until I'm ready to write the first record
    num_bytes = None     # Will need to get (or at least check) the fingerprint byte length

    # Decoded encoded fingerprints, yielding (id, fp, num_bits)
    
    def decode_fingerprints(encoded_fp_reader, error_handler):
        expected_num_bits = -1
        expected_fp_size = None
        
        for id, encoded_fp in encoded_fp_reader:
            if not id:
                msg = MISSING_ID % dict(id=id, where=location.where(),
                                        tag=args.id_tag)
                error_handler(msg)
                continue
            
            if not encoded_fp:
                msg = MISSING_FP % dict(id=id, where=location.where(),
                                        tag=args.fp_tag)
                error_handler(msg)
                continue

            # Decode the fingerprint, and complain if it isn't decodeable.
            try:
                num_bits, fp = fp_decoder(encoded_fp)
            except TypeError, err:
                msg = ("Could not %(decoder_name)s decode %(tag)r value %(encoded_fp)r: %(err)s %(where)s" %
                       dict(decoder_name=fp_decoder_name, tag=args.fp_tag,
                            where=location.where(), err=err, encoded_fp=encoded_fp))
                error_handler(msg)
                continue

            if num_bits != expected_num_bits:
                if expected_num_bits == -1:
                    expected_num_bits = num_bits
                else:
                    msg = ("Tag %(tag)r value %(encoded_fp)r has %(got)d bits but expected %(expected)d %(where)s" %
                           dict(tag=args.fp_tag, encoded_fp=encoded_fp,
                                got=num_bits, expected=expected_num_bits,
                                where=location.where()))
                    error_handler(msg)
                    continue

            if len(fp) != expected_fp_size:
                if expected_fp_size is None:
                    expected_fp_size = len(fp)
                else:
                    msg = ("Tag %(tag)r value %(encoded_fp)r has %(got)d bytes but expected %(expected)d %(where)s" %
                           dict(tag=args.fp_tag, encoded_fp=encoded_fp,
                                got=len(fp), expected=expected_fp_size,
                                where=location.where()))
                    error_handler(msg)
                    continue

            yield id, fp, num_bits



    sdf_iters = get_sdf_iters()
    encoded_fps = iter_encoded_fingerprints(sdf_iters)
    decoded_fps = decode_fingerprints(encoded_fps, error_handler)

    try:
        id, fp, num_bits = next(decoded_fps)
    except ParseError, err:
        sys.stderr.write("ERROR: %s. Exiting." % (err,))
        raise SystemExit(1)
    except StopIteration:
        # No fingerprints? Make a new empty stream
        metadata = Metadata(date = io.utcnow())
        chained_reader = iter([])

    else:
        # Got the first fingerprint
        expected_num_bytes = len(fp)

        # Verify that they match
        expected_num_bits = _check_num_bits(args.num_bits, num_bits, expected_num_bytes, parser)
        

        chained_reader = itertools.chain( [(id, fp)], (x[:2] for x in decoded_fps) )
        metadata = Metadata(num_bits = expected_num_bits,
                            software = args.software,
                            type = args.type,
                            sources = args.filenames,
                            date = io.utcnow())

    try:
        io.write_fps1_output(chained_reader, args.output, metadata)
    except ParseError, err:
        sys.stderr.write("ERROR: %s. Exiting." % (err,))
        raise SystemExit(1)

if __name__ == "__main__":
    main()
