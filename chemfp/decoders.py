"""chemfp.decoders - decode different fingerprint representations into chemfp form

The chemfp fingerprints are stored as byte strings, with the bytes in
least-significant bit order (bit #0 is stored in the first/left-most
byte) and with the bits in most-significant bit order (bit #0 is
stored in the first/right-most bit of the first byte).

Other systems use different encodings. These include:
  - the '0 and '1' characters, as in '00111101'
  - hex encoding, like '3d'
  - base64 encoding, like 'SGVsbG8h'
  - CACTVS's variation of base64 encoding

plus variations of different LSB and MSB orders.

This module decodes most of the fingerprint encodings I have come
across. The fingerprint decoders return a 2-ple of the bit length and
the chemfp fingerprint. The bit length is None unless the bit length
is known exactly, which currently is only the case for the binary and
CACTVS fingerprints. (The hex and other encoders must round the
fingerprints up to a multiple of 8 bits.)

"""
import string
import binascii

_lsb_bit_table = {} # "10000000" -> 1
_msb_bit_table = {} # "00000001" -> 1

_reverse_bits_in_a_byte_transtable = None

# These are in lsb order; 
_lsb_4bit_patterns = (
    "0000", "1000", "0100", "1100",
    "0010", "1010", "0110", "1110",
    "0001", "1001", "0101", "1101",
    "0011", "1011", "0111", "1111")

# Generate '00000000', '10000000', '01000000', ... , '01111111', '11111111'
def _lsb_8bit_patterns():
    for right in _lsb_4bit_patterns:
        for left in _lsb_4bit_patterns:
            yield left + right

def _init():
    to_trans = [None]*256
    for value, bit_pattern in enumerate(_lsb_8bit_patterns()):
        # Each pattern maps to the byte
        byte_value = chr(value)
        to_trans[value] = chr(int(bit_pattern, 2))

        _lsb_bit_table[bit_pattern] = byte_value
        # Include the forms with trailing 0s
        # 10000000, 1000000, 100000, 10000, 1000, 100, 10 and 1 are all 0x01
        # (RDKit fingerprint lengths don't need to be a multiple of 8)
        lsb_pattern = bit_pattern
        while lsb_pattern[-1:] == "0":
            lsb_pattern = lsb_pattern[:-1]
            _lsb_bit_table[lsb_pattern] = byte_value

        msb_pattern = bit_pattern[::-1]
        _msb_bit_table[msb_pattern] = byte_value
        while msb_pattern[:1] == "0":
            msb_pattern = msb_pattern[1:]
            _msb_bit_table[msb_pattern] = byte_value
    global _reverse_bits_in_a_byte_transtable
    _reverse_bits_in_a_byte_transtable = string.maketrans(
        "".join(chr(i) for i in range(256)),
        "".join(to_trans))
    

_init()
assert _lsb_bit_table["10000000"] == "\x01", _lsb_bit_table["10000000"]
assert _lsb_bit_table["1000000"] == "\x01", _lsb_bit_table["1000000"]
assert _lsb_bit_table["100000"] == "\x01"
assert _lsb_bit_table["10000"] == "\x01"
assert _lsb_bit_table["1"] == "\x01"
assert _lsb_bit_table["1111111"] == "\x7f"

assert _msb_bit_table["00000001"] == "\x01"
assert _msb_bit_table["0000001"] == "\x01"
assert _msb_bit_table["000001"] == "\x01"
assert _msb_bit_table["00001"] == "\x01"
assert _msb_bit_table["1"] == "\x01"
assert _msb_bit_table["00000011"] == "\x03"
assert _msb_bit_table["00000011"] == "\x03"
assert _msb_bit_table["10000000"] == "\x80"
assert _msb_bit_table["1000000"] == "\x40"


def from_binary_lsb(text):
    """Convert a string like '00010101' (bit 0 here is off) into '\\xa8'

    The encoding characters '0' and '1' are in LSB order, so bit 0 is the left-most field.
    The result is a 2-ple of the fingerprint length and the decoded chemfp fingerprint

    >>> from_binary_lsb('00010101')
    (8, '\\xa8')
    >>> from_binary_lsb('11101')
    (5, '\\x17')
    >>> from_binary_lsb('00000000000000010000000000000')
    (29, '\\x00\\x80\\x00\\x00')
    >>>
    """
    try:
        N = len(text)
        bytes = []
        for i in range(0, N, 8):
            bytes.append(_lsb_bit_table[text[i:i+8]])
        return (N, "".join(bytes))
    except KeyError:
        raise TypeError("Not a binary string")

def from_binary_msb(text):
    """Convert a string like '10101000' (bit 0 here is off) into '\\xa8'

    The encoding characters '0' and '1' are in MSB order, so bit 0 is the right-most field.

    >>> from_binary_msb('10101000')
    (8, '\\xa8')
    >>> from_binary_msb('00010101')
    (8, '\\x15')
    >>> from_binary_msb('00111')
    (5, '\\x07')
    >>> from_binary_msb('00000000000001000000000000000')
    (29, '\\x00\\x80\\x00\\x00')
    >>>
    """
    # It feels like there should be a faster, more elegant way to do this.
    # While close,
    #   hex(int('00010101', 2))[2:].decode("hex")
    # does not keep the initial 0s
    try:
        N = len(text)
        bytes = []
        end = N
        start = N-8
        while start > 0:
            bytes.append(_msb_bit_table[text[start:end]])
            end = start
            start -= 8
        bytes.append(_msb_bit_table[text[0:end]])
        return (N, "".join(bytes))
    except KeyError:
        raise TypeError("Not a binary string")


def from_base64(text):
    """Decode a base64 encoded fingerprint string

    The encoded fingerprint must be in chemfp form, with the bytes in
    LSB order and the bits in MSB order.

    >>> from_base64("SGk=")
    (None, 'Hi')
    >>> from_base64("SGk=")[1].encode("hex")
    '4869'
    >>> 
    """
    try:
        # This is the same as doing text.decode("base64") but since I
        # need to catch the exception, I might as well work with the
        # underlying implementation code.
        return (None, binascii.a2b_base64(text))
    except binascii.Error, err:
        raise TypeError(str(err))

#def from_base64_msb(text):
#    return (None, text.decode("base64")[::-1], None)

#def from_base64_lsb(text):
#    return (None, text.decode("base64").translate(_reverse_bits_in_a_byte_transtable), None)
    
def from_hex(text):
    """Decode a hex encoded fingerprint string

    The encoded fingerprint must be in chemfp form, with the bytes in
    LSB order and the bits in MSB order.

    >>> from_hex('10f2')
    (None, '\\x10\\xf2')
    >>>

    Raises a TypeError if the hex string is not a multiple of 2 bytes long
    or if it contains a non-hex character.
    """
    return (None, text.decode("hex"))

def from_hex_msb(text):
    """Decode a hex encoded fingerprint string where the bits and bytes are in MSB order

    >>> from_hex_msb('10f2')
    (None, '\\xf2\\x10')
    >>>

    Raises a TypeError if the hex string is not a multiple of 2 bytes long
    or if it contains a non-hex character.
    """
    return (None, text.decode("hex")[::-1])

def from_hex_lsb(text):
    """Decode a hex encoded fingerprint string where the bits and bytes are in LSB order

    >>> from_hex_lsb('102f')
    (None, '\\x08\\xf4')
    >>> 

    Raises a TypeError if the hex string is not a multiple of 2 bytes long
    or if it contains a non-hex character.
    """
    return (None, text.decode("hex").translate(_reverse_bits_in_a_byte_transtable))


# ftp://ftp.ncbi.nlm.nih.gov/pubchem/specifications/pubchem_fingerprints.txt

#   This comes from cid:11 which is 1,2-dichloroethane
# AAADcYBAAAAGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGAIAAAAAAAOAAEAAAAA
# AAAAAAAAAAAAAAAAAAAAAAEAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA==
# That's simple enough to check the bit ordering by eye. Here's the decoded start
#  80-40-00-00-06-00-00 ... 
# We know it has to match the bits (starting with bit 0)
#  1000 0000 0100 0000 0000 0000 0000 0000 0000 0110
# and it does, perfectly. That means CACTVS is pure little endian.
# chem-fp has little-endian byte order but big endian bit order.


# 0111 1000 0100 0000 0000 0101 0000 0000 0000 0000 0000 0000

def from_cactvs(text):
    """Decode a 881-bit CACTVS-encoded fingerprint used by PubChem

    >>> from_cactvs("AAADceB7sQAEAAAAAAAAAAAAAAAAAWAAAAAwAAAAAAAAAAABwAAAHwIYAAAADA" +
    ...             "rBniwygJJqAACqAyVyVACSBAAhhwIa+CC4ZtgIYCLB0/CUpAhgmADIyYcAgAAO" +
    ...             "AAAAAAABAAAAAAAAAAIAAAAAAAAAAA==")
    (881, '\\x07\\xde\\x8d\\x00 \\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x80\\x06\\x00\\x00\\x00\\x0c\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x80\\x03\\x00\\x00\\xf8@\\x18\\x00\\x00\\x000P\\x83y4L\\x01IV\\x00\\x00U\\xc0\\xa4N*\\x00I \\x00\\x84\\xe1@X\\x1f\\x04\\x1df\\x1b\\x10\\x06D\\x83\\xcb\\x0f)%\\x10\\x06\\x19\\x00\\x13\\x93\\xe1\\x00\\x01\\x00p\\x00\\x00\\x00\\x00\\x00\\x80\\x00\\x00\\x00\\x00\\x00\\x00\\x00@\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00')
    >>>

    For format details, see
      ftp://ftp.ncbi.nlm.nih.gov/pubchem/specifications/pubchem_fingerprints.txt
    """
    fp = text.decode("base64")
    # first 4 bytes are the length (struct.unpack(">I"))
    if fp[:4] != '\x00\x00\x03q':
        raise TypeError("This implementation is hard-coded for 881 bit CACTVS fingerprints")
    return 881, fp[4:].translate(_reverse_bits_in_a_byte_transtable)


def import_decoder(path):
    """Find a decoder function given its full name, as in 'chemfp.decoders.from_cactvs'

    This function imports any intermediate modules, which may be a security concern.
    """
    terms = path.split(".")
    if not terms:
        raise TypeError("missing import name")
    if "" in terms:
        raise TypeError("Empty module name in %r" % (path,))

    # It's impossible to tell if the dotted terms corresponds to
    # module or class/instance attribute lookups, so I don't know
    # which fields are imports and which fields are getattrs. To get
    # around that, I'll import everything, and if that fails I'll
    # remove the deepest term and try again.
    tmp_terms = terms[:]
    while tmp_terms:
        try:
            __import__(".".join(tmp_terms), level=0)
        except ImportError:
            del tmp_terms[-1]
        else:
            break
    # I've imported as deep as possible.
    # Now start from the top and work down with getattr calls
    obj = __import__(terms[0], level=0)
    for i, subattr in enumerate(terms[1:]):
       obj = getattr(obj, subattr, None)
       if obj is None:
           failure_path = ".".join(terms[:i+2])
           raise TypeError(("Unable to import a decoder: "
                            "Could not find %(attr)r from %(path)r") %
                           dict(attr=failure_path, path=path))

    return obj


##### Helper code for dealing with common command-line parameters

_decoding_args = []
_decoder_table = {}
def _A(arg, action, decoder, help):
    _decoding_args.append ( ((arg,), dict(action=action, help=help)) )
    _decoder_table[arg.lstrip("-").replace("-","_")] = decoder

_A("--binary", "store_true", from_binary_lsb,
   "Encoded with the characters '0' and '1'. Bit #0 comes first. Example: 00100000 encodes the value 4")
_A("--binary-msb", "store_true", from_binary_msb,
   "Encoded with the characters '0' and '1'. Bit #0 comes last. Example: 00000100 encodes the value 4")
_A("--hex", "store_true", from_hex,
   "Hex encoded. Bit #0 is the first bit (1<<0) of the first byte. Example: 01f2 encodes the value \\x01\\xf2 = 498")
_A("--hex-lsb", "store_true", from_hex_lsb,
   "Hex encoded. Bit #0 is the eigth bit (1<<7) of the first byte. Example: 804f encodes the value \\x01\\xf2 = 498")
_A("--hex-msb", "store_true", from_hex_msb,
   "Hex encoded. Bit #0 is the first bit (1<<0) of the last byte. Example: f201 encodes the value \\x01\\xf2 = 498")
_A("--base64", "store_true", from_base64,
   "Base-64 encoded. Bit #0 is first bit (1<<0) of first byte. Example: AfI= encodes value \\x01\\xf2 = 498")
_A("--cactvs", "store_true", from_cactvs,
   help="CACTVS encoding, based on base64 and includes a version and bit length")
_A("--decoder", "store", None,
    help="import and use the DECODER function to decode the fingerprint")

def _add_decoding_group(parser):
    decoding_group = parser.add_argument_group("Fingerprint decoding options")
    for (args, kwargs) in _decoding_args:
        decoding_group.add_argument(*args, **kwargs)

def _extract_decoder(parser, namespace):
    """An internal helper function for the command-line programs"""
    # Were any command-line decoder arguments specified?
    # Make sure that multiple decoders were not specified
    decoder_name = None
    for arg in _decoder_table:
        if getattr(namespace, arg):
            if decoder_name is not None:
                parser.error("Cannot decode with both --%(old_arg)s and --%(arg)s" % 
                             dict(old_arg=decoder_name, arg=arg))
            decoder_name = arg
    # When in doubt, assume a hex decoder
    if decoder_name is None:
        decoder_name = "hex"

    # If --decoder was specified, do the import and return (name, decoder)
    if decoder_name == "decoder":
        function_name = getattr(namespace, "decoder")
        fp_decoder = import_decoder(function_name)
        return function_name, fp_decoder

    # Otherwise it's in the decoder table
    fp_decoder = _decoder_table[decoder_name]
    return decoder_name, fp_decoder
