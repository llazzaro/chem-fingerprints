import string
import binascii

_lsb_bit_table = {}
_msb_bit_table = {}

_reverse_bits_in_a_byte_transtable = None

def _init():
    # This creates the 256 bit patterns:
    #    00000000, 10000000, 01000000, 11000000 ... 11111111
    def _C(it):
        for i in it:
            yield "0" + i
            yield "1" + i
    to_trans = [None]*256
    for value, bit_pattern in enumerate(_C(_C(_C(_C(_C(_C(_C("01")))))))):
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
    try:
        N = len(text)
        bytes = []
        for i in range(0, N, 8):
            bytes.append(_lsb_bit_table[text[i:i+8]])
        return (N, "".join(bytes))
    except KeyError:
        raise TypeError("Not a binary string")

def from_binary_msb(text):
    try:
        N = len(text)
        bytes = []
        end = N
        start = N-8
        while start > 8:
            bytes.append(_msb_bit_table[text[start:end]])
            end = start
            start -= 8
        bytes.append(_msb_bit_table[text[start:end]])
        return (N, "".join(bytes))
    except KeyError:
        raise TypeError("Not a binary string")


def from_base64(text):
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
    return (None, text.decode("hex"))

def from_hex_msb(text):
    return (None, text.decode("hex")[::-1])

def from_hex_lsb(text):
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
    fp = text.decode("base64")
    # first 4 bytes are the length (struct.unpack(">I"))
    assert fp[:4] == '\x00\x00\x03q', "not 881 bits?"
    return 881, fp[4:].translate(_reverse_bits_in_a_byte_transtable)


def import_decoder(path):
    terms = path.split(".")
    if not terms:
        raise TypeError("missing import name")
    tmp_terms = terms[:]
    while tmp_terms:
        try:
            __import__(".".join(tmp_terms), level=0)
        except ImportError:
            del tmp_terms[-1]
        else:
            break
    obj = __import__(terms[0], level=0)
    for subattr in terms[1:]:
       obj = getattr(obj, subattr)
       if obj is None:
           raise TypeError(("Unable to import a decoder: "
                           "stopped at {attr!r} in {path!r}").format(
                   attr=subattr, path=path))
    return obj


##### Helper code for dealing with common command-line parameters

_decoding_args = []
_decoder_table = {}
def _A(arg, action, decoder, help):
    _decoding_args.append ( ((arg,), dict(action=action, help=help)) )
    _decoder_table[arg.lstrip("-").replace("-","_")] = decoder

_A("--binary", "store_true", from_binary_lsb,
   "encoded using 0s and 1s. #0 is first")
_A("--binary-msb", "store_true", from_binary_msb,
   "encoded using 0s and 1s. #0 is last")
_A("--hex", "store_true", from_hex,
   "hex-encoded bytes. #0 is first bit of first byte")
_A("--hex-lsb", "store_true", from_hex_lsb,
   "hex-encoded bytes. #0 is eigth bit of first byte")
_A("--hex-msb", "store_true", from_hex_msb,
   "hex-encoded bytes. #0 is first bit of last byte")
_A("--base64", "store_true", from_base64,
   "base64 encoded bytes. #0 is first bit of first byte")
_A("--cactvs", "store_true", from_cactvs,
   help="CACTVS encoding, based on base64")
_A("--decoder", "store", None,
    help="import DECODER and use it to decode the fingerprint")

def _add_decoding_group(parser):
    decoding_group = parser.add_argument_group("Fingerprint decoding options")
    for (args, kwargs) in _decoding_args:
        decoding_group.add_argument(*args, **kwargs)

def _extract_decoder(parser, namespace):
    decoder_name = None
    for arg in _decoder_table:
        if getattr(namespace, arg):
            if decoder_name is not None:
                parser.error("Cannot decode with both --{old_arg} and --{arg}".format(
                        old_arg=decoder_name, arg=arg))
            decoder_name = arg
    if decoder_name is None:
        decoder_name = "hex"

    if decoder_name == "decoder":
        function_name = getattr(namespace, "decoder")
        fp_decoder = import_decoder(function_name)
        return function_name, fp_decoder

    fp_decoder = _decoder_table[decoder_name]
    return decoder_name, fp_decoder
