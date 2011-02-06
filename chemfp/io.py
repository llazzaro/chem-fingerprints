import re
import os

class Header(object):
    def __init__(self, num_bits=None, software=None, type=None,
                 source=None, date=None):
        assert num_bits is None or isinstance(num_bits, int)
        self.num_bits = num_bits
        self.software = software
        self.type = type
        self.source = source
        self.date = date

    @property
    def num_bytes_per_fp(self):
        if self.num_bits is None:
            return None
        return (self.num_bits+7)//8

class FPIterator(object):
    def __init__(self, header, iterator):
        self.header = header
        self.iterator = iterator

    def iter_fingerprints(self):
        return self.iterator

    def __iter__(self):
        return self.iterator
    

_compression_extensions = {
    ".gz": ".gz",
    ".gzip": ".gz",
    ".bz2": ".bz2",
    ".bzip": ".bz2",
    ".bzip2": ".bz2",
    }
_compression_regex = "|".join(re.escape(ext) for ext in _compression_extensions)

_format_pat = re.compile("^([a-zA-Z0-9]+)(" + _compression_regex + ")?$")

def normalize_format(source, format, default=("fps", None)):
    if source is None:
        # Read from stdin
        filename = None
    elif isinstance(source, basestring):
        # This is a filename
        filename = source
    elif hasattr(source, "read"):
        # This is a Python file object
        filename = getattr(source, "name", None)
    else:
        raise TypeError("Unknown source type %r" % (source,))
    
    if format is not None:
        # This must be of the form <name> [. <compression> ]
        m = _format_pat.match(format)
        if m is None:
            raise TypeError("Unknown format %r" % (format,))

        if m.group(2):
            # This is a compression
            compression = _compression_extensions[m.group(2)]
        else:
            compression = ""

        format_name = m.group(1)
        return (format_name, compression)

    if filename is None:
        # Reading from stdin or an unnamed file-like object with no
        # specified format Not going to sniff the input. Instead, just
        return default

    

    # The filename could have 0, 1 or 2 extensions
    base, ext = os.path.splitext(filename)
    if ext == "":
        # No extensions, use the default
        return default
    ext = ext.lower()

    # If it's not a compression extension then it's a format indicator
    if ext not in _compression_extensions:
        # the [1:] is to remove the leading "."
        return (ext[1:], "")

    # Found a compression, now look for the format
    compression = ext

    base, ext = os.path.splitext(base)
    if ext == "":
        # compression found but not the actual format type
        return (default[0], compression)

    # The [1:] is to remove the leading "."
    format_name = ext[1:].lower()
    
    return (format_name, compression)
