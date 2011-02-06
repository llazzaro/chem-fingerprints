import re
import os
import sys
import binascii
import contextlib

from datetime import datetime

def utcnow():
    return datetime.utcnow().isoformat().split(".", 1)[0]

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
    
####
    
_compression_extensions = {
    ".gz": ".gz",
    ".gzip": ".gz",
    ".bz2": ".bz2",
    ".bzip": ".bz2",
    ".bzip2": ".bz2",
    ".xz": ".xz",
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
    compression = _compression_extensions[ext]

    base, ext = os.path.splitext(base)
    if ext == "":
        # compression found but not the actual format type
        return (default[0], compression)

    # The [1:] is to remove the leading "."
    format_name = ext[1:].lower()
    
    return (format_name, compression)

####

def open_output(destination):
    if destination is None:
        return sys.stdout
    if not isinstance(destination, basestring):
        return destination
    base, ext = os.path.splitext(destination)
    ext = ext.lower()
    if ext not in _compression_extensions:
        return open(destination, "w")
    else:
        return open_compressed_output(destination, ext)
    
def open_compressed_output(destination, compression):
    if not compression:
        if destination is None:
            return sys.stdout
        elif isinstance(source, basestring):
            return open(source, "w")
        else:
            return destination

    if compression == ".gz":
        import gzip
        if destination is None:
            return gzip.GzipFile(mode="w", fileobj=sys.stdout)
        elif isinstance(source, basestring):
            return gzip.open(source, "w")
        else:
            return gzip.GzipFile(mode="w", fileobj=destination)

    if compression == ".bz2":
        raise NotImplementedError("bzip2 compression not supported")

    if compression == ".xz":
        raise NotImplementedError("xz compression not supported")

    raise TypeError("unknown compression type %r" % (compression,))

def open_compressed_input_universal(source, compression):
    if not compression:
        if source is None:
            return sys.stdin
        elif isinstance(source, basestring):
            return open(source, "rU")
        else:
            return source

    if compression == ".gz":
        import gzip
        if source is None:
            return gzip.GzipFile(fileobj=sys.stdin)
        elif isinstance(source, basestring):
            return gzip.open(source, "r")
        else:
            return source

    if compression == ".bz2":
        raise NotImplementedError("bzip2 decompression not supported")

    if compression == ".xz":
        raise NotImplementedError("xz decompression not supported")

    raise TypeError("unknown compression typ %r" % (compression,))



def write_fps1_magic(outfile):
    outfile.write("#FPS1\n")

def write_fps1_header(outfile, header):
    lines = []
    if header.num_bits is not None:
        lines.append("#num_bits=%d\n" % header.num_bits)

    if header.software is not None:
        assert "\n" not in header.software
        lines.append("#software=" + header.software.encode("utf8")+"\n")

    if header.type is not None:
        assert "\n" not in header.type
        lines.append("#type=" + header.type.encode("ascii")+"\n") # type cannot contain non-ASCII characters

    if header.source is not None:
        # Ignore newlines in the source filename, if present
        source = header.source.replace("\n", "")
        lines.append("#source=" + source.encode("utf8")+"\n")

    if header.date is not None:
        lines.append("#date=" + header.date.encode("ascii")+"\n") # date cannot contain non-ASCII characters
        
    outfile.writelines(lines)

class _IgnorePipeErrors(object):
    def __enter__(self):
        return None
    
    def __exit__(self, type, value, tb):
        if type is None:
            return False
        import errno
        # Catch any pipe errors, like when piping the output to "head -10"
        if issubclass(type, IOError) and value[0] == errno.EPIPE:
            return True

        return False

ignore_pipe_errors = _IgnorePipeErrors()


def write_fps1_fingerprint(outfile, fp, title):
    outfile.write("%s %s\n" % (binascii.hexlify(fp), title))

def write_fps1_output(reader, destination):
    hexlify = binascii.hexlify
    outfile = open_output(destination)
    with contextlib.closing(open_output(destination)) as outfile:
        with ignore_pipe_errors:
            write_fps1_magic(outfile)
            write_fps1_header(outfile, reader.header)

            for (fp, title) in reader:
                outfile.write("%s %s\n" % (hexlify(fp), title))
