from __future__ import with_statement

import re
import os
import sys
import binascii

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

def normalize_format(source, format, default=("fps", "")):
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

def get_filename(source):
    if source is None:
        return None
    elif isinstance(source, basestring):
        return source
    else:
        return getattr(infile, "name", None)

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
        elif isinstance(destination, basestring):
            return gzip.open(destination, "w")
        else:
            return gzip.GzipFile(mode="w", fileobj=destination)

    if compression == ".bz2":
        import bz2
        if destination is None:
            if not exists("/dev/stdout"):
                raise NotImplementedError("Cannot write bz2 compressed data to stdout on this platform")
            return bz2.BZ2File("/dev/stdout", "w")
        elif isinstance(destination, basestring):
            return bz2.BZ2File(destination, "w")
        else:
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
            # GzipFile doesn't have a "U"/universal file mode?
            return gzip.GzipFile(fileobj=sys.stdin)
        elif isinstance(source, basestring):
            return gzip.open(source, "r")
        else:
            return gzip.GzipFile(fileobj=source)

    if compression == ".bz2":
        import bz2
        if source is None:
            # bz2 doesn't support Python objects. On some platforms
            # I can work around the problem
            if not os.path.exists("/dev/stdin"):
                raise NotImplementedError("Cannot read bz2 compressed data from stdin on this platform")
            return bz2.BZ2File("/dev/stdin", "rU")
        elif isinstance(source, basestring):
            return bz2.BZ2File(source, "rU")
        else:
            # Well, I could emulate it, but I'm not going to
            raise NotImplementedError("Cannot read bz2 compress data from a Python stream")

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


# This is a bit of a hack. If I open a file then I want to close it,
# but if I use stdout then I don't want to close it.

class _closing_output(object):
    def __init__(self, destination):
        self.output = open_output(destination)
    def __enter__(self):
        return self.output
    def __exit__(self, *exec_info):
        if self.output is not sys.stdout:
            self.output.close()

def write_fps1_output(reader, destination):
    hexlify = binascii.hexlify
    with _closing_output(destination) as outfile:
        with ignore_pipe_errors:
            write_fps1_magic(outfile)
            write_fps1_header(outfile, reader.header)

        
            for i, (fp, title) in enumerate(reader):
                outfile.write("%s %s\n" % (hexlify(fp), title))
                if i == 100:
                    return
