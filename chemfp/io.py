from __future__ import with_statement

import re
import os
import sys
import binascii

from datetime import datetime

if sys.platform.startswith("win"):
    DEV_STDIN = "CON"
else:
    if os.path.exists("/dev/stdin"):
        DEV_STDIN = "/dev/stdin"
    else:
        DEV_STDIN = None

def utcnow():
    return datetime.utcnow().isoformat().split(".", 1)[0]

# XXX should this be here?
def remove_special_characters_from_id(id):
    if "\n" in id:
        id = id.splitlines()[0]
    if "\t" in id:
        id = id.replace("\t", "")
    if " " in id:
        id = id.strip()
    return id

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
        raise ValueError("Unknown source type %r" % (source,))
    
    if format is not None:
        # This must be of the form <name> [. <compression> ]
        m = _format_pat.match(format)
        if m is None:
            if "." in format:
                if _format_pat.match(format.split(".")[0]):
                    raise ValueError("Unsupported compression in format %r" % (format,))
            raise ValueError("Incorrect format syntax %r" % (format,))

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
        return getattr(source, "name", None)

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
        elif isinstance(destination, basestring):
            return open(destination, "w")
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
            if not os.path.exists("/dev/stdout"):
                raise NotImplementedError("Cannot write bz2 compressed data to stdout on this platform")
            return bz2.BZ2File("/dev/stdout", "w")
        elif isinstance(destination, basestring):
            return bz2.BZ2File(destination, "w")
        else:
            raise NotImplementedError("bzip2 compression to file-like objects is not supported")

    if compression == ".xz":
        raise NotImplementedError("xz compression is not supported")

    raise ValueError("Unknown compression type %r" % (compression,))

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
                raise NotImplementedError("Cannot compressed bzip2 data from stdin on this platform")
            return bz2.BZ2File("/dev/stdin", "rU")
        elif isinstance(source, basestring):
            return bz2.BZ2File(source, "rU")
        else:
            # Well, I could emulate it, but I'm not going to
            raise NotImplementedError("bzip decompression from file-like objects is not supported")

    if compression == ".xz":
        raise NotImplementedError("xz decompression is not supported")

    raise ValueError("Unknown compression type %r" % (compression,))



def write_fps1_magic(outfile):
    outfile.write("#FPS1\n")

def write_fps1_header(outfile, metadata):
    lines = []
    if metadata.num_bits is not None:
        lines.append("#num_bits=%d\n" % metadata.num_bits)

    if metadata.type is not None:
        assert "\n" not in metadata.type
        lines.append("#type=" + metadata.type.encode("ascii")+"\n") # type cannot contain non-ASCII characters

    if metadata.software is not None:
        assert "\n" not in metadata.software
        lines.append("#software=" + metadata.software.encode("utf8")+"\n")

    if metadata.aromaticity is not None:
        lines.append("#aromaticity=" + metadata.aromaticity.encode("ascii") + "\n")

    for source in metadata.sources:
        # Ignore newlines in the source filename, if present
        source = source.replace("\n", "")
        lines.append("#source=" + source.encode("utf8")+"\n")

    if metadata.date is not None:
        lines.append("#date=" + metadata.date.encode("ascii")+"\n") # date cannot contain non-ASCII characters

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


def write_fps1_fingerprint(outfile, fp, id):
    if "\t" in id:
        raise ValueError("Fingerprint ids must not contain a tab: %r" % (id,))
    if "\n" in id:
        raise ValueError("Fingerprint ids must not contain a newline: %r" % (id,))
    if not id:
        raise ValueError("Fingerprint ids must not be the empty string")
    
    outfile.write("%s\t%s\n" % (binascii.hexlify(fp), id))


# This is a bit of a hack. If I open a file then I want to close it,
# but if I use stdout then I don't want to close it.

class _closing_output(object):
    def __init__(self, destination):
        self.destination = destination
        self.output = open_output(destination)
    def __enter__(self):
        return self.output
    def __exit__(self, *exec_info):
        if isinstance(self.destination, basestring):
            self.output.close()

def write_fps1_output(reader, destination, metadata=None):
    if metadata is None:
        metadata = reader.metadata
    hexlify = binascii.hexlify
    with _closing_output(destination) as outfile:
        with ignore_pipe_errors:
            write_fps1_magic(outfile)
            write_fps1_header(outfile, metadata)

            for i, (id, fp) in enumerate(reader):
                if "\t" in id:
                    raise ValueError("Fingerprint ids must not contain a tab: %r in record %d" %
                                     (id, i+1))
                if "\n" in id:
                    raise ValueError("Fingerprint ids must not contain a newline: %r in record %d" %
                                     (id, i+1))
                if not id:
                    raise ValueError("Fingerprint ids must not be the empty string in record %d" %
                                     (i+1,))
                outfile.write("%s\t%s\n" % (hexlify(fp), id))
