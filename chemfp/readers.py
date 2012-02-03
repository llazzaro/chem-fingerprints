from __future__ import absolute_import

from cStringIO import StringIO
from __builtin__ import open as _builtin_open
import binascii
import _chemfp
import re
import sys
import heapq
import itertools
import ctypes

from . import load_fingerprints, Metadata, ParseError
from . import fps_search
from . import io

# I tried a wide range of sizes for my laptop, with both compressed
# and uncompressed files, and found that the best size was around
# 2**17. Actually, 2**16.8 was the absolute best, which gives
BLOCKSIZE=11400
# (BTW, the compressed time took 1.3x the uncompressed time)

class FPSParseError(ParseError):
    def __init__(self, errcode, lineno, filename):
        self.errcode = errcode
        self.lineno = lineno
        self.filename = filename
    def __repr__(self):
        return "FPSParseError(%d, %d, %s)" % (self.errcode, self.lineno, self.filename)
    def __str__(self):
        msg = _chemfp.strerror(self.errcode)
        msg += " at line %d" % (self.lineno,)
        if self.filename is not None:
            msg += " of %r" % (self.filename,)
        return msg


def open_fps(source, format=None):
    format_name, compression = io.normalize_format(source, format)
    if format_name != "fps":
        raise ValueError("Unknown format %r" % (format_name,))

    infile = io.open_compressed_input_universal(source, compression)
    filename = io.get_filename(source)

    metadata, lineno, block = read_header(infile, filename)
    return FPSReader(infile, metadata, lineno, block)


# This never buffers
def _read_blocks(infile):
    while 1:
        block = infile.read(BLOCKSIZE)
        if not block:
            break
        if block[-1:] == "\n":
            yield block
            continue
        line = infile.readline()
        if not line:
            # Note: this might not end with a newline!
            yield block
            break
        yield block + line

            

class FPSReader(object):
    _search = fps_search
    def __init__(self, infile, metadata, first_fp_lineno, first_fp_block):
        self._infile = infile
        self._filename = getattr(infile, "name", "<unknown>")
        self.metadata = metadata
        self._first_fp_lineno = first_fp_lineno
        self._first_fp_block = first_fp_block
        self._expected_hex_len = 2*metadata.num_bytes
        self._hex_len_source = "size in header"

        self._at_start = True
        self._it = None
        self._block_reader = None

# Not sure if this is complete. Also, should have a context manager
#    def close(self):
#        self._infile.close()
        
    def iter_blocks(self):
        if self._block_reader is None:
            self._block_reader = iter(self._iter_blocks())
        return self._block_reader

    def _iter_blocks(self):
        if not self._at_start:
            raise TypeError("Already iterating")
        
        self._at_start = False

        if self._first_fp_block is None:
            return
        
        block_stream = _read_blocks(self._infile)
        yield self._first_fp_block
        for block in block_stream:
            yield block

    def iter_arenas(self, arena_size = 1000):
        id_fps = iter(self)
        if arena_size is None:
            yield load_fingerprints(id_fps,
                                    metadata = self.metadata,
                                    reorder = False)
            return
        while 1:
            arena = load_fingerprints(itertools.islice(id_fps, 0, arena_size),
                                      metadata = self.metadata,
                                      reorder = False)
            if not arena:
                break
            yield arena
        
    def iter_rows(self):
        unhexlify = binascii.unhexlify
        lineno = self._first_fp_lineno
        expected_hex_len = self._expected_hex_len
        for block in self.iter_blocks():
            for line in block.splitlines(True):
                err = _chemfp.fps_line_validate(expected_hex_len, line)
                if err:
                    raise FPSParseError(err, lineno, self._filename)
                yield line[:-1].split("\t")
                lineno += 1

    def __iter__(self):
        unhexlify = binascii.unhexlify
        lineno = self._first_fp_lineno
        expected_hex_len = self._expected_hex_len
        for block in self.iter_blocks():
            for line in block.splitlines(True):
                err, id_fp = _chemfp.fps_parse_id_fp(expected_hex_len, line)
                if err:
                    # Include the line?
                    raise FPSParseError(err, lineno, self._filename)
                yield id_fp
                lineno += 1

    def _check_at_start(self):
        if not self._at_start:
            raise TypeError("FPS file is not at the start of the file; cannot search")


    def count_tanimoto_hits_fp(self, query_fp, threshold=0.7):
        self._check_at_start()
        return fps_search.count_tanimoto_hits_fp(query_fp, self, threshold)

    def count_tanimoto_hits_arena(self, queries, threshold=0.7, arena_size=100):
        self._check_at_start()
        return fps_search.count_tanimoto_hits_arena(queries, self, threshold)

    def threshold_tanimoto_search_fp(self, query_fp, threshold=0.7):
        self._check_at_start()
        return fps_search.threshold_tanimoto_search_fp(query_fp, self, threshold)

    def threshold_tanimoto_search_arena(self, queries, threshold=0.7, arena_size=100):
        self._check_at_start()
        return fps_search.threshold_tanimoto_search_arena(queries, self, threshold)

    def knearest_tanimoto_search_fp(self, query_fp, k=3, threshold=0.7):
        self._check_at_start()
        return fps_search.knearest_tanimoto_search_fp(query_fp, self, k, threshold)

    def knearest_tanimoto_search_arena(self, queries, k=3, threshold=0.7, arena_size=100):
        self._check_at_start()
        return fps_search.knearest_tanimoto_search_arena(queries, self, k, threshold)

def _where(filename, lineno):
    if filename is None:
        return "line %d" % (lineno,)
    else:
        return "%r line %d" % (filename, lineno)

# XXX Use Python's warning system
def warn_to_stderr(filename, lineno, message):
    where = _where(filename, lineno)
    sys.stderr.write("WARNING: %s at %s\n" % (message, where))

def read_header(f, filename, warn=warn_to_stderr):
    metadata = Metadata()

    lineno = 1
    for block in _read_blocks(f):
        # A block must be non-empty
        start = 0
        while 1:
            c = block[start:start+1]
            if c == "":
                # End of the block; get the next one
                break
            if c != '#':
                # End of the header. This block contains the first fingerprint line
                block = block[start:]
                if metadata.num_bits is None:
                    # We can figure this out from the fingerprint on the first line
                    err = _chemfp.fps_line_validate(-1, block)
                    if err:
                        raise FPSParseError(err, lineno, filename)
                    i = block.index("\t")
                    # If you don't specify the number of bits then I'll do it for you.
                    metadata.num_bits = i * 4
                    metadata.num_bytes = i // 2
                    
                return metadata, lineno, block

            start += 1 # Skip the '#'
            end = block.find("\n", start)
            if end == -1:
                # Only happens when the last line of the file contains
                # no newlines. In that case, we're at the last block.
                line = block[start:]
                start = len(block)
            else:
                line = block[start:end]
                start = end+1

            # Right! We've got a line. Check if it's magic
            # This is the only line which cannot contain a '='
            if lineno == 1:
                if line.rstrip() == "FPS1":
                    lineno += 1
                    continue

            if line.startswith("x-") or line.startswith("X-"):
                # Completely ignore the contents of 'experimental' lines
                continue

            if "=" not in line:
                raise TypeError("header line must contain an '=': %r at %s" %
                                (line, _where(filename, lineno)))
            key, value = line.split("=", 1)
            key = key.strip()
            value = value.strip()

            if key == "num_bits":
                try:
                    metadata.num_bits = int(value)
                    metadata.num_bytes = (metadata.num_bits + 7)//8
                    if not (metadata.num_bits > 0):
                        raise ValueError
                except ValueError:
                    raise TypeError(
                        "num_bits header must be a positive integer, not %r: %s" %
                        (value, _where(filename, lineno)))
                metadata.num_bytes = (metadata.num_bits+7)//8
            elif key == "software":
                metadata.software = value.decode("utf8")
            elif key == "type":
                # Should I have an auto-normalization step here which
                # removes excess whitespace?
                #metadata.type = normalize_type(value)
                metadata.type = value
            elif key == "source":
                metadata.sources.append(value)
            elif key == "date":
                metadata.date = value
            elif key == "aromaticity":
                metadata.aromaticity = value
            elif key.startswith("x-"):
                pass
            else:
                #print "UNKNOWN", repr(line), repr(key), repr(value)
                #warn(filename, lineno, "Unknown header %r" % (value,))
                pass
            lineno += 1

    # Reached the end of file. No fingerprint lines and nothing left to process.
    if metadata.num_bits is None:
        metadata.num_bits = 0
        metadata.num_bytes = 0
    return metadata, lineno, None

