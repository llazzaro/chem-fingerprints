from __future__ import division

from cStringIO import StringIO
from __builtin__ import open as _builtin_open
import binascii
import _chemfp
import re
import heapq

from chemfp import bitops
import ctypes

from . import search

BLOCKSIZE = 20000

class FPHeader(object):
    def __init__(self, num_bits=None, software=None, params=None,
                 source=None, date=None):
        self.num_bits = num_bits
        self.software = software
        self.params = params
        self.source = source
        self.date = date

    def _get_num_bytes_per_fp(self):
        if self.num_bits is None:
            return None
        return (self.num_bits+7)//8
    num_bytes_per_fp = property(_get_num_bytes_per_fp)



def _open_filename_or_stream(source, mode):
    if isinstance(source, basestring):
        if source[-3:].lower() == ".gz":
            import gzip
            filename = source[:-3]
            return source, gzip.open(source)

        return source, _builtin_open(source, mode)

    name = getattr(source, "name", "<unknown>")
    return name, source

#def open_fps_mmap(source):
    

def open_fps(source):
    filename, infile = _open_filename_or_stream(source, "rU")
    header, lineno, block = read_header(infile, filename)
    
    return FPSReader(infile, header, lineno, block)

def open_fpb(source):
    name, infile = _open_filename_or_stream(source, "rB")
    return FPBReader(infile)

def open(filename):
    s = filename.lower()
    if s[-3:] == ".gz":
        s = s[:-3]
    if s[-4:] == ".fps":
        return open_fps(filename)
    if s[-4:] == ".fpb":
        return open_fpb(filename)
    raise TypeError("Unknown extension")

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
    def __init__(self, infile, header, first_fp_lineno, first_fp_block):
        self._infile = infile
        self._filename = getattr(infile, "name", "<unknown>")
        self.header = header
        self._first_fp_lineno = first_fp_lineno
        self._first_fp_block = first_fp_block
        self._expected_hex_len = 2*header.num_bytes_per_fp
        self._hex_len_source = "size in header"

        try:
            self._seekpos = infile.tell()
        except IOError:
            self._seekpos = None
        self._at_start = True
        self._it = None
        
    def reset(self):
        self._it = None
        if self._seekpos is None:
            raise IOError("ASDFASDF")
        self._infile.seek(self._seekpos)
        self._at_start = True
        
    def iter_blocks(self):
        if self._at_start:
            self._at_start = False
            if self._first_fp_block is not None:
                yield self._first_fp_block
            else:
                return
        else:
            raise TypeError("Already iterating")
        for block in _read_blocks(self._infile):
            yield block

    def iter_rows(self):
        unhexlify = binascii.unhexlify
        lineno = self._first_fp_lineno
        expected_hex_len = self._expected_hex_len
        for block in self.iter_blocks():
            for line in block.splitlines():
                errcode = _chemfp.fps_line_validate(expected_hex_len, line)
                if errcode:
                    raise Error
                fields = line.split()
                fp = unhexlify(fields[0])
                yield (fp,) + tuple(fields[1:])
                lineno += 1

    def iter_fingerprints(self):
        unhexlify = binascii.unhexlify
        lineno = self._first_fp_lineno
        expected_hex_len = self._expected_hex_len
        for block in self.iter_blocks():
            for line in block.splitlines(True):
                errcode = _chemfp.fps_line_validate(expected_hex_len, line)
                if errcode:
                    raise Exception(errcode, line)
                fields = line.split(None, 2)
                yield unhexlify(fields[0]), fields[1]
                lineno += 1
        
    def __iter__(self):
        return self.iter_fingerprints()

    def _chemfp_tanimoto_knearest_search_batch(self, queries, n, threshold):
        return search.block_tanimoto_knearest_search_batch(queries, self, n, threshold)

    def _chemfp_tanimoto_count_batch(self, queries, threshold):
        return search.block_tanimoto_count_batch(queries, self, threshold)


# XXX Use Python's warning system
def warn_to_stderr(filename, lineno, message):
    where = _where(filename, lineno)
    sys.stderr.write("WARNING: %s at %s\n" % (message, where))

_whitespace = re.compile(r"[ \t\n]")
def read_header(f, filename, warn=warn_to_stderr):
    header = FPHeader()

    lineno = 0
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
                if header.num_bits is None:
                    # We can figure this out from the fingerprint on the first line
                    m = _whitespace.search(block)
                    if m is None:
                        raise TypeError(block)
                    i = m.end()-1 # Back up from the whitespace
                    if i % 2 == 1:
                        raise TypeError(block)
                    header.num_bits = i * 4
                    
                return header, lineno, block

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
            lineno += 1

            # Right! We've got a line. Check if it's magic
            # This is the only line which cannot contain a '='
            if lineno == 1:
                if line.rstrip() == "FPS1":
                    continue
                assert "=" not in line, line
                
            assert "=" in line, line
            key, value = line.split("=", 1)
            key = key.strip()
            value = value.strip()

            if key == "num_bits":
                try:
                    header.num_bits = int(value)
                    if not (header.num_bits > 0):
                        raise ValueError
                except ValueError:
                    raise TypeError(
                        "num_bits header must be a positive integer, not %r: %s" %
                        (value, _where(filename, lineno)))
            elif key == "software":
                header.software = value
            elif key == "type":
                #header.params = normalize_params(value)
                header.params = value
            elif key == "source":
                header.source = value
            elif key == "date":
                header.date = value
            else:
                print "UNKNOWN", repr(line), repr(key), repr(value)
                #warn(filename, lineno, "Unknown header %r" % (value,))

    # Reached the end of file. No fingerprint lines and nothing left to process.
    return header, lineno, None

def fps_to_in_memory(fps_reader, header=None):
    if header is None:
        header = fps_reader.header
    num_bits = header.num_bits
    assert num_bits
    id_bins = [[] for i in xrange(num_bits+1)]
    popcount_bins = [StringIO() for i in xrange(num_bits+1)]
    for i, (fp, id) in enumerate(fps_reader):
        popcount = bitops.byte_popcount(fp)
        id_bins[popcount].append(id)
        popcount_bins[popcount].write(fp)

    # Do this so we don't always double the memory use
    popcount_arenas = []
    for i in xrange(num_bits+1):
        popcount_arenas.append(popcount_bins[i].getvalue())
        popcount_bins[i] = None
    
    return InMemory(header, popcount_arenas, id_bins)

class InMemory(object):
    def __init__(self, header, arenas, arena_ids):
        self.header = header
        self.arenas = arenas
        self.arena_ids = arena_ids

    def reset(self):
        pass

    def iter_fingerprints(self):
        target_fp_size = self.header.num_bytes_per_fp
        for arena, arena_id in zip(self.arenas, self.arena_ids):
            for index, offset in enumerate(xrange(0, len(arena), target_fp_size)):
                yield arena[offset:offset+target_fp_size], arena_id[index]

    def __iter__(self):
        return self.iter_fingerprints()
    
    def _chemfp_tanimoto_knearest_search_batch(self, queries, k, threshold):
        return search.arena_tanimoto_knearest_search_batch(queries, self.arenas, self.arena_ids, k, threshold)
    
    def _chemfp_tanimoto_count_batch(self, queries, threshold):
        return search.arena_tanimoto_count_batch(queries, self.arenas, threshold)
