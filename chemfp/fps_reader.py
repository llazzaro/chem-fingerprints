import sys
import re
import itertools
import ctypes

from chemfp import bitops
import _chemfp

_fps_header = re.compile(r"#FPS(\d+)$")

def normalize_params(params):
    return " ".join(params.split())

class TanimotoHeap(ctypes.Structure):
    _fields_ = [("size", ctypes.c_int),
                ("k", ctypes.c_int),
                ("unique_idx", ctypes.c_int),
                ("_reserved", ctypes.c_int),
                ("threshold", ctypes.c_double),
                ("indicies", ctypes.c_void_p),
                ("scores", ctypes.c_void_p),
                ("id_starts", ctypes.c_void_p),
                ("id_lens", ctypes.c_void_p),
                ]

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

_open = open

def split_header(lineno, line):
    words = line[1:].split("=", 1)
    if not words:
        # Unknown format - skip
        return None, None
    return words[0].strip(), words[1].strip()

def _where(filename, lineno):
    if filename is not None:
        return "file %r, line %d" % (filename, lineno)
    return "line %d" % (lineno,)

# XXX Use Python's warning system
def warn_to_stderr(filename, lineno, message):
    where = _where(filename, lineno)
    sys.stderr.write("WARNING: %s at %s\n" % (message, where))

def _read_header(f, filename, warn=warn_to_stderr):
    header = FPHeader()

    lineno = 0
    line = None
    first_fp_line = None
    try:
        seekpos = f.tell()
    except IOError:
        seekpos = None

    while 1:
        if seekpos is not None:
            seekpos = f.tell()
        line = f.readline()
        if not line:
            # End of file. No fingerprint data
            break
        lineno += 1

        if lineno == 1:
            # See if this is the magic line
            m = _fps_header.match(line)
            if m is not None:
                version = m.group(1)
                if version != "1":
                    warn(filename, lineno, "Unexpected FPS version number %s" % (version,))
                continue

            # Not the magic line. Either a header or fingerprint data
            warn(filename, lineno, "Expected '#FPS1'")
            # Fall through

        if line[:1] != "#":
            # At a fingerprint line
            # XXX Make sure the line is valid
            if seekpos is not None:
                f.seek(seekpos)
                return header, None, lineno, seekpos
            else:
                return header, line, lineno, None
            
        key, value = split_header(lineno, line)
        if key is None:
            # How future compatible should this be?
            raise TypeError("Unexpected header %r at %s" % (line, _where(filename, lineno)))

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
        elif key == "params":
            header.params = normalize_params(value)
        elif key == "source":
            header.source = value
        elif key == "date":
            header.date = value
        else:
            warn(filename, lineno, "Unknown header %r" % (value,))

    # Reached the end of file. No fingerprint lines and nothing left to process.
    return header, None, lineno, seekpos

def open(filename_or_stream, warn=warn_to_stderr):
    if isinstance(filename_or_stream, basestring):
        # Yes, read binary. XXX how important is that?
        f = _open(filename_or_stream, "rb")
    else:
        f = filename_or_stream
    
    header, first_fp_line, lineno, seekpos  = _read_header(f, warn)
    fps_io = FileFPS(f, first_fp_line, seekpos)
    return FPSReader(fps_io, header, lineno)

def open_in_memory(filename_or_stream, warn=warn_to_stderr):
    if isinstance(filename_or_stream, basestring):
        f = _open(filename_or_stream, "rb")
    else:
        f = filename_or_stream
    header, first_fp_line, lineno, seekpos = _read_header(f, warn)
    data = ""
    if first_fp_line is not None:
        data += first_fp_line
    data += f.read()
    assert data[-1:] == "\n"
    fps_io = InMemoryFPS(data)
    return FPSReader(fps_io, header, lineno)



class Location(object):
    def __init__(self, filename=None, lineno=0):
        self.filename = filename
        self.lineno = lineno
    def where(self):
        return _where(self.filename, self.lineno)

class HexFPRecord(object):
    __slots__ = ["hex_fp", "id"]
    def __init__(self, hex_fp, id):
        self.hex_fp = hex_fp
        self.id = id
    def _byte_fp(self):
        return self.hex_fp.decode("hex")
    byte_fp = property(_byte_fp)

class ByteFPRecord(object):
    __slots__ = ["byte_fp", "id"]
    def __init__(self, byte_fp, id):
        self.byte_fp = byte_fp
        self.id = id
    def _hex_fp(self):
        return self.byte_fp.encode("hex")
    hex_fp = property(_hex_fp)


def _bad_fp_line(msg, filename, lineno):
    raise TypeError("%s: %s" % (msg, _where(filename, lineno)))

BLOCKSIZE = 2000000

def find_block_boundaries(s):
    boundaries = []
    start = 0
    n = len(s)
    while start < n:
        end = start + BLOCKSIZE
        if end > n:
            boundaries.append( (start, n) )
            break
        end = s.index("\n", end)+1
        boundaries.append( (start, end) )
        start = end
    return boundaries

class InMemoryFPS(object):
    def __init__(self, fps_data):
        self.fps_data = fps_data
        assert fps_data[-1:] == "\n"
        self._block_boundaries = find_block_boundaries(fps_data)
    def reset(self):
        pass
    def read_blocks(self):
        yield self.fps_data
    def read_lines(self):
        #yield self.fps_data.splitlines(True)
        #return
        for (start, end) in self._block_boundaries:
            yield self.fps_data[start:end].splitlines(True)
            
class FileFPS(object):
    def __init__(self, f, first_fp_line=None, seekpos=None):
        self.f = f
        self.first_fp_line = first_fp_line
        self.seekpos = seekpos
    def reset(self):
        if self.seekpos is None:
            raise AssertionError("can not seek")
        self.f.seek(self.seekpos)
    def read_blocks(self):
        if self.first_fp_line:
            yield self.first_fp_line
        while 1:
            block = self.f.read(BLOCKSIZE) + self.f.readline()
            if not block:
                break
            yield block
    def read_lines(self):
        if self.first_fp_line:
            yield [self.first_fp_line]
        while 1:
            lines = self.f.readlines(BLOCKSIZE)
            if not lines:
                break
            yield lines

class FPSReader(object):
    def __init__(self, fps_io, header, lineno=None):
        self._fps_io = fps_io
        self.header = header
        if lineno is None:
            lineno = 1
        self._lineno = lineno

    def reset(self):
        self._fps_io.reset()

    def fp_iter(self, loc=None):
        lineno = self._lineno
        if loc is not None:
            loc.filename = self._f.name
            loc.lineno = lineno
        hex_len = self.header.num_bytes_per_fp * 2
        for lines in self._fps_io.read_lines():
            for line in lines:
                # Because this is reference code, I'm doing an extra
                # layer of checking to make sure the content is valid.
                err = _chemfp.fps_line_validate(hex_len, line)
                if err == 0:
                    words = line.split()
                    if loc is not None:
                        loc.lineno = lineno
                    yield words[0].decode("hex"), words[1]
                    lineno += 1
                else:
                    raise AssertionError(_chemfp.strerror(err))

    def __iter__(self):
        return self.fp_iter()

    # calc_tanimoto?
    def _chemfp_tanimoto(self, query, threshold=0.0):
        hex_query = query.encode("hex")
        min_line_width = len(hex_query) + 3 # 3 because of space + id + newline

        num_found_ptr = ctypes.c_int()
        lineno_ptr = ctypes.c_int()
        lineno_ptr.value = self._lineno
        id_starts = None
        id_lens = None
        scores = None

        num_allocated = -1
        for block in self._fps_io.read_blocks():
            if len(block) // min_line_width > num_allocated:
                num_allocated = len(block) // min_line_width + 20
                id_starts = (ctypes.c_void_p*num_allocated)()
                id_lens = (ctypes.c_int*num_allocated)()
                scores = (ctypes.c_double*num_allocated)()

            n = _chemfp.fps_tanimoto(hex_query, block, threshold,
                                     num_found_ptr, id_starts, id_lens, scores, lineno_ptr)
            if n < 0:
                raise AssertionError("ERER %d %s" % (n, _chemfp.strerror(n)))
                raise _error(n, self._f, lineno_ptr.value)
            for i in range(num_found_ptr.value):
                yield ctypes.string_at(id_starts[i], id_lens[i]), scores[i]

    # count_tanimoto?
    def _chemfp_tanimoto_count(self, query, threshold=0.0):
        hex_query = query.encode("hex")
        lineno_ptr = (ctypes.c_int)()
        lineno_ptr.value = self._lineno
        count_ptr = (ctypes.c_int)()
        count_ptr.value = 0

        for block in self._fps_io.read_blocks():
            err = _chemfp.fps_tanimoto_count(hex_query, block, threshold,
                                             count_ptr, lineno_ptr)
            if err < 0:
                raise _error(err, self._f, lineno_ptr.value)
        return count_ptr.value

    # knearest_tanimoto ?
    def _chemfp_tanimoto_knearest(self, query, k, threshold=0.0):
        assert k >= 0
        if k == 0:
            return []

        hex_query = query.encode("hex")
        indicies = (ctypes.c_int*k)()
        scores = (ctypes.c_double*k)()
        id_starts = (ctypes.c_void_p*k)()
        id_lens = (ctypes.c_int*k)()

        lineno_ptr = (ctypes.c_int)()
        lineno_ptr.value = self._lineno

        identifiers = {}
        # Use ctypes to alloc raw space
        # XXX add class methods?
        heap = TanimotoHeap()
        _chemfp.fps_heap_init(heap, k, threshold, indicies, scores,
                              id_starts, id_lens)
        num_allocated = -1
        for target_block in self._fps_io.read_blocks():
            # Returns the current size of the heap
            err = _chemfp.fps_heap_update_tanimoto(heap, hex_query, target_block)
            if err < 0:
                raise _error(err, self._f, lineno_ptr.value)
            # This is an O(k * N) operation
            # If k depends on N, this becomes O(N**2)
            # The algorithm can be made more efficient for those cases.
            # (Eg, have a max id length, or allow memory allocation)
            for i in range(heap.size):
                if id_starts[i] > 0:
                    identifiers[indicies[i]] = ctypes.string_at(id_starts[i], id_lens[i])
                    id_starts[i] = 0

        _chemfp.fps_heap_finish_tanimoto(heap)
        return [(identifiers[indicies[i]], scores[i]) for i in range(heap.size)]


def test():
    reader = open("nci.fps")
    print reader.header
    loc = Location()
    for i, rec in enumerate(reader.fp_iter(loc)):
        print i, len(rec[0]), rec[1], loc.lineno
        if i == 10:
            break
    reader.reset()
    for i, rec in enumerate(reader.fp_iter(loc)):
        print i, len(rec[0]), rec[1], loc.lineno
        if i == 10:
            break


FP = "200206000000000402800e00040140010100014008206200000200c0082200080200500201c9804100270538000000402000a2040080c1240001c2c2004600200c200c04020800200410a0001490000200a803c018005400c80c00000000810100840000880064a0124010000000080102060142400110200a00000004800000".decode("hex")

FP = "600601000200008002000620000600082c0000000000000004000080010400023c00402aa10000070002063820000080016200000000e000010000020020a01008900a0000003880012008020280000000800060002020084008080000200000400000000000020000000000880008000007b1050001100100188001d0000000".decode("hex")

import threading
import Queue
import itertools

request_queue = Queue.Queue()

def add_task(response, f, *args, **kwargs):
    request_queue.put( (response, f, args, kwargs) )

class WorkerThread(threading.Thread):
    def __init__(self):
        threading.Thread.__init__(self)
    def run(self):
        while 1:
            response, f, args, kwargs = request_queue.get()
            response.put(f(*args, **kwargs))

thread_pool = [WorkerThread() for i in range(2)]
for thread in thread_pool:
    thread.setDaemon(True)
    thread.start()

class MTSearch(object):
    def __init__(self, fps_data, header):
        NUM_READERS = 2
        start = 0
        self.readers = readers = []
        for i in range(1, NUM_READERS):
            end = (i*len(fps_data))//NUM_READERS 
            end = fps_data.index("\n", end)+1
            readers.append(FPSReader(InMemoryFPS(fps_data[start:end]), header))
            start = end
        readers.append(FPSReader(InMemoryFPS(fps_data[start:]), header))

    def reset(self):
        return True
    def fp_iter(self):
        return itertools.chain(*self.readers)
    __iter__ = fp_iter

    def _chemfp_tanimoto_count(self, query, threshold=0.0):
        response = Queue.Queue()
        for reader in self.readers:
            add_task(response, reader._chemfp_tanimoto_count, query, threshold)
        total = 0
        for reader in self.readers:
            total += response.get()
        return total

def test():
    import time
    t1 = time.time()
    # Times [16.6, 8.90, 3.25, 3.30] (no thread release)
    # Times [16.7, 8.93, 3.17, 3.22] (thread release)
    #   OpenMP times were horrible, not recorded
    # Times [0.0, 0.0, 3.88, 3.9]  (loop unroll*4)
    # Times [0.0, 0.0, 3.21, 3.2]  (verify normal)
    # Times [0.0, 0.0, 3.17, 0.0]  (loop unroll*2 - best)
    # Times [0.0, 0.0, 4.23, 4.2]  (unroll*4 a different way)

    #Times [0.0, 0.0, 3.297, 0.0] basic
    #Times [0.0, 0.0, 3.307, 0.0]
    #Times [0.0, 0.0, 3.304, 0.0]
    #Times [0.0, 0.0, 3.327, 0.0]
    #Times [0.0, 0.0, 3.300, 0.0]

    #Times [0.0, 0.0, 3.174, 0.0] unroll*2
    #Times [0.0, 0.0, 3.247, 0.0]
    #Times [0.0, 0.0, 3.198, 0.0]
    #Times [0.0, 0.0, 3.220, 0.0]
    #Times [0.0, 0.0, 3.245, 0.0]

    #Times [0.0, 0.0, 3.161, 0.0] correct unrolling
    #Times [0.0, 0.0, 3.180, 0.0]
    #Times [0.0, 0.0, 3.180, 0.0]
    #Times [0.0, 0.0, 3.197, 0.0]


    reader = open_in_memory("nci.fps")

    # Times [16.9, 9.35, 3.78, 3.86]
    #reader = open("nci.fps")

    # Times [0.0, 0.0, 1.9, 0.0] (if I use 2 threads)
    # Times [0.0, 0.0, 2.1, 0.0] (if I use 3 threads)

    #reader = open_in_memory("nci.fps")
    #reader = MTSearch(reader._fps_io.fps_data, reader.header)

    times = [0.0, 0.0, 0.0, 0.0]
    for i in range(100):
        """
        reader.reset()
        t1 = time.time()
        n = 0
        for fp, id in reader:
            n += 1
        print "iter", n, "last", len(fp), repr(id)
        dt = time.time() - t1
        times[0] += dt

        t1 = time.time()
        reader.reset()
        n = 0
        for id, score in reader._chemfp_tanimoto(FP, 0.0):
            n += 1
            #print id, score
        dt = time.time() - t1
        times[1] += dt
        print "computed", n, "in", dt
"""
        reader.reset()
        t1 = time.time()
        count = reader._chemfp_tanimoto_count(FP, 0.85)
        dt = time.time()-t1
        print "count", count, "in", dt
        times[2] += dt
        reader.reset()

        """
        t1 = time.time()
        x = reader._chemfp_tanimoto_knearest(FP, 10, 0.0)
        dt = time.time() - t1
        print "knearest", len(x), "results in", dt
        times[3] += dt
"""

    print "Times", times
if __name__ == "__main__":
    test()
