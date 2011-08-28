"""This is an internal module. Do not call these functions directly"""
from __future__ import division

import operator
import heapq
import ctypes
import array
import itertools

from chemfp import bitops
import _chemfp

# This code is so much simpler than the C extension code.
# And about 6 times slower, even with bitops.byte_tanimoto in C!


#
# Support code for fast FPS searches
# 


"""
class BlockTanimotoCounter(object):
    __slots__ = ["hex_query", "threshold", "lineno_ptr", "count_ptr"]
    def __init__(self, query, threshold, lineno):
        self.hex_query = query.encode("hex")
        self.threshold = threshold
        self.lineno = lineno
        self.count_ptr = (ctypes.c_int)()
        self.count_ptr.value = 0

    def feed(self, block):
        err = _chemfp.fps_tanimoto_count(self.hex_query, block, self.threshold,
                                         self.count_ptr, self.lineno_ptr)
        if err < 0:
            raise _error(err, self.lineno_ptr.value)

    def finish(self):
        return self.count_ptr.value

def block_tanimoto_count(query_fp, targets, threshold):
    result = block_tanimoto_count_all( [(None, query_fp)], targets, threshold )
    return result[0][1]

def block_tanimoto_count_all(queries, targets, threshold):
    lineno = targets._first_fp_lineno
    query_ids = []
    counters = []
    for query_id, query_fp in queries:
        query_ids.append(query_id)
        counters.append(BlockTanimotoCounter(query_id, threshold, lineno))

    for block in targets.iter_blocks():
        for counter in counters:
            counter.feed(block)
    return zip(query_ids, (counter.finish() for counter in counters))

"""

##############

def _make_knearest_search(num_queries, k):
    class TanimotoHeap(ctypes.Structure):
        _fields_ = [("size", ctypes.c_int),
                    ("heap_state", ctypes.c_int),
                    ("indicies", POINTER(ctypes.c_int*k)),
                    ("ids", POINTER(ctypes.c_int*k)),
                    ("scores", POINTER(ctypes.c_double*k))]

    class KNearestSearch(ctypes.Structure):
        _fields_ = [("num_queries", ctypes.c_int),
                    ("fp_size", ctypes.c_int),
                    ("k", ctypes.c_int),
                    ("search_state", ctypes.c_int),
                    ("threshold", ctypes.c_double),
                    ("queries", ctypes.c_str),
                    ("heaps", POINTER(TanimotoHeap*num_queries)),
                    ("num_targets_processed", ctypes.c_int),
                    ("_all_ids", ctypes.c_void_ptr),
                    ("_all_scores", ctypes.c_void_ptr)]

    return KNearestSearch()


def knearest_tanimoto_search(queries, targets, k, threshold):
    id, fp = queries[0]
    fp_size = len(fp)
    query_str = "".join(fp for (id, fp) in queries)

    blah = _make_knearest_search(len(queries), k)
    
    err = _chemfp.fps_knearest_search_alloc(blah, fp_size, query_str, k, threshold)
    if err < 0:
        error
    assert blah.fp_size == fp_size
    try:
        for block in targets.iter_blocks():
            _chemfp.fps_knearest_search_feed(blah, block)

        _chemfp.fps_knearest_finish(blah)

        results = []
        for heap in blah.heaps:
            results.append( [(heap.ids[i], heap.scores[i]) for i in xrange(heap.size)] )
        return results
    finally:
            fps_knearest_search_free(blah)


class TanimotoCell(ctypes.Structure):
    _fields_ = [("score", ctypes.c_double),
                ("query_index", ctypes.c_int),
                ("id_start", ctypes.c_int),
                ("id_end", ctypes.c_int)]


def threshold_tanimoto_search_fp(query_fp, target_reader, threshold):
    assert query_arena.header.num_bits == target_reader.header.num_bits
    hits = []

    fp_size = len(query_fp)
    if target_reader.num_bytes_per_fp != fp_size:
        ERROR
    num_bits = fp_size * 8
        
    NUM_CELLS = 1000
    cells = (TanimotoCell*num_cells)()

    lineno = target_reader._first_fp_lineno
    
    for block in targets.iter_blocks():
        start = 0
        while 1:
            start, num_lines, num_cells = _chemfp.fps_threshold_tanimoto_search(
                num_bits, fp_size, query_fp, 0, -1,
                block, start, -1,
                threshold, cells)
            if num_cells < 0:
                raise TypeError("Error at or after line %r of %r: %s" %
                                (lineno, target_reader, _chemfp.strerror(num_cells)))
                
            for cell in itertools.islice(cells, 0, num_cells):
                    id = block[cell.id_start:cell.id_end]
                    hits.append( (id, cell.score) )
            lineno += num_lines
            
            if start == end:
                break
    return hits

def threshold_tanimoto_search_all(queries, target_reader, threshold):
    if not queries:
        return []
    all_hits = [[] for i in xrange(len(queries))]

    total_num_lines = 0

    
    # Compute at least 100 tanimotos per query, but at most 10,000 at a time
    # (That's about 200K of memory)
    NUM_CELLS = max(10000, len(queries) * 100)
    cells = (TanimotoCell*NUM_CELLS)()

    lineno = array.array("i", (0,))

    for block in target_reader.iter_blocks():
        start = 0
        end = len(block)
        while 1:
            print start, end
            start, num_lines, num_cells = _chemfp.fps_threshold_tanimoto_search(
                queries.header.num_bits, queries.storage_size, queries.arena, 0, -1,
                block, start, end,
                threshold, cells)
            print " =>", num_cells, start, end
            if num_cells < 0:
                ERROR
                
            for cell in itertools.islice(cells, 0, num_cells):
                id = block[cell.id_start:cell.id_end]
                all_hits[cell.query_index].append( (id, cell.score) )

            if start == end:
                break

    return all_hits
            
# queries must be a list
def threshold_tanimoto_search_all(queries, targets, threshold):
    if not queries:
        return []
    query_ids, query_fps = zip(*queries)
    all_hits = threshold_tanimoto_search_fps(query_fps)
    return zip(query_ids, all_hits)
        
    
    
# This does not make sense as targets is consumed
#def threshold_tanimoto_search(queries, targets, threshold):
    
    
    


##############

"""

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


class BlockTanimotoSearcher(object):
    __slots__ = ["hex_query", "indicies", "scores", "id_starts", "id_lens", "lineno_ptr", "identifiers", "heap"]
    def __init__(self, query, k, threshold, lineno):
        self.hex_query = query.encode("hex")
        self.indicies = (ctypes.c_int*k)()
        self.scores = (ctypes.c_double*k)()
        self.id_starts = (ctypes.c_void_p*k)()
        self.id_lens = (ctypes.c_int*k)()

        self.lineno_ptr = (ctypes.c_int)()
        self.lineno_ptr.value = lineno

        self.identifiers = {}
        # Use ctypes to alloc raw space
        # XXX add class methods?
        self.heap = TanimotoHeap()
        _chemfp.fps_heap_init(self.heap, k, threshold,
                              self.indicies, self.scores,
                              self.id_starts, self.id_lens)
    def feed(self, target_block):
        err = _chemfp.fps_heap_update_tanimoto(self.heap, self.hex_query, target_block)
        if err < 0:
            # XXX lineno_ptr.value ?
            raise _error(err, 0)

        id_starts = self.id_starts
        id_lens  = self.id_lens
        indicies = self.indicies
        for i in range(self.heap.size):
            if id_starts[i] > 0:
                self.identifiers[indicies[i]] = ctypes.string_at(id_starts[i], id_lens[i])
                id_starts[i] = 0

    def finish(self):
        _chemfp.fps_heap_finish_tanimoto(self.heap)
        return [(self.identifiers[self.indicies[i]], self.scores[i]) for i in range(self.heap.size)]


def block_knearest_tanimoto_search_fp(query_fp, targets, k, threshold):
    results = block_knearest_tanimoto_search([(None, query_fp)], targets, k, threshold)
    return results[0][1]

def block_knearest_tanimoto_search_all(queries, targets, k, threshold):
    if k == 0:
        return [(query_id, []) for (query_id, query_fp) in queries]

    lineno = targets._first_fp_lineno # XXX hide?
    query_ids = []
    searchers = []
    for query_id, query_fp in queries:
        query_ids.append(query_id)
        searchers.append(BlockTanimotoSearcher(query_fp, k, threshold, lineno))

    for block in targets.iter_blocks():
        for searcher in searchers:
            searcher.feed(block)

    return zip(query_ids, (searcher.finish() for searcher in searchers))

# This makes no sense, because targets is a consumable iterator.
# def block_knearest_tanimoto_search_iter(): raise AssertionError

####

"""
