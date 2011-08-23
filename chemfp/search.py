"""This is an internal module. Do not call these functions directly"""
from __future__ import division

import operator
import heapq
import ctypes

from chemfp import bitops
import _chemfp

# This code is so much simpler than the C extension code.
# And about 6 times slower, even with bitops.byte_tanimoto in C!



#
# Support code for fast FPS searches
# 

# XXX What is this?
def _error(errno, lineno):
    raise AssertionError( (errno, _chemfp.strerror(errno), lineno) )

import ctypes

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



class BlockTanimotoCounter(object):
    __slots__ = ["hex_query", "threshold", "lineno_ptr", "count_ptr"]
    def __init__(self, query, threshold, lineno):
        self.hex_query = query.encode("hex")
        self.threshold = threshold
        self.lineno_ptr = (ctypes.c_int)()
        self.lineno_ptr.value = lineno
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



##############
