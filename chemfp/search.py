"""This is an internal module. Do not call these functions directly"""
from __future__ import division

import operator
import heapq
import ctypes

from chemfp import bitops
import _chemfp

# This code is so much simpler than the C extension code.
# And about 6 times slower, even with bitops.byte_tanimoto in C!

def generic_tanimoto_search(query, targets, threshold):
    for fp, id in targets:
        score = bitops.byte_tanimoto(query, fp)
        if score == -1:
            raise AssertionError("should not happen")
        if score >= threshold:
            yield score, id
    
def generic_tanimoto_count_batch(queries, targets, threshold):
    results = []
    for query in queries:
        results.append( sum(1 for (fp, id) in targets if bitops.byte_tanimoto(query, fp) >= threshold) )
    return results
                        
def generic_tanimoto_knearest_search(queries, targets, k, threshold):
    results = []
    for query in queries:
        results.append(heapq.nlargest(k, generic_tanimoto_search(query, targets, threshold),
                                      key=operator.itemgetter(1)))
    return results

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


def block_tanimoto_knearest_search_batch(queries, targets, k, threshold):
    if k == 0:
        return [[] for q in queries]

    lineno = targets._first_fp_lineno # XXX hide?
    searchers = [BlockTanimotoSearcher(query, k, threshold, lineno) for query in queries]

    for block in targets.iter_blocks():
        for searcher in searchers:
            searcher.feed(block)

    return [searcher.finish() for searcher in searchers]

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

def block_tanimoto_count_batch(queries, targets, threshold):
    lineno = targets._first_fp_lineno
    counters = [BlockTanimotoCounter(query, threshold, lineno) for query in queries]

    for block in targets.iter_blocks():
        for counter in counters:
            counter.feed(block)
    return [counter.finish() for counter in counters]

# T = A&B / A|B
# T = number of bits in common / (number A + number B - number of bits in common)
# T = C / (A+B-C)
# T*(A+B-C) = C
# T*(A+B)-T*C = C
# T*(A+B) = C+T*C = C(1+T)
# C = T*(A+B) / (1+T)

def arena_tanimoto_count_batch(queries, arenas, threshold):
    if threshold <= 0.0:
        # Everything is close enough
        n = sum(len(arena) // len(query) for arena in arenas)
        return [n for query in queries]

    results = []
    for query in queries:
        query_popcount = bitops.byte_popcount(query)
        if query_popcount == 0:
            # By definition, everything has a similarity of 0.0 to this query
            # (including other cases with no bits set). Since I handled
            # the threshold <= 0.0 case earlier, the answer here must be 0.
            results.append(0)
            continue
            
        count = 0

        # XXX I should have something which rounds up
        if threshold == 0.0:
            range_ = xrange(0, len(arenas))
        else:
            range_ = xrange(int(query_popcount*threshold), int(query_popcount/threshold)+1)

        for target_popcount in range_:
            in_common = threshold * (query_popcount+target_popcount)/(1+threshold)
            min_equal_bits = int(in_common)
            if min_equal_bits != in_common:
                min_equal_bits += 1
            count += bitops.byte_intersect_popcount_count(query, arenas[target_popcount], min_equal_bits)

        results.append(count)
    return results


# queries iterates over fingerprints
# arena is a block of target fingerprints as a byte string
# storage_len is the fingerprint size in bytes (may include trailing 0 padding)
# popcount_offsets contains index information for each popcount in the arena
def x_arena_tanimoto_knearest_search_batch(queries, arena, storage_len,
                                         popcount_offsets, k, threshold):
    assert k >= 1
    assert 0.0 <= threshold <= 1.0
    results = []

    # Not yet handling padded zeros. Don't know if I should have popcount_offsets
    # include the 0s for the pad characters.
    max_popcount = len(popcount_offsets) - 1
    end = popcount_offsets[-1]
    if max_popcount != (storage_len*8 + 1):
        raise TypeError("Size mismatch between max popcount and storage_len")

    hit_indicies = (ctypes.c_int*k)()
    hit_scores = (ctypes.c_double*k)()
    
    for query in queries:
        num_hits = _chemfp.arena_tanimoto_knearest_search(query,
                                                          arena, storage_len, 0, end,
                                                          popcount_offsets, k, threshold,
                                                          hit_indicies, hit_scores)
        assert num_hits >= 0, num_hits
        
        # Copy over the hits 
        results.append( (hit_indicies[:hits], hit_scores[:hits]) )
        
    return results


def threshold_tanimoto_arena_arena(threshold, queries, targets):
    assert queries.header.num_bits == targets.header.num_bits
    num_bits = queries.header.num_bits

    num_queries = len(queries)

    offsets = (ctypes.c_int * (num_queries+1))()
    offsets[0] = 0

    num_cells = min(100, len(queries)) * len(targets)
    indicies = (ctypes.c_int * num_cells)()
    scores = (ctypes.c_double * num_cells)()
    
    query_start = 0
    query_end = len(queries)


    def add_rows(query_start):
        return _chemfp.threshold_tanimoto_arena(
            threshold, num_bits,
            queries.storage_size, queries.arena, query_start, query_end,
            targets.storage_size, targets.arena, 0, -1,
            targets.popcount_indicies,
            offsets, query_start,
            indicies, scores)

    return _search(query_end, offsets, indicies, scores, add_rows)

def _search(query_end, offsets, indicies, scores, add_rows):
    num_added = add_rows(0)
    if num_added == query_end:
        return SearchResults(offsets, indicies, scores)

    query_start = num_added

    last = offsets[num_added]
    all_indicies = indicies[:last]
    all_scores = scores[:last]

    while query_start < query_end:
        num_added = add_rows(query_start)
        assert num_added > 0

        prev_last = offsets[query_start]
        all_indicies[prev_last:] = indicies
        all_scores[prev_last:] = scores

        query_start += num_added

    return SearchResults(offsets, all_indicies, all_scores)


def knearest_tanimoto_arena_arena(k, threshold, queries, targets):
    assert queries.header.num_bits == targets.header.num_bits
    num_bits = queries.header.num_bits

    num_queries = len(queries)

    offsets = (ctypes.c_int * (num_queries+1))()
    offsets[0] = 0

    num_cells = min(100, len(queries))*k

    indicies = (ctypes.c_int * num_cells)()
    scores = (ctypes.c_double * num_cells)()

    query_start = 0
    query_end = len(queries)

    def add_rows(query_start):
        return _chemfp.klargest_tanimoto_arena(
            k, threshold, num_bits,
            queries.storage_size, queries.arena, query_start, query_end,
            targets.storage_size, targets.arena, 0, -1,
            targets.popcount_indicies,
            offsets, query_start,
            indicies, scores)

    return _search(query_end, offsets, indicies, scores, add_rows)


class SearchResults(object):
    def __init__(self, offsets, indicies, scores):
        assert len(offsets) > 0
        if offsets:
            assert indicies[offsets[-2]]
            assert scores[offsets[-2]]
        self.offsets = offsets
        self.indicies = indicies
        self.scores = scores
    def __len__(self):
        return len(offsets)-1

    def size(self, i):
        return self.offsets[i+1]-self.offsets[i]
    
    def __getitem__(self, i):
        start, end = self.offsets[i:i+2]
        return zip(self.indicies[start:end], self.scores[start:end])

    def __iter__(self):
        start = self.offsets[0]
        for end in self.offsets[1:]:
            yield zip(self.indicies[start:end], self.scores[start:end])
            start = end
