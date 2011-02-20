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
            raise _error(err, self._f, lineno_ptr.value)

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
            raise _error(err, self._f, self.lineno_ptr.value)

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



def arena_tanimoto_knearest_search_batch(queries, arenas, arena_ids, k, threshold):
    results = []
    num_arenas = len(arenas)
    for query in queries:
        query_popcount = bitops.byte_popcount(query)
        if query_popcount == 0:
            # By definition this will be 0.
            results.append([])
            continue

        # This is the same algorithm in Swamidass and Baldi.  I tried
        # the obvious merge sort but it was slower.
        # (Probably because of Python function call overhead and
        # because the built-in timsort handles this case well.)
        scores = []
        for target_popcount in range(0, len(arenas)):
            if target_popcount < query_popcount:
                best_possible_score = target_popcount/query_popcount
            else:
                best_possible_score = query_popcount/target_popcount
            if best_possible_score >= threshold:
                scores.append( (best_possible_score, target_popcount) )
        scores.sort(reverse=True)

        best = []
        best_len = 0
        for negative_best_possible_score, target_popcount in scores:
            best_possible_score = negative_best_possible_score

            if best_len == k and best[0][0] > best_possible_score:
                break

            arena = arenas[target_popcount]
            target_ids = arena_ids[target_popcount]
            for (index, score) in bitops.byte_nlargest_tanimoto_block(k, query, arena, 0, threshold=threshold):
                if best_len == k:
                    if best[0][0] < score:
                        new_hit = (score, target_ids[index])
                        heapq.heapreplace(best, new_hit)
                        if best[0][0] >= best_possible_score:
                            break
                else:
                    new_hit = (score, target_ids[index])
                    best_len += 1
                    heapq.heappush(best, new_hit)

        # XXX Should I use a bunch of heap pops here?
        # Probably end up with Python function call overhead.
        results.append( [(id, score) for (score, id) in sorted(best, reverse=True)] )
    return results

