import array
import ctypes
import functools
import heapq
import operator

import _chemfp
from _chemfp import (hex_isvalid, hex_popcount, hex_intersect_popcount,
                     hex_tanimoto, hex_contains)

from _chemfp import (byte_popcount, byte_intersect_popcount,
                     byte_tanimoto, byte_contains)

from _chemfp import nlargest_tanimoto_block as _nlargest_tanimoto_block
from _chemfp import hex_nlargest_tanimoto_block as _hex_nlargest_tanimoto_block
from _chemfp import intersect_popcount_count as _intersect_popcount_count

def _score_iterator(calculate_score, target_fps):
    for i, (target_fp, target_id) in enumerate(target_fps):
        #print "Calc for", repr(target_fp), target_fps
        yield calculate_score(target_fp), -i, target_id


def _nlargest_tanimoto(n, calculate_score, target_fps):
    result = heapq.nlargest(n, _score_iterator(calculate_score, target_fps))
    return [(-i, target_id, target_score) for (target_score, i, target_id) in result]

def byte_nlargest_tanimoto(n, query_fp, target_fps):
    return _nlargest_tanimoto(n, functools.partial(byte_tanimoto, query_fp), target_fps)


def byte_nlargest_tanimoto_block(n, query_fp, target_block, offset=0,
                                 storage_len=None, threshold=0.0):
    if storage_len is None:
        storage_len = len(query_fp)
    indicies = array.array('i', (0,)*n)
    scores = array.array('d', (0.0,)*n)
    k = _nlargest_tanimoto_block(n, query_fp, target_block, offset,
                                 storage_len, threshold, indicies, scores)
    return [(indicies[i], scores[i]) for i in xrange(k)]

def byte_intersect_popcount_count(query, arena, min_overlap, offset=0,
                                  storage_len=None):
    if storage_len is None:
        storage_len = len(query)
    response = _intersect_popcount_count(query, arena, offset, storage_len, min_overlap)
    assert response >= 0
    return response

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
    
class FPSTanimotoSearch(object):
    def __init__(self, hex_query, k, threshold=0.0, lineno=1):
        self.hex_query = hex_query
        self.k = k
        self.threshold = threshold

        # Let ctypes worry about memory management
        self._indicies = (ctypes.c_int*k)()
        self._scores = (ctypes.c_double*k)()
        self._id_start_offsets = (ctypes.c_void_p*k)()
        self._id_lens = (ctypes.c_int*k)()

        self._lineno_ptr = (ctypes.c_int)()
        self._lineno_ptr.value = lineno

        self._identifiers = {}
        self._heap = TanimotoHeap()
        _chemfp.fps_heap_init(self._heap, k, threshold, self._indicies, self._scores,
                               self._id_start_offsets, self._id_lens)

    def search_block(self, target_block):
        if self._heap is None:
            raise TypeError("Cannot call 'search_block' after calling 'finish'")
        err = _chemfp.fps_heap_update_tanimoto(self._heap, self.hex_query, target_block)
        if err < 0:
            XXX
        # This is an O(k * N) operation
        # If k depends on N, this becomes O(N**2)
        # The algorithm can be made more efficient for those cases.
        # (Eg, have a max id length, or allow memory allocation)
        id_start_offsets = self._id_start_offsets
        id_lens = self._id_lens
        indicies = self._indicies
        identifiers = self._identifiers
        
        for i in range(self._heap.size):
            # Look for any new entries. These will have a non-zero
            # position since the first id must be at least at offset 3
            # (two hex characters for a fingerprint byte, and then a
            # space and then the id).
            if id_start_offsets[i] > 0:
                identifiers[indicies[i]] = ctypes.string_at(id_start_offsets[i], id_lens[i])

                # Reset the start offset to 0 so that next time around
                # we know it's been seen. The sort doesn't look at this
                # field so this change won't affect the results.
                id_start_offsets[i] = 0

    def finish(self):
        heap = self._heap
        _chemfp.fps_heap_finish_tanimoto(heap)
        self._heap = None

        identifiers = self._identifiers
        indicies = self._indicies
        scores = self._scores

        return [(identifiers[indicies[i]], scores[i]) for i in range(heap.size)]
        
        
class FPSTanimotoCount(object):
    def __init__(self, hex_query, threshold=0.0, lineno=1):
        self.hex_query = hex_query
        self.threshold = threshold

        self._lineno_ptr = (ctypes.c_int)()
        self._lineno_ptr.value = _lineno
        self._count_ptr = (ctypes.c_int)()
        self._count_ptr.value = 0
        

    # Or should I call this 'count_block'?
    def search_block(self, target_block):
        err = _chemfp.fps_tanimoto_count(self.hex_query, target_block,
                                         self.threshold, self._count_ptr, self._lineno_ptr)
        if err < 0:
            raise _error(err, source, self._lineno_ptr.value)

    def get_current_count(self):
        return self._count_ptr.value

## Do I want closer API similarity? But note that FPSTanimotoSearch does not
## allow searching once finish() is called, while FPSTanimotoCount allows it.
#    def finish(self):
#        return self._count_ptr.value
